function dump_reference(condDir, chNames, outDir, varargin)
% dump_reference(condDir, chNames, outDir, 'NumImageMovies', n, 'Seed', 1)
%
% Reproduces  data = loadConditionData(condDir, chNames, markers); rng(seed); runDetection(data)
% step by step with the ORIGINAL cmeAnalysis functions and writes every
% intermediate in the same layout as `cme_detect --dump-dir`, so that
% tests/compare_reference.py can diff the two.  The only modified copy of a
% cmeAnalysis function is psd_dump.m (pointSourceDetection with dump hooks);
% the final results are additionally cross-checked by run_matlab_reference.m,
% which calls the untouched runDetection.
%
% Temporary reference tooling for the C++ port (cpp_port/); GPL-3.0 like cmeAnalysis.

ip = inputParser;
ip.addParameter('NumImageMovies', 2);
ip.addParameter('Seed', 1);
ip.addParameter('Markers', []);
ip.parse(varargin{:});
nImgMovies = ip.Results.NumImageMovies;

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, '..', '..', 'cmeAnalysis-master', 'software')));
addpath(here);

markers = ip.Results.Markers;
if isempty(markers), markers = repmat({'cy5'}, 1, numel(chNames)); end
data = loadConditionData(condDir, chNames, markers);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rng(ip.Results.Seed);

nd = numel(data);
nCh = numel(data(1).channels);
mCh = 1;

%% ---------------- sigma estimation (runDetection + getGaussianPSFsigmaFromData) ----------------
sigmaRaw = zeros(1, nCh);
nf = round(40/nd);
for c = 1:nCh
    frames = cell(nd, nf);
    for i = 1:nd
        fidx = round(linspace(1, data(i).movieLength, nf));
        if iscell(data(i).framePaths{c})
            frames(i,:) = arrayfun(@(f) double(imread(data(i).framePaths{c}{f})), fidx, 'unif', 0);
        else
            frames(i,:) = arrayfun(@(f) double(readtiff(data(i).framePaths{c}, f)), fidx', 'unif', 0);
        end
    end
    imageList = vertcat(frames(:));
    n = numel(imageList);
    svect = cell(1, n);
    for i = 1:n
        img = double(imageList{i});
        pstruct = pointSourceDetection(img, 1.5, 'Mode', 'xyac');
        if ~isempty(pstruct)
            pstruct = fitGaussians2D(img, pstruct.x, pstruct.y, pstruct.A, 1.5*ones(1,length(pstruct.x)), pstruct.c, 'xyasc');
            isPSF = ~[pstruct.hval_AD] & [pstruct.pval_Ar] < 0.05;
            svect{i} = pstruct.s(~isnan(pstruct.s) & isPSF);
            write_pstruct(fullfile(outDir, sprintf('sigma_ch%d_refit_img%03d.tsv', c, i)), pstruct);
        end
    end
    svect_cells = svect;
    svect = [svect{:}];
    write_col(fullfile(outDir, sprintf('sigma_ch%d_svect.tsv', c)), 's', svect);
    cnt = cellfun(@numel, svect_cells);
    fid = fopen(fullfile(outDir, sprintf('sigma_ch%d_svect_per_image.tsv', c)), 'w');
    fprintf(fid, 'image\tcount\n');
    fprintf(fid, '%d\t%d\n', [1:n; cnt]);
    fclose(fid);

    opts = statset('maxIter', 200);
    fid = fopen(fullfile(outDir, sprintf('sigma_ch%d_gmm.tsv', c)), 'w');
    fprintf(fid, 'k\tcomponent\tmu\tSigma\tPComponents\tNlogL\tBIC\titers\tconverged\n');
    gmmFailed = 0; chosenK = 0; chosenComponent = 0;
    try
        w = warning('off', 'stats:gmdistribution:FailedToConverge');
        obj = cell(1,3);
        for k = 1:3
            obj{k} = gmdistribution.fit(svect', k, 'Options', opts);
            for j = 1:k
                fprintf(fid, '%d\t%d\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%d\t%d\n', k, j, obj{k}.mu(j), obj{k}.Sigma(1,1,j), ...
                    obj{k}.PComponents(j), obj{k}.NlogL, obj{k}.BIC, obj{k}.Iters, obj{k}.Converged);
            end
        end
        [~,idx] = min(cellfun(@(i) i.BIC, obj));
        chosenK = idx;
        obj = obj{idx};
        [mu,idx] = sort(obj.mu);
        svec = sqrt(squeeze(obj.Sigma(:,:,idx)));
        amp = obj.PComponents(idx)';
        [~,idx] = max(amp./(sqrt(2*pi)*svec));
        chosenComponent = idx;
        sigmaRaw(c) = mu(idx);
        warning(w);
    catch
        gmmFailed = 1;
        sigmaRaw(c) = mean(svect);
    end
    fprintf(fid, '# gmmFailed=%d chosenK=%d chosenComponent=%d\n', gmmFailed, chosenK, chosenComponent);
    fclose(fid);
end
sigma = sigmaRaw;
sigma(sigma<1.1) = 1.1;
fid = fopen(fullfile(outDir, 'sigma.tsv'), 'w');
fprintf(fid, 'channel\tsigma_raw\tsigma\n');
for c = 1:nCh, fprintf(fid, '%d\t%.17g\t%.17g\n', c, sigmaRaw(c), sigma(c)); end
fclose(fid);
fprintf('sigma = %s\n', mat2str(sigma, 17));

%% ---------------- detection per movie / frame (runDetection main) ----------------
dfields = {'x', 'y', 'A', 'c', 'x_pstd', 'y_pstd', 'A_pstd', 'c_pstd', 'sigma_r', 'SE_sigma_r', 'RSS', 'pval_Ar'};
lfields = {'hval_Ar', 'hval_AD', 'isPSF'};
sfields = [dfields {'hval_Ar', 'hval_AD'}];
rmfields = [dfields lfields {'x_init', 'y_init', 'maskA', 'maskN', 'mask_Ar'}];

for i = 1:nd
    d = data(i);
    cp = d.channels{1};
    idxs = regexp(cp, filesep);
    cellName = cp(idxs(end-2)+1:idxs(end-1)-1);
    mdir = fullfile(outDir, sprintf('movie%03d_%s', i, cellName));
    if ~exist(mdir, 'dir'), mkdir(mdir); end
    writeImages = i <= nImgMovies;
    for k = 1:d.movieLength
        if ~iscell(d.framePaths{mCh})
            img = double(readtiff(d.framePaths{mCh}, k));
        else
            img = double(imread(d.framePaths{mCh}{k}));
        end
        pre = fullfile(mdir, sprintf('frame%04d_', k));
        [pstruct, mask, dbg] = psd_dump(img, sigma(mCh), 'Alpha', 0.05, 'Mask', [], 'RemoveRedundant', true);
        if writeImages
            write_bin([pre 'imgLoG.bin'], dbg.imgLoG);
            write_bin([pre 'A_est.bin'], dbg.A_est);
            write_bin([pre 'c_est.bin'], dbg.c_est);
            write_bin([pre 'pval_prefilter.bin'], dbg.pval);
            write_bin([pre 'mask_prefilter.bin'], uint8(dbg.maskPrefilter));
            write_bin([pre 'mask_combined.bin'], uint8(dbg.maskCombined));
            write_bin([pre 'mask_final.bin'], uint8(mask));
        end
        fid = fopen([pre 'lm.tsv'], 'w'); fprintf(fid, 'lmx\tlmy\n'); fprintf(fid, '%.17g\t%.17g\n', [dbg.lmx(:) dbg.lmy(:)]'); fclose(fid);
        fid = fopen([pre 'logThreshold.txt'], 'w'); fprintf(fid, '%.17g\n', dbg.logThreshold); fclose(fid);
        write_pstruct([pre 'master_fitall.tsv'], dbg.fitAll);
        write_col([pre 'master_keep.tsv'], 'keepNonNaN', dbg.keepNonNaN);
        write_col([pre 'master_keepfinal.tsv'], 'keepFinal', dbg.keepFinal);
        if ~isempty(pstruct)
            write_pstruct([pre 'master_final.tsv'], pstruct);
        else
            write_pstruct([pre 'master_final.tsv'], []);
        end

        dRange = cell(1, nCh);
        dRange{mCh} = [min(img(:)) max(img(:))];
        if ~isempty(pstruct)
            pstruct.s = sigma;
            pstruct = rmfield(pstruct, 's_pstd');
            np = numel(pstruct.x);
            for f = 1:length(dfields)
                tmp = NaN(nCh, np); tmp(mCh,:) = pstruct.(dfields{f}); pstruct.(dfields{f}) = tmp;
            end
            for f = 1:length(lfields)
                tmp = false(nCh, np); tmp(mCh,:) = pstruct.(lfields{f}); pstruct.(lfields{f}) = tmp;
            end
            CC = bwconncomp(mask);
            labels = labelmatrix(CC);
            loclabels = labels(sub2ind(size(img), pstruct.y_init, pstruct.x_init));
            compSize = cellfun(@(i) numel(i), CC.PixelIdxList);
            pstruct.maskN = compSize(loclabels);
            compInt = cellfun(@(i) sum(img(i))/numel(i), CC.PixelIdxList);
            pstruct.maskA = compInt(loclabels);

            for ci = setdiff(1:nCh, mCh)
                if ~iscell(d.framePaths{ci})
                    simg = double(readtiff(d.framePaths{ci}, k));
                else
                    simg = double(imread(d.framePaths{ci}{k}));
                end
                dRange{ci} = [min(simg(:)) max(simg(:))];
                pstructSlave = fitGaussians2D(simg, pstruct.x(mCh,:), pstruct.y(mCh,:), [], sigma(ci)*ones(1,np), [], 'Ac');
                pstructSlaveLoc = fitGaussians2D(simg, pstruct.x(mCh,:), pstruct.y(mCh,:), pstructSlave.A, sigma(ci)*ones(1,np), pstructSlave.c, 'xyAc');
                idx = sqrt((pstruct.x(mCh,:)-pstructSlaveLoc.x).^2 + (pstruct.y(mCh,:)-pstructSlaveLoc.y).^2) < 3*sigma(mCh) & pstructSlaveLoc.A > pstructSlave.A;
                write_pstruct([pre sprintf('slave_ch%d_fixed.tsv', ci)], pstructSlave);
                write_pstruct([pre sprintf('slave_ch%d_loc.tsv', ci)], pstructSlaveLoc);
                write_col([pre sprintf('slave_ch%d_useloc.tsv', ci)], 'useLoc', idx);
                for f = 1:length(sfields)
                    pstruct.(sfields{f})(ci,~idx) = pstructSlave.(sfields{f})(~idx);
                    pstruct.(sfields{f})(ci,idx) = pstructSlaveLoc.(sfields{f})(idx);
                end
                nanIdx = isnan(pstructSlave.x);
                for f = 1:length(rmfields)
                    pstruct.(rmfields{f})(:,nanIdx) = [];
                end
                np = size(pstruct.x,2);
                pstruct.isPSF(ci,:) = ~pstruct.hval_AD(ci,:);
            end
            write_frameinfo([pre 'frameinfo.tsv'], pstruct, nCh);
        else
            for ci = setdiff(1:nCh, mCh)
                if ~iscell(d.framePaths{ci})
                    simg = double(readtiff(d.framePaths{ci}, k));
                else
                    simg = double(imread(d.framePaths{ci}{k}));
                end
                dRange{ci} = [min(simg(:)) max(simg(:))];
            end
            write_frameinfo([pre 'frameinfo.tsv'], [], nCh);
        end
        fid = fopen([pre 'dRange.tsv'], 'w'); fprintf(fid, 'channel\tmin\tmax\n');
        for c = 1:nCh, fprintf(fid, '%d\t%.17g\t%.17g\n', c, dRange{c}(1), dRange{c}(2)); end
        fclose(fid);
    end
    fprintf('dumped movie %d/%d\n', i, nd);
end
end

function write_bin(path, img)
fid = fopen(path, 'w');
if isa(img, 'uint8')
    fwrite(fid, img, 'uint8'); cls = 'uint8';
else
    fwrite(fid, double(img), 'double'); cls = 'float64';
end
fclose(fid);
fid = fopen([path '.hdr'], 'w'); fprintf(fid, '%s %d %d column-major\n', cls, size(img,1), size(img,2)); fclose(fid);
end

function write_col(path, name, v)
fid = fopen(path, 'w');
fprintf(fid, '%s\n', name);
if islogical(v), fprintf(fid, '%d\n', v); else, fprintf(fid, '%.17g\n', v); end
fclose(fid);
end

function write_pstruct(path, P)
fid = fopen(path, 'w');
fprintf(fid, 'x\ty\tA\ts\tc\tx_pstd\ty_pstd\tA_pstd\ts_pstd\tc_pstd\tx_init\ty_init\tsigma_r\tSE_sigma_r\tRSS\tpval_Ar\tmask_Ar\thval_Ar\thval_AD\n');
if ~isempty(P)
    M = [P.x(:) P.y(:) P.A(:) P.s(:) P.c(:) P.x_pstd(:) P.y_pstd(:) P.A_pstd(:) P.s_pstd(:) P.c_pstd(:) ...
        P.x_init(:) P.y_init(:) P.sigma_r(:) P.SE_sigma_r(:) P.RSS(:) P.pval_Ar(:) P.mask_Ar(:) double(P.hval_Ar(:)) double(P.hval_AD(:))];
    fprintf(fid, [repmat('%.17g\t', 1, 17) '%d\t%d\n'], M');
end
fclose(fid);
end

function write_frameinfo(path, P, nCh)
dn = {'x', 'y', 'A', 'c', 'x_pstd', 'y_pstd', 'A_pstd', 'c_pstd', 'sigma_r', 'SE_sigma_r', 'RSS', 'pval_Ar'};
ln = {'hval_Ar', 'hval_AD', 'isPSF'};
fid = fopen(path, 'w');
fprintf(fid, 'idx');
for c = 1:nCh
    for f = 1:numel(dn), fprintf(fid, '\t%s_%d', dn{f}, c); end
    for f = 1:numel(ln), fprintf(fid, '\t%s_%d', ln{f}, c); end
end
fprintf(fid, '\tx_init\ty_init\tmaskA\tmaskN\tmask_Ar\n');
if ~isempty(P)
    np = size(P.x, 2);
    for p = 1:np
        fprintf(fid, '%d', p);
        for c = 1:nCh
            for f = 1:numel(dn), fprintf(fid, '\t%.17g', P.(dn{f})(c,p)); end
            for f = 1:numel(ln), fprintf(fid, '\t%d', P.(ln{f})(c,p)); end
        end
        fprintf(fid, '\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n', P.x_init(p), P.y_init(p), P.maskA(p), P.maskN(p), P.mask_Ar(p));
    end
end
fclose(fid);
end
