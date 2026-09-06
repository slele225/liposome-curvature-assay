function run_matlab_reference(condDir, chNames, varargin)
% run_matlab_reference(condDir, chNames, 'Seed', 1, 'Markers', {...})
%
% Runs the UNMODIFIED cmeAnalysis workflow
%     data = loadConditionData(condDir, chNames, markers);
%     rng(seed);
%     runDetection(data, 'Overwrite', true);
% with timing, then exports every detection_v2.mat as a TSV
% (<master>/Detection/detection_matlab.tsv) with the same columns as the C++
% tool's detection_cpp.tsv for end-to-end comparison.

ip = inputParser;
ip.addParameter('Seed', 1);
ip.addParameter('Markers', []);
ip.parse(varargin{:});

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, '..', '..', 'cmeAnalysis-master', 'software')));

markers = ip.Results.Markers;
if isempty(markers), markers = repmat({'cy5'}, 1, numel(chNames)); end

t0 = tic;
data = loadConditionData(condDir, chNames, markers);
tLoad = toc(t0);
t1 = tic;
rng(ip.Results.Seed);
runDetection(data, 'Overwrite', true);
tDet = toc(t1);
fprintf('TIMING loadConditionData %.3f s, runDetection %.3f s, total %.3f s\n', tLoad, tDet, tLoad + tDet);

nCh = numel(chNames);
dn = {'x', 'y', 'A', 'c', 'x_pstd', 'y_pstd', 'A_pstd', 'c_pstd', 'sigma_r', 'SE_sigma_r', 'RSS', 'pval_Ar', 'hval_Ar', 'hval_AD', 'isPSF', 's', 'dRange_min', 'dRange_max'};
for i = 1:numel(data)
    matPath = [data(i).channels{1} 'Detection' filesep 'detection_v2.mat'];
    S = load(matPath);
    fi = S.frameInfo;
    out = [data(i).channels{1} 'Detection' filesep 'detection_matlab.tsv'];
    fid = fopen(out, 'w');
    fprintf(fid, 'movie\tframe\tindex\tsource_image');
    for c = 1:nCh
        for f = 1:numel(dn), fprintf(fid, '\t%s_%s', dn{f}, chNames{c}); end
    end
    fprintf(fid, '\tx_init\ty_init\tmaskA\tmaskN\tmask_Ar\n');
    cp = data(i).channels{1};
    idxs = regexp(cp, filesep);
    cellName = cp(idxs(end-2)+1:idxs(end-1)-1);
    for k = 1:numel(fi)
        F = fi(k);
        if iscell(data(i).framePaths{1}), src = data(i).framePaths{1}{k}; else, src = data(i).framePaths{1}; end
        np = size(F.x, 2);
        for p = 1:np
            fprintf(fid, '%s\t%d\t%d\t%s', cellName, k, p, src);
            for c = 1:nCh
                fprintf(fid, '\t%.17g', F.x(c,p), F.y(c,p), F.A(c,p), F.c(c,p), F.x_pstd(c,p), F.y_pstd(c,p), F.A_pstd(c,p), F.c_pstd(c,p), ...
                    F.sigma_r(c,p), F.SE_sigma_r(c,p), F.RSS(c,p), F.pval_Ar(c,p));
                fprintf(fid, '\t%d\t%d\t%d', F.hval_Ar(c,p), F.hval_AD(c,p), F.isPSF(c,p));
                fprintf(fid, '\t%.17g\t%.17g\t%.17g', F.s(c), F.dRange{c}(1), F.dRange{c}(2));
            end
            fprintf(fid, '\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n', F.x_init(p), F.y_init(p), F.maskA(p), F.maskN(p), F.mask_Ar(p));
        end
    end
    fclose(fid);
    fprintf('exported %s\n', out);
end
end
