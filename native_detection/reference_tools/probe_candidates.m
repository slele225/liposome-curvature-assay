function probe_candidates(imgPath, dumpPrefix, sigma, candIdx, mode, outfile)
% Re-run the fitGaussian2D MEX on selected candidates of a dumped frame
% (dumpPrefix = '<dump dir>/frame0001_'), reconstructing the window exactly as
% fitGaussians2D.m does (label mask -> NaN), and write windows + init + MEX
% result in the 'C' line format understood by tests/fit_probe.cpp.
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'cmeAnalysis-master', 'software')));
img = double(imread(imgPath));
[ny, nx] = size(img);
A_est = read_bin([dumpPrefix 'A_est.bin'], ny, nx, 'double');
c_est = read_bin([dumpPrefix 'c_est.bin'], ny, nx, 'double');
mask = read_bin([dumpPrefix 'mask_combined.bin'], ny, nx, 'uint8') ~= 0;
lm = readmatrix([dumpPrefix 'lm.tsv'], 'FileType', 'text', 'Delimiter', '\t');
lmx = lm(:,1); lmy = lm(:,2);
labels = bwlabel(mask);
w4 = ceil(4*sigma);
fid = fopen(outfile, 'w');
for q = 1:numel(candIdx)
    p = candIdx(q);
    xi = lmx(p); yi = lmy(p);
    if ~(xi>w4 && xi<=nx-w4 && yi>w4 && yi<=ny-w4), continue; end
    maskWindow = labels(yi-w4:yi+w4, xi-w4:xi+w4);
    maskWindow(maskWindow==maskWindow(w4+1,w4+1)) = 0;
    window = img(yi-w4:yi+w4, xi-w4:xi+w4);
    window(maskWindow~=0) = NaN;
    init = [0 0 A_est(yi, xi) sigma c_est(yi, xi)];
    prm = fitGaussian2D(window, init, mode);
    fprintf(fid, 'C %d init %.17g %.17g %.17g %.17g %.17g prm %.17g %.17g %.17g %.17g %.17g\n', p, init, prm);
    fprintf(fid, '%.17g ', window(:)); fprintf(fid, '\n');
    fprintf('cand %d: prm = %s\n', p, mat2str(prm, 10));
end
fclose(fid);
end

function a = read_bin(path, ny, nx, cls)
fid = fopen(path, 'r');
a = fread(fid, [ny nx], ['*' cls]);
fclose(fid);
a = double(a);
end
