function probe_mex_options(imgPath, outfile)
% Identify the default [maxIter eAbs eRel] of the fitGaussian2D MEX by
% re-running the free-sigma refit of getGaussianPSFsigmaFromData with
% explicit option vectors and counting exact matches with the default call.
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'cmeAnalysis-master', 'software')));
img = double(imread(imgPath));
pstruct = pointSourceDetection(img, 1.5, 'Mode', 'xyac');
np = numel(pstruct.x);
xi = round(pstruct.x); yi = round(pstruct.y);
w4 = 6; [ny, nx] = size(img);
opts = {[], [500 1e-8 1e-8], [1000 1e-8 1e-8], [200 1e-8 1e-8], [100 1e-8 1e-8], [500 1e-6 1e-6], [500 1e-10 1e-10], ...
        [500 1e-8 1e-6], [500 1e-6 1e-8], [500 0 1e-8], [500 1e-8 0], [500 1e-9 1e-9], [500 1e-7 1e-7], [50 1e-8 1e-8], [10000 1e-8 1e-8]};
names = {'default', '500,1e-8,1e-8', '1000,1e-8,1e-8', '200,1e-8,1e-8', '100,1e-8,1e-8', '500,1e-6,1e-6', '500,1e-10,1e-10', ...
        '500,1e-8,1e-6', '500,1e-6,1e-8', '500,0,1e-8', '500,1e-8,0', '500,1e-9,1e-9', '500,1e-7,1e-7', '50,1e-8,1e-8', '10000,1e-8,1e-8'};
res = cell(1, numel(opts));
for o = 1:numel(opts)
    R = NaN(np, 5);
    for p = 1:np
        if (xi(p)>w4 && xi(p)<=nx-w4 && yi(p)>w4 && yi(p)<=ny-w4)
            window = img(yi(p)-w4:yi(p)+w4, xi(p)-w4:xi(p)+w4);
            init = [pstruct.x(p)-xi(p) pstruct.y(p)-yi(p) pstruct.A(p) 1.5 pstruct.c(p)];
            if isempty(opts{o})
                prm = fitGaussian2D(window, init, 'xyasc');
            else
                prm = fitGaussian2D(window, init, 'xyasc', opts{o});
            end
            R(p,:) = prm;
        end
    end
    res{o} = R;
end
fid = fopen(outfile, 'w');
for o = 2:numel(opts)
    same = all(res{o} == res{1} | (isnan(res{o}) & isnan(res{1})), 2);
    fprintf(fid, '%-18s exact matches with default: %d / %d\n', names{o}, sum(same), np);
    fprintf('%-18s exact matches with default: %d / %d\n', names{o}, sum(same), np);
end
% dump windows + init + default result of every candidate for the C++ side
fprintf(fid, 'CANDIDATES %d\n', np);
for p = 1:np
    if (xi(p)>w4 && xi(p)<=nx-w4 && yi(p)>w4 && yi(p)<=ny-w4)
        window = img(yi(p)-w4:yi(p)+w4, xi(p)-w4:xi(p)+w4);
        init = [pstruct.x(p)-xi(p) pstruct.y(p)-yi(p) pstruct.A(p) 1.5 pstruct.c(p)];
        fprintf(fid, 'C %d init %.17g %.17g %.17g %.17g %.17g prm %.17g %.17g %.17g %.17g %.17g\n', p, init, res{1}(p,:));
        fprintf(fid, '%.17g ', window(:)); fprintf(fid, '\n');
    end
end
fclose(fid);
end
