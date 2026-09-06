% Generates reference values from MATLAB (incl. the fitGaussian2D MEX) for cpp_port/tests/unit_tests.cpp
function gen_unit_reference(outfile)
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'cmeAnalysis-master', 'software')));
fid = fopen(outfile, 'w');
p = @(varargin) fprintf(fid, varargin{:});
p('norminv_0975 %.17g\n', norminv(1-0.05/2));
xs = [-3.2 -1.5 -0.7 0 0.4 2.5 8];
nus = [0.5 1 2.7 10.5 168 1e8];
for a = 1:numel(xs), for b = 1:numel(nus)
    p('tcdf %.17g %.17g %.17g\n', xs(a), nus(b), tcdf(xs(a), nus(b)));
end, end
p('tcdf_nan0 %d\n', isnan(tcdf(0,0)));
p('normcdf %.17g\n', normcdf(0.3, 0.1, 1.7));
x = sin((1:60)*1.7) + 0.3*cos((1:60)*0.4);
[H,~,A2] = adtest1(x, 0.05);
p('adtest_H %d\nadtest_A2 %.17g\n', H, A2);
x2 = x.^3;
[H,~,A2] = adtest1(x2, 0.05);
p('adtest2_H %d\nadtest2_A2 %.17g\n', H, A2);
% synthetic window for fitGaussian2D
w = 6; [X,Y] = meshgrid(-w:w);
true = [0.3 -0.2 100 1.5 50];
noise = 4*sin(1000*(X*13+Y)) + 2*cos(37*X.*Y);
win = true(3)*exp(-((X-true(1)).^2+(Y-true(2)).^2)/(2*true(4)^2)) + true(5) + noise;
p('win %d\n', numel(win));
p('%.17g ', win(:)); p('\n');
init = [0 0 max(win(:))-min(win(:)) 1.5 min(win(:))];
modes = {'xyAc', 'xyasc', 'Ac', 'xyac'};
for m = 1:numel(modes)
    [prm, prmStd, C, res] = fitGaussian2D(win, init, modes{m});
    p('fit %s prm %.17g %.17g %.17g %.17g %.17g\n', modes{m}, prm);
    p('fit %s prmStd', modes{m}); p(' %.17g', prmStd); p('\n');
    p('fit %s res %.17g %.17g %.17g %d\n', modes{m}, res.mean, res.std, res.RSS, res.hAD);
    p('fit %s C', modes{m}); p(' %.17g', C(:)); p('\n');
end
% masked window (NaNs)
win2 = win; win2(1:3, 1:4) = NaN; win2(10, 12) = NaN;
[prm, prmStd, C, res] = fitGaussian2D(win2, init, 'xyAc');
p('fitnan xyAc prm %.17g %.17g %.17g %.17g %.17g\n', prm);
p('fitnan xyAc prmStd'); p(' %.17g', prmStd); p('\n');
p('fitnan xyAc res %.17g %.17g %.17g %d\n', res.mean, res.std, res.RSS, res.hAD);
% asymmetric window to pin down x/y convention: shift only in x
win3 = true(3)*exp(-((X-1.2).^2+(Y+0.0).^2)/(2*1.5^2)) + 50 + noise;
[prm] = fitGaussian2D(win3, init, 'xyAc');
p('fitx xyAc prm %.17g %.17g %.17g %.17g %.17g\n', prm);
% fitGaussians2D on a synthetic image with two spots
img = 100*ones(40, 50);
[XI, YI] = meshgrid(1:50, 1:40);
img = img + 300*exp(-((XI-20.3).^2+(YI-15.6).^2)/(2*1.4^2)) + 200*exp(-((XI-35.7).^2+(YI-28.2).^2)/(2*1.6^2)) + 5*sin(1000*(XI*40+YI));
ps = fitGaussians2D(img, [20 36], [16 28], [], [1.5 1.5], [], 'xyAc');
p('fg2d x %.17g %.17g\nfg2d y %.17g %.17g\nfg2d A %.17g %.17g\nfg2d c %.17g %.17g\n', ps.x, ps.y, ps.A, ps.c);
p('fg2d A_pstd %.17g %.17g\nfg2d sigma_r %.17g %.17g\nfg2d pval_Ar %.17g %.17g\nfg2d hval_AD %d %d\nfg2d mask_Ar %.17g %.17g\n', ps.A_pstd, ps.sigma_r, ps.pval_Ar, ps.hval_AD, ps.mask_Ar);
ps = fitGaussians2D(img, [20 36], [16 28], [], [1.5 1.5], [], 'Ac');
p('fg2dAc A %.17g %.17g\nfg2dAc c %.17g %.17g\nfg2dAc pval_Ar %.17g %.17g\n', ps.A, ps.c, ps.pval_Ar);
% pointSourceDetection on the synthetic image
[pstruct, mask, imgLM, imgLoG] = pointSourceDetection(img, 1.5);
p('psd n %d\n', numel(pstruct.x));
p('psd x'); p(' %.17g', pstruct.x); p('\n');
p('psd y'); p(' %.17g', pstruct.y); p('\n');
p('psd A'); p(' %.17g', pstruct.A); p('\n');
p('psd c'); p(' %.17g', pstruct.c); p('\n');
p('psd RSS'); p(' %.17g', pstruct.RSS); p('\n');
p('psd masksum %d\n', sum(mask(:)));
p('psd loG %.17g %.17g %.17g\n', imgLoG(1,1), imgLoG(16,20), imgLoG(40,50));
% img values for C++ side
p('img %d %d\n', size(img,1), size(img,2));
p('%.17g ', img(:)); p('\n');
% padarrayXT
a = reshape(1:12, 3, 4);
b = padarrayXT(a, [2 2], 'symmetric');
p('padxt %d %d\n', size(b,1), size(b,2)); p('%g ', b(:)); p('\n');
% gmdistribution.fit reference on deterministic data
rng(1);
xg = [1.4 + 0.2*sin((1:150)*2.1), 2.6 + 0.3*cos((1:80)*1.3), 1.9+0.05*sin((1:30)*5)];
for k = 1:3
    obj = gmdistribution.fit(xg', k, 'Options', statset('maxIter', 200));
    p('gmm %d mu', k); p(' %.17g', obj.mu); p('\n');
    p('gmm %d Sigma', k); p(' %.17g', squeeze(obj.Sigma)); p('\n');
    p('gmm %d P', k); p(' %.17g', obj.PComponents); p('\n');
    p('gmm %d NlogL %.17g BIC %.17g Iters %d Conv %d\n', k, obj.NlogL, obj.BIC, obj.Iters, obj.Converged);
end
p('gmm_data %d\n', numel(xg)); p('%.17g ', xg); p('\n');
p('rand_after %.17g\n', rand);
fclose(fid);
end
