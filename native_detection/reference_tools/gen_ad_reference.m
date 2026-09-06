% Residuals + hAD decisions from the original fitGaussian2D MEX for many
% synthetic windows, to identify the exact Anderson-Darling variant used.
function gen_ad_reference(outfile)
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'cmeAnalysis-master', 'software')));
fid = fopen(outfile, 'w');
w = 6; [X,Y] = meshgrid(-w:w);
for k = 1:400
    amp = 1 + 15*mod(k*0.37, 1);
    sk = mod(k*0.61, 1);
    ph = k*0.7;
    base = sin(1000*(X*13+Y)+ph) + 0.5*cos(37*X.*Y+ph);
    noise = amp*(base + sk*base.^3 + 0.3*sk*abs(base));
    if mod(k, 7) == 0, noise = amp*sign(base).*abs(base).^0.3; end
    if mod(k, 11) == 0, noise = amp*(base > 0.6); end
    win = 100*exp(-((X-0.3).^2+(Y+0.2).^2)/(2*1.5^2)) + 50 + noise;
    init = [0 0 max(win(:))-min(win(:)) 1.5 min(win(:))];
    [prm, prmStd, C, res] = fitGaussian2D(win, init, 'xyAc');
    fprintf(fid, '%d %d %.17g %.17g %.17g\n', k, res.hAD, res.mean, res.std, res.RSS);
    fprintf(fid, '%.17g ', res.data(:)); fprintf(fid, '\n');
end
fclose(fid);
end
