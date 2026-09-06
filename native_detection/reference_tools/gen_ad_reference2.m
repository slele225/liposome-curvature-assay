% Fine-grained MEX Anderson-Darling decisions near the threshold.
function gen_ad_reference2(outfile)
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'cmeAnalysis-master', 'software')));
fid = fopen(outfile, 'w');
w = 6; [X,Y] = meshgrid(-w:w);
n = 0;
for k = 1:6000
    sk = 0.02 + 1.5*mod(k*0.0137, 1);
    ph = k*0.7;
    fam = mod(k, 5);
    base = sin(1000*(X*13+Y)+ph) + 0.5*cos(37*X.*Y+ph);
    switch fam
        case 0, noise = base + sk*base.^3;
        case 1, noise = base + sk*abs(base);
        case 2, noise = sign(base).*abs(base).^(0.2+sk);
        case 3, noise = base + sk*(base > 0.5);
        case 4, noise = tanh(sk*3*base);
    end
    win = 100*exp(-((X-0.3).^2+(Y+0.2).^2)/(2*1.5^2)) + 50 + 5*noise;
    init = [0 0 max(win(:))-min(win(:)) 1.5 min(win(:))];
    [prm, prmStd, C, res] = fitGaussian2D(win, init, 'xyAc');
    r = res.data(:); r = r(~isnan(r));
    % quick pre-screen using the case-2 statistic so that only cases near the
    % threshold are stored
    m = mean(r); s = std(r); xs = sort(r); z = normcdf(xs, m, s); nn = numel(r); i = (1:nn)';
    A2 = -nn - sum((2*i-1).*(log(z) + log(1-z(nn+1-i))))/nn;
    if A2 > 2.0 && A2 < 2.7
        n = n + 1;
        fprintf(fid, '%d %d %.17g %.17g %.17g\n', k, res.hAD, res.mean, res.std, res.RSS);
        fprintf(fid, '%.17g ', res.data(:)); fprintf(fid, '\n');
    end
end
fclose(fid);
fprintf('stored %d cases\n', n);
end
