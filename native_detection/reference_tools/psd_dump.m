function [pstruct, mask, dbg] = psd_dump(img, sigma, varargin)
% psd_dump: verbatim copy of cmeAnalysis pointSourceDetection.m (Francois Aguet,
% GPL-3.0, Danuser Lab) with the mixture branch removed and debug hooks added
% (struct 'dbg').  Used only to produce golden intermediates for cpp_port.

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('img', @isnumeric);
ip.addRequired('sigma', @isscalar);
ip.addParamValue('Mode', 'xyAc', @ischar);
ip.addParamValue('Alpha', 0.05, @isscalar);
ip.addParamValue('Mask', [], @(x) isnumeric(x) || islogical(x));
ip.addParamValue('RemoveRedundant', true, @islogical);
ip.addParamValue('RedundancyRadius', 0.25, @isscalar);
ip.addParamValue('Prefilter', true, @islogical);
ip.addParamValue('RefineMaskLoG', true, @islogical);
ip.addParamValue('RefineMaskValid', true, @islogical);
ip.addParamValue('ConfRadius', []);
ip.addParamValue('WindowSize', []);
ip.KeepUnmatched = true;
ip.parse(img, sigma, varargin{:});
mode = ip.Results.Mode;
alpha = ip.Results.Alpha;

dbg = struct('imgLoG', [], 'A_est', [], 'c_est', [], 'pval', [], 'maskPrefilter', [], 'maskCombined', [], ...
    'logThreshold', NaN, 'lmx', [], 'lmy', [], 'fitAll', [], 'keepNonNaN', [], 'keepFinal', []);

if ~isa(img, 'double')
    img = double(img);
end

w = ceil(4*sigma);
x = -w:w;
g = exp(-x.^2/(2*sigma^2));
u = ones(1,length(x));

imgXT = padarrayXT(img, [w w], 'symmetric');
fg = conv2(g', g, imgXT, 'valid');
fu = conv2(u', u, imgXT, 'valid');
fu2 = conv2(u', u, imgXT.^2, 'valid');

gx2 = g.*x.^2;
imgLoG = 2*fg/sigma^2 - (conv2(g, gx2, imgXT, 'valid')+conv2(gx2, g, imgXT, 'valid'))/sigma^4;
imgLoG = imgLoG / (2*pi*sigma^2);

g = g'*g;
n = numel(g);
gsum = sum(g(:));
g2sum = sum(g(:).^2);

A_est = (fg - gsum*fu/n) / (g2sum - gsum^2/n);
c_est = (fu - A_est*gsum)/n;
dbg.imgLoG = imgLoG; dbg.A_est = A_est; dbg.c_est = c_est;

if ip.Results.Prefilter
    J = [g(:) ones(n,1)];
    C = inv(J'*J);
    f_c = fu2 - 2*c_est.*fu + n*c_est.^2;
    RSS = A_est.^2*g2sum - 2*A_est.*(fg - c_est*gsum) + f_c;
    RSS(RSS<0) = 0;
    sigma_e2 = RSS/(n-3);
    sigma_A = sqrt(sigma_e2*C(1,1));
    sigma_res = sqrt(RSS/(n-1));
    kLevel = norminv(1-alpha/2.0, 0, 1);
    SE_sigma_c = sigma_res/sqrt(2*(n-1)) * kLevel;
    df2 = (n-1) * (sigma_A.^2 + SE_sigma_c.^2).^2 ./ (sigma_A.^4 + SE_sigma_c.^4);
    scomb = sqrt((sigma_A.^2 + SE_sigma_c.^2)/n);
    T = (A_est - sigma_res*kLevel) ./ scomb;
    pval = tcdf(-T, df2);
    mask = pval < 0.05;
    dbg.pval = pval;
else
    mask = true(size(img));
    dbg.pval = zeros(size(img));
end
dbg.maskPrefilter = mask;

allMax = locmax2d(imgLoG, 2*ceil(sigma)+1);
imgLM = allMax .* mask;

pstruct = [];
dbg.maskCombined = mask;
if sum(imgLM(:))~=0
    if ip.Results.RefineMaskLoG
        logThreshold = min(imgLoG(imgLM~=0));
        logMask = imgLoG >= logThreshold;
        mask = mask | logMask;
        dbg.logThreshold = logThreshold;
    end
    imgLM = allMax .* mask;
    if ~isempty(ip.Results.Mask)
        imgLM(ip.Results.Mask==0) = 0;
    end
    dbg.maskCombined = mask;
    [lmy, lmx] = find(imgLM~=0);
    lmIdx = sub2ind(size(img), lmy, lmx);
    dbg.lmx = lmx; dbg.lmy = lmy;
    if ~isempty(lmIdx)
        pstruct = fitGaussians2D(img, lmx, lmy, A_est(lmIdx), sigma*ones(1,length(lmIdx)),...
            c_est(lmIdx), mode, 'mask', mask, 'alpha', alpha,...
            'ConfRadius', ip.Results.ConfRadius, 'WindowSize', ip.Results.WindowSize);
        dbg.fitAll = pstruct;
        idx = ~isnan([pstruct.x]);
        dbg.keepNonNaN = idx;
        if sum(idx)~=0
            fnames = fieldnames(pstruct);
            for k = 1:length(fnames)
                pstruct.(fnames{k}) = pstruct.(fnames{k})(idx);
            end
            idx = [pstruct.hval_Ar] == 1;
            if ip.Results.RemoveRedundant
                pM = [pstruct.x' pstruct.y'];
                idxKD = KDTreeBallQuery(pM, pM, ip.Results.RedundancyRadius*ones(numel(pstruct.x),1));
                idxKD = idxKD(cellfun(@numel, idxKD)>1);
                for k = 1:length(idxKD)
                    RSS = pstruct.RSS(idxKD{k});
                    idx(idxKD{k}(RSS ~= min(RSS))) = 0;
                end
            end
            dbg.keepFinal = idx;
            if sum(idx)>0
                fnames = fieldnames(pstruct);
                for k = 1:length(fnames)
                    pstruct.(fnames{k}) = pstruct.(fnames{k})(idx);
                end
                pstruct.hval_Ar = logical(pstruct.hval_Ar);
                pstruct.hval_AD = logical(pstruct.hval_AD);
                pstruct.isPSF = ~pstruct.hval_AD;
            else
                pstruct = [];
            end
        else
            pstruct = [];
        end
    end
end

if ~isempty(pstruct) && ip.Results.RefineMaskValid
    CC = bwconncomp(mask);
    labels = labelmatrix(CC);
    loclabels = labels(sub2ind(size(img), pstruct.y_init, pstruct.x_init));
    idx = setdiff(1:CC.NumObjects, loclabels);
    CC.PixelIdxList(idx) = [];
    CC.NumObjects = length(CC.PixelIdxList);
    mask = labelmatrix(CC)~=0;
end
