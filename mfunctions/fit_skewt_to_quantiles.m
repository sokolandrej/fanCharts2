function [fittedParams, fittedQuantiles] = fit_skewt_to_quantiles(lpQuantiles,quantiles,quantilesToFit)
% fits skew-t distribution to quantiles (see
% https://www.newyorkfed.org/medialibrary/media/research/staff_reports/sr794.pdf?la=en,
% p.9 and http://azzalini.stat.unipd.it/SN/ for the matlab code and other
% resources. This function allows both exactly-identified and
% over-identified fitting.

if length(quantilesToFit)<4
    
    error('Not enough quantiles to fit parameters')
    
end

options = optimoptions('fmincon','Display','notify','UseParallel',true);


horizon = size(lpQuantiles,1);
nQ = length(quantiles);
quantileIdxToFit = find(ismember(quantiles,quantilesToFit));
fittedParams = zeros(horizon,4);
fittedQuantiles = zeros(horizon, nQ);

lb = [-inf 1e-6 -1e6 1];
ub = [inf inf 1e6 inf];

for hh = 1:horizon
    
    minimand = @(x) sum((lpQuantiles(hh,quantileIdxToFit) - qskt(quantiles(quantileIdxToFit),x(1),x(2),x(3),round(x(4)))).^2);
    fittedParams(hh,:) = fmincon(minimand,[0 1 0 1],[],[],[],[],lb,ub,[],options);
    fittedParams(hh,4) = round(fittedParams(hh,4));
    fittedQuantiles(hh,:) = qskt(quantiles,fittedParams(hh,1),fittedParams(hh,2),fittedParams(hh,3),fittedParams(hh,4));
    
end


end