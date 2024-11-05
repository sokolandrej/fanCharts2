function [ fittedDist ] = fit_kernel_to_quantiles( targetQuantiles, tau,st,varargin)
% fits kernel density to quantiles

if ~isempty(varargin)
    
    krnl = varargin{1};
    
else
    
    krnl = 'normal';
    
end

objective = @(q)  sum((targetQuantiles-icdf(fitdist(q,'Kernel','Kernel',krnl),tau)).^2);
solution = fmincon(objective,targetQuantiles(1:st:end)',[],[],[],[],[],[],[],optimoptions('fmincon','Display','off'));
fittedDist = fitdist(solution,'Kernel');

end

