function [ LL ] = skewt_relative_entropy( condParams,uncondParams, cutoff, varargin)
% Relative entropy
%   computes relative entropy of unconditinoal vs conditinoal distribution
%   (something like the information loss associated with using the
%   unconditional rather than the conditional distribution to compute the
%   probability of tail events).

if isempty(varargin)
    
    upside = 0;
    
else
    
    upside = varargin{1};
    
end

f = @(x) dskt(x,condParams(1),condParams(2),condParams(3),condParams(4));
g = @(x) dskt(x,uncondParams(1),uncondParams(2),uncondParams(3),uncondParams(4));
fun = @(x) (log(g(x))-log(f(x))).*f(x);

intRange = -1e3:.1:1e3;
lb = find(~isnan(fun(intRange)),1,'first')+1;
ub = find(~isnan(fun(intRange)),1,'last');    


if upside
    
    %LL = -integral(fun,qskt(cutoff,condParams(1),condParams(2),condParams(3),condParams(4)),10e3);
    LL = -integral(fun,qskt(cutoff,uncondParams(1),uncondParams(2),uncondParams(3),uncondParams(4)),intRange(ub));
    
else
    
    %LL = -integral(fun,-10e3,qskt(cutoff,condParams(1),condParams(2),condParams(3),condParams(4)));
    LL = -integral(fun,intRange(lb),qskt(cutoff,uncondParams(1),uncondParams(2),uncondParams(3),uncondParams(4)));
    
end


end

