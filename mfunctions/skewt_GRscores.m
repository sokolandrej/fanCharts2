function out = skewt_GRscores(y,params)
%skewt_GRscores Computes Gneiting and Ranjan (JBES, 2011) quantile-weighted
%scores.
%   based on eq. (8) and weighting functions in Table 2

% weighting functions
v0 = @(q) 1;
v1 = @(q) q*(1-q);
v2 = @(q) (2*q-1)^2;
v3 = @(q) q^2;
v4 = @(q) (1-q)^2;

%icdf
f = @(q) qskt(q,params(1),params(2),params(3),params(4));

% integrands
s0 = @(q) rho(y-f(q),q)*v0(q);
s1 = @(q) rho(y-f(q),q)*v1(q);
s2 = @(q) rho(y-f(q),q)*v2(q);
s3 = @(q) rho(y-f(q),q)*v3(q);
s4 = @(q) rho(y-f(q),q)*v4(q);

unweighted = integral(s0,1e-3,.999,'ArrayValued',1);
center = integral(s1,1e-3,.999,'ArrayValued',1);
tails = integral(s2,1e-3,.999,'ArrayValued',1);
rTail = integral(s3,1e-3,.999,'ArrayValued',1);
lTail = integral(s4,1e-3,.999,'ArrayValued',1);

out = [unweighted,center,tails,rTail,lTail];

end
