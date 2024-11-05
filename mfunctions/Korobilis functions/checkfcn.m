function [y] = checkfcn(x,p)

% Check function for Bayesian quantile regression. Inputs are the data x
% and p and the quantile level.
y = (x < 0)*-1*x*(1-p) + (x > 0)*p*x;