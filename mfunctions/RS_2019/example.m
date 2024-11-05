% This script demonstrates how to use the tests proposed in Rossi and
% Sekhposyan "Alternative Tests for Correct Specification of Conditional 
% Predictive Densities," 2019, Journal of Econometrics 208(2), 638-657.
% Link here: https://www.sciencedirect.com/science/article/abs/pii/S0304407618302197
% Last tested on November 3, 2023, with Matlab R2022a

% How to use our codes, an empirical example. 
% This codes generates data from a standard normal distribution and then
% evaluates the standard normal predictive density for that realization.

% *************************************************************************
% Preliminaries
% *************************************************************************
clear all;
rng("default")

% Discretizes the grid for r, the support of the uniform distribution used
% to evaluate the PITs.
rvec = [0:0.001:1]; 

% *************************************************************************
% Generating the sample and the predictive density
% *************************************************************************
% Creates the sample. In this example the data is generated from a standard
% normal distribution and the predictive distribution is assumed to be a
% standard normal as well.
T = 200;
ytrue = randn(T,1);
forecast_mean = 0;
forecast_std = 1;

% Creates the probability integral transform (PITs), given the realization
% and the predictive density. The PITs are stored in a vector z.
z = [];
for i = 1:T
    z = [z; normcdf(ytrue(i,1),forecast_mean,forecast_std)];
end

% Implementing Rossi-Sekhposyan test
[crit_value boot_crit_value] = rs_test(z, rvec);
% pick the critical value, i.e. the last argument of plotsofuniformity by
% indicating either 1, 5, or 10.
plotsofuniformity(z,'Illustrative Example for One-Step-Ahead', crit_value,5);
plotsofuniformity(z,'Illustrative Example for Multi-Step-Ahead', boot_crit_value,5);
