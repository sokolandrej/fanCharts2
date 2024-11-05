% This script demonstrates how to use the tests proposed in Rossi and
% Sekhposyan "Alternative Tests for Correct Specification of Conditional 
% Predictive Densities," 2019, Journal of Econometrics 208(2), 638-657.
% Link here: https://www.sciencedirect.com/science/article/abs/pii/S0304407618302197


% example.m generates the data forecast densities and realizations based on a
% standard normal distrbution

% rs_test.m implements the two tests proposed in the paper, the one proposed 
% in Theorem 2 and Theorem 3, implemented based on the bootstrap proposed in
% Theorem 4.

% CVfinalbootstrapInoueExample.m is implementing the bootstrap proposed in
% Theorem 4, only needed if evaluating multiple-step-ahead forecasts