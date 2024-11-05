%% HLN test -  Harvey, Leybourne, and Newbold (1997) correction of Diebold-Mariano test

function [tHLN,pHLN,tDM,pDM]=DM_HLNtest(e1,e2,h,squares)
%HLNtest: Retrieves the Diebold-Mariano test statistic (1995) for the 
% equality of forecast accuracy of two forecasts under general assumptions 
% with Harvey, Leybourne, and Newbold (1997) correction.
%
%   [tstat,pvalue]=HLNtest(e1,e2,h) calculates the D-M test statistic 
%   on the base of the loss differential which is defined as the difference
%   of the squared forecast errors
%
%   In particular, with the DM statistic one can test the null hypothesis: 
%   H0: E(d) = 0. The Diebold-Mariano test assumes that the loss 
%   differential process 'd' is stationary and defines the statistic as:
%   DM = mean(d) / sqrt[ (1/T) * VAR(d) ]  ~ N(0,1),
%   where VAR(d) is a HAC estimate of the unconditional variance of 'd'.
%
%   This function also corrects for the autocorrelation that multi-period 
%   forecast errors usually exhibit. Note that an efficient h-period 
%   forecast will have forecast errors following MA(h-1) processes. 
%   Diebold-Mariano use a Newey-West type estimator for sample variance of
%   the loss differential to account for this concern.
%
%   'e1' is a 'T1-by-1' vector of the forecast errors from the first model
%   'e2' is a 'T2-by-1' vector of the forecast errors from the second model
%
%   It should hold that T1 = T2 = T.
%
%   [tstat,pvalue]=HLNtest(e1,e2,h) allows you to specify an additional 
%   parameter value 'h' to account for the autocorrelation in the loss 
%   differential for multi-period ahead forecasts.   
%       'h'        the forecast horizon, initially set equal to 1
%       'squares'  use squares of 'e1' and 'e2' when computing differential
%
%   [tstat,pvalue]=HLNtest(e1,e2,h) returns a value of the test statistic
%   and corresponding p-value:
%       'tstat'      value of the Diebold-Mariano (1995) test statistic
%                    with Harvey, Leybourne, and Newbold (1997) correction 
%       'pvalue'     p-value of the Diebold-Mariano (1995) test statistic
%                    with Harvey, Leybourne, and Newbold (1997) correction 

if nargin < 2
    error('HLNtest:TooFewInputs','At least two arguments are required');
end
if nargin < 3
    h = 1;
end
if size(e1,1) ~= size(e2,1) || size(e1,2) ~= size(e2,2)
    error('HLNtest:InvalidInput','Vectors should be of equal length');
end
if size(e1,2) > 1 || size(e2,2) > 1
    error('HLNtest:InvalidInput','Input should have T rows and 1 column');
end

% Define the loss differential
if squares
    differential = e1.^2 - e2.^2;
else
    differential = e1 - e2;
end
differential=differential(~isnan(differential));
t = size(differential,1);
  

% Calculate the variance of the loss differential, taking into account
% autocorrelation.
varD=NeweyWest(differential,h);

tDM=sqrt(t)*mean(differential)/sqrt(varD);
tHLN=sqrt((t+1-2*h+h*(h-1)/t)/t)*tDM;
x=normcdf(tDM,0,1);
pDM=1-2*abs(x-.5);
pHLN= 2*tcdf(-abs(tHLN),t-1);
