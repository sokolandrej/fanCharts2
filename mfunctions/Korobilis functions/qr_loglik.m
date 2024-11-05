function logLikelihood = qr_loglik(ytilde,x,bb,z,t)
%QR_LOGLIK Computes logLikelihood of quantile regression
%   Detailed explanation goes here

logLikelihood = -.5*sum(log(z))-.5*(sum(((ytilde-x*bb).^2)./(t*sqrt(z)).^2)); 

end

