function res = qbvar_alt(X,lags,stationary,lambda0,theta0,miu0,alpha0,ar,varargin)
% Estimates a blocked BVAR (BBVAR)

% Preliminaries
if ~isempty(varargin)
    
    doDensity = 1;
    nDraws = varargin{1};
    nSave = varargin{2};
    
else
    
    doDensity = 0;
    nDraws = 10000;
    
end

% Prepare data
endEstimT = find(~isnan(sum(X,2)),1,'last');

if ar
    
    yQ = X(1:endEstimT-1,:);
    
else
    
    yQ = X(1:endEstimT,:);

end

optNaN.method = 1;
optNaN.k = 3;
yQ = remNaNs_spline(yQ,optNaN);
    
estimRes = bvarGLP_alt(yQ,lags,'mcmc',doDensity,'pos',stationary,'MCMCfcast',0,'Fcast',0,'MNalpha',0,'Ndraws',2*nDraws,'MCMCconst',1,'lambda0',lambda0,'theta0',theta0,'miu0',miu0,'alpha0',alpha0,'sur',0,'MNPsi',1);
% estimRes = bvarGLP_fixedhyp(yQ,lags,'mcmc',doDensity,'pos',stationary,'MCMCfcast',0,'Fcast',0,'MNalpha',1,'Ndraws',2*nDraws,'MCMCconst',1,'lambda0',lambda0,'theta0',theta0,'miu0',miu0,'alpha0',alpha0,'sur',1,'MNpsi',0,'hyperpriors',0);
betaHat = estimRes.postmax.betahat;
sigmaHat = estimRes.postmax.sigmahat;

% disp('Running Kalman smoother')
% [X_sm,logLik] = run_Ksmoother(betaHat, sigmaHat, lags,X,yQ);
[X_sm,logLik] = run_Ksmoother_FM(betaHat, sigmaHat, lags,X(endEstimT-lags:end,:));
X_sm = [X(lags+1:endEstimT-1,:);X_sm];

if doDensity
    
    yMCMC = NaN([size(X(lags+1:end,:)) nSave]);
    auxT = size(X(endEstimT-lags:end,:),1);
    yMCMC(end-auxT+lags+1:end,:,:) = run_Ssmoother_block(estimRes.mcmc.beta,estimRes.mcmc.sigma,nDraws/nSave,lags,X(endEstimT-lags:end,:));
    
    for ss = 1:size(yMCMC,3)
        
        yMCMC(1:end-auxT+lags,:,ss) = X(lags+1:endEstimT-1,:);
        
    end
    
%     disp('Running simulation smoother')
%     yMCMC = run_Ssmoother(estimRes.mcmc.beta,estimRes.mcmc.sigma,nDraws/nSave,lags,X,yQ);
    res.yMCMC = yMCMC;
    res.beta = estimRes.mcmc.beta;
    res.sigma = estimRes.mcmc.sigma;
    
end

res.X_sm = X_sm;
res.yQ = yQ;
res.betaHat = betaHat;
res.sigmaHat = sigmaHat;
res.logLik = logLik;
res.lambda = estimRes.postmax.lambda;
res.theta = estimRes.postmax.theta;
res.miu = estimRes.postmax.miu;
res.alpha = estimRes.postmax.alpha;

end