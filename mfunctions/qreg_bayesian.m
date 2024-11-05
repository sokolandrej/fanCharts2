function [beta_draws, gamma_draws] = qreg_bayesian(y,x,quant,nsave,nburn,chainStep,varargin)

[T,p] = size(x);
ntot = nsave + nburn;   % Number of total draws

if ~isempty(varargin)
    
bQRoptions = varargin{1};

end

% prior for sigma2 ~ IG(a1,a2)
if exist('bQRoptions.a1')
    
    a1 = bQRoptions.a1;
    
else
    
    a1 = .01;
    
end

if exist('bQRoptions.a2')
    
    a2 = bQRoptions.a2;
    
else
    
    a2 = .01;
    
end

if exist('bQRoptions.d1')
    
    d1 = bQRoptions.d1;
    
else
    
    d1 = .01;
    
end

if exist('bQRoptions.d2')
    
    d2 = bQRoptions.d2;
    
else
    
    d2 = .01;
    
end


% option to shrink constant coefficients to unconditional quantiles rather
% than 0
if exist('bQRoptions.constPrior')
    
    constPrior = bQRoptions.constPrior;
    
else
    
    constPrior = 1;
    
end

a_1 = a1 + 3*T/2;
d_1 = d1+1/2;

% ==============| Initialize vectors
nq = length(quant);
beta = repmat(x\y,1,nq);
betaPrior = zeros(size(beta));

if all(x(:,1) == ones(size(x,1),1)) && constPrior
    
    betaPrior(1,:) = prctile(y,quant*100); %makes only sense with constant in x
    
end

z = ones(T,nq);
scaleInv = zeros(1,nq);
invd2 = .1*ones(p,nq);
theta = zeros(1,nq);
tau2 = zeros(1,nq);
k2 = tau2;
k1 = k2;
xTransp = x';

% ==============| Storage matrices
beta_draws = zeros(ntot,p,nq);

for q = 1:nq   % sample for each quantile
    
    tau2(q) = 2/(quant(q)*(1-quant(q)));
    theta(q) = (1-2*quant(q))/(quant(q)*(1-quant(q)));
    k1(q) = sqrt(theta(:,q).^2 + 2*tau2(:,q));
    k2(q) = (theta(:,q).^2 + 2*tau2(:,q))/tau2(:,q);
    
end

for irep = 1:ntot
    
    for q = 1:nq   % sample for each quantile
        
        % Sample regression variance sigma2
        sse = (y-x*beta(:,q) - theta(q)*z(:,q)).^2;
        a_2 = a2 + sum(sse./(2*z(:,q)*tau2(q))) + sum(z(:,q));
        scaleInv(q) = gamrnd(a_1,1/a_2);
              
        % Sample regression coefficients beta
        U =  scaleInv(q)./(tau2(q).*z(:,q)') ; % this is presumably what xtilde in the paper means? see 3.1 in Kozumi and Kobayashi
        y_tilde = y - theta(q)*z(:,q);
        xsq = U.*xTransp;
        V_beta = inv2(xsq*x + diag(invd2(:,q)));
        miu_beta = V_beta*(xsq*y_tilde + betaPrior(:,q).*invd2(:,q));
        draw = randn(p,1);
        beta(:,q) = miu_beta + chol(V_beta,'lower')*draw;
                               
        % Sample invd2
        d_2 = ((beta(:,q)-betaPrior(:,q)).^2)/2+d2;
        invd2(:,q) = gamrnd(d_1,1./d_2);
                
        % Sample latent variables z_{t}
        k1Temp = k1(q)./abs(y-x*beta(:,q));
        k2Temp = k2(q)*scaleInv(q);
        z(:,q) = 1./igrnd(k1Temp,k2Temp);
        
    end
    
    z = min(max(z,1e-4),1e10);
    
    beta_draws(irep,:,:) = beta;
    
end

beta_draws = permute(beta_draws(nburn+1:chainStep:end,:,:),[2 3 1]);

end
