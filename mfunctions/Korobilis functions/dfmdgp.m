function [X]=dfmdgp(L,Sigma,T,N,K)
%BFMDGP - Bayesian Factor Model DGP

if nargin==0
    T=200;
    N=9;
    K=3;
    lag_f=1;

    L = [1  0  0;
         0  1  0;
         0  0  1;
         0.99  0.65  0.34;
         0.99  0  0;
         0.78  0  0.54;
         0.35  0.85  0.78;
         0  0.33  0.90;
         0.78  0  0.90];

    Sigma = diag([0.2;0.19;0.36;0.2;0.2;0.19;0.19;0.36;0.36]);
    
    PHI = [0.5    0    0 ;
           0    0.5    0;
           0    0    0.5];
    
    PSI = [1   0.5   0.5;
           0   1     0.5;
           0   0     1];
    
     Q = inv(PSI*PSI');
end

f =[rand(lag_f,K) ; zeros(T,K)];
% Now generate f from VAR (L,PHI,PSI)
for nn = lag_f:T+lag_f
    u = chol(Q)'*randn(K,1);
    flag = mlag(f,lag_f);
    f(nn,:) = flag(nn,:)*PHI + u';
end
f=f(lag_f+1:T+lag_f,:);
X = f*L' + randn(T,N)*chol(Sigma)';
