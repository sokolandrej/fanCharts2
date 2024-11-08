function [X,f]=bfmdgp(L,Sigma,T,N,K)
%BFMDGP - Bayesian Factor Model DGP

if nargin==0
    T=100;
    N=9;
    K=3;

    L = [1     0     0;
         0     1     0;
         0     0     1;
         0.99  0     0;
         0.99  0     0;
         0     0.95  0;
         0     0.95  0;
         0     0     0.90;
         0     0     0.90;];

    Sigma = diag([0.02;0.19;0.36;0.02;0.02;0.19;0.19;0.36;0.36]);
end

f = randn(T,K);
X = f*L' + randn(T,N)*chol(Sigma);