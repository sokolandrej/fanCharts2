function [X]=bfmdgp2(L,Sigma,T,N,K)
%BFMDGP - Bayesian Factor Model DGP

if nargin==0
    T=300;
    N=7;
    K=1;

    L = [0.995; 0.975; 0.949; 0.922; 0.894; 0.866; 0.837];

    Sigma = diag([0.01; 0.05; 0.10; 0.15; 0.20; 0.25; 0.30]);
end

f=randn(T,K);
X = f*L' + randn(T,N)*chol(Sigma);