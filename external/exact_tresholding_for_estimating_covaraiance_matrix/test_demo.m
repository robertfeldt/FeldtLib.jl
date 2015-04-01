%create data set for testing code compute_GL_exact.m

% Creates a sample covariance matrix as described in our paper experiment:
% Exact Covariance Thresholding into Connected Components for large-scale Graphical Lasso
% by Mazumder and Hastie, 2011


addpath('SMOOTH_COVSEL_Matlab/')
clear;
clc;

% p1: size of each block 
p1=20; 
no_blocks=2; % number of blocks

sigma_cov=ones(p1,p1);
Sigma_cov=[];

min_value=min(sigma_cov(:));

for i=1:no_blocks;
Sigma_cov=blkdiag(Sigma_cov,sigma_cov);
end

p=size(Sigma_cov,1);
zero_ids=(Sigma_cov==0); 

% add noise to suppress the block-diagonal structure

NSE=randn(size(Sigma_cov)); NSE=NSE*NSE'/(10*p);
max_value=max(NSE(:)); NSE=NSE*min_value*0.8/max_value;
Sigma_cov = Sigma_cov + NSE;

tmp=Sigma_cov.*(1-zero_ids); tmp=tmp(tmp~=0);
rho_val_max=min(tmp(:))
clear tmp;
tmp=Sigma_cov.*zero_ids; tmp=tmp(tmp~=0);
rho_val_min=max(tmp(:));

% any rho in [rho_val_min, rho_val_max] will have atleast two connected
% components


% we take the following convex combination.

rho=rho_val_min*.01 + rho_val_max*.99;

%% calling sequence

% no specification for convergence of smoothcovsel_main specified
% they will be at the default values, see 
%>> help compute_GL_exact
% for details on specifying the convergence criterion 

[Time, X, S,C] = compute_GL_exact(Sigma_cov,rho); 

 




