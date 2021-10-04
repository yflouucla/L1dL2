%=============================================================
% demo codes for the L1/L2 constrained model 
%
% Solves
%           min  norm(x,1)/norm(x,2)
%           s.t. Ax = b
%
% Reference: 
% 1. Yaghoub Rahimi, Chao Wang, Hongbo Dong, and Yifei Lou, 
%    "A Scale Invariant Approach for Sparse Signal Recovery," SISC (2019) 
%              
% 2. Chao Wang, Ming Yan, Yaghoub Rahimi, and Yifei Lou, 
%    "Accelerated Schemes for the L1/L2 Minimization," TSP (2020) 
%
% 
% Author: Chao Wang  
% Date: May 5, 2020
%============================================================= 

Please install Gurobi (https://www.gurobi.com/) for L1 minimization or replace mL1_constrained_LP_Gurobi.m with any L1 solver as an initial condition for the L1/L2 model. 


Simply run

demo.m 		% for sparse recovery
demoMRI.m 	% for MRI reconstruction
