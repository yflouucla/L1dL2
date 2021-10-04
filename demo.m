%=============================================================
% demo ---- ADMM for the L1/L2 constrained model 
%
% Solves
%           min  norm(x,1)/norm(x,2)
%           s.t. Ax = b
%
% Reference: "A Scale Invariant Approach for Sparse Signal Recovery" 
%             Yaghoub Rahimi, Chao Wang, Hongbo Dong, Yifei Lou 
% Available at: 
%             https://arxiv.org/abs/1812.08852
% 
% Author: Chao Wang  
% Date: Feb. 21 2018
%============================================================= 

clear; close all

%% parameter settings
pm.M = 64;  pm.N = 1024; % matrix dimension m-by-n
pm.sen_mat = 'Oversampled_DCT';
pm.restol = 1e-3;
pm.cond1 = 0;
pm.dynamic = 0;
pm.normalized = 0;  
pm.F = 15;
pm.sparsity = 15;
pm.range_factor = 3;
pm.min_seperation = 2*pm.F; 
nt = 12; % 12 for dynmaic 0
%% -------- Simulation -----------------------
rng(nt*10)
A = data_generator_A(pm);
xg = data_generator_xg(pm); pm.xg = xg;
b = A*xg;
xL1 = mL1_constrained_LP_Gurobi(A,b);
pm.xr = xL1; % initial guess
%% -------- L1/L2 via ADMM (SISC paper)----------
[~,Result_p] = mL1dL2_constrained_ADMM_projection(A,b,pm);

% plot residual errors (SISC paper Fig. 2)
figure
semilogy(1:length(Result_p.disy_list),Result_p.disy_list,1:length(Result_p.disy_list),Result_p.disz_list,'--','LineWidth',2);
axis([1 length(Result_p.obj_list) 1e-10 1e-1])
xlabel('iteration')
legend('distance y','distance z')

%% -------- L1/L2 via DCA (TSP paper) ----------
tStart = tic;
[~,Result_a1] = mL1dL2_constrained_A1(A,b,pm);
[~,Result_b] = mL1dL2_constrained_BS(A,b,pm);
pm.beta = 1e-5; pm.rho =.3;
[~,Result_a2] = mL1dL2_constrained_A2(A,b,pm);

%% plot alpha in each iteration (TSP paper Fig. 1)
figure;
semilogy(Result_b.alpha_list,'-bo','LineWidth',2);
hold on
semilogy(Result_a1.alpha_list,'dr-','LineWidth',2);
hold on
semilogy(Result_a2.alpha_list,'-k>','LineWidth',2);
xlabel('\it k','fontsize',12)
LEG = legend('BS','A1','A2', 'Location','southeast');set(LEG,'FontSize',12);
axis([1 10 0 3.6]) 


