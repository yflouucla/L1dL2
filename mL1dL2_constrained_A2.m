function [x,result] = mL1dL2_constrained_A2(A,b,pm)
%=============================================================
% the L1/L2 adaptive 2 (L1/L2-A2)
%
% Solves
%           min  norm(x,1)/norm(x,2)
%           s.t. Ax = b
%
% Reference: "Accelerated Schemes for the L1/L2 minimization"
%             Chao Wang, Ming Yan, Yifei Lou
% Available at:
%             xxxx
%
% Author: Chao Wang
% Date: May 19. 2019
%=============================================================

result = pm;
[~,N] = size(A);
AAt = A'*((A*A')\b);
B = eye(N) - A'*((A*A')\A);
xg = pm.xg;
N_inner = N;
reletive_error = 1e-8;
if isfield(pm,'reletive_error'); reletive_error = pm.reletive_error; end
if isfield(pm,'N_inner'); N_inner = pm.N_inner; end
iterDCA = 30;
if isfield(pm,'iterDCA'); iterDCA = pm.iterDCA; end

% L1 solution
if isfield(pm,'xr')
    xr = pm.xr;
else
    xr = mL1_constrained_LP_Gurobi(A,b);
end
alpha = norm(xr,1)/norm(xr);
alpha_list = alpha;
%  L1-L2 DCA
x = xr;
u = x;
y = x;
eps = 1e-14;
obj_list = [];
rho = 20;
beta = 1;
if isfield(pm,'beta'); beta = pm.beta; end
if isfield(pm,'rho'); rho = pm.rho; end
disy_list = [];

for i = 1:iterDCA
    xold = x;
    c = alpha*x/(norm(x) + eps);
    for j = 1:N_inner
        xold_inner = x;
        % x update
        x = B*((y*rho + xold*beta - u +c)/(rho+beta)) + AAt;
        % y update
        y = shrink(x+u/rho,1/rho);
        % u update
        u = u + rho*(x-y);
        if norm(x-xold_inner)/(norm(xold_inner)+eps) < 1e-6
            break;
        end
        disy_list = [disy_list, norm(x-y)];
    end
    obj_list = [obj_list, norm(x,1)-alpha*norm(x)];
    alpha = norm(x,1)/norm(x);
    alpha_list = [alpha_list,alpha];
    if norm(x-xold)/norm(x) < reletive_error
        break;
    end
end
result.i = i;
result.j = j;
result.alpha_list = alpha_list;
result.disy_list = disy_list;
result.xnorm2_r = (norm(x)-norm(xg))/norm(xg);
result.xnorm2 = norm(x)-norm(xg);
result.norm2 = norm(x);
result.obj_list = obj_list;
fx = norm(x,1)-norm(x);
result.fx = fx;
fxg = norm(xg,1)-norm(xg);
result.fxg = fxg;
result.alpha = alpha;

%% Evaluation
result.error = norm(x-xg)/norm(xg);
if result.error < pm.restol % success
    result.rate = 1;
elseif fx+ eps < fxg % model failure
    result.rate = -1;
else    % algorithm failure
    result.rate = -2;
end


end
