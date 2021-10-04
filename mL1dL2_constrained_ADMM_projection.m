function [x,result] = mL1dL2_constrained_ADMM_projection(A,b,pm)
%=============================================================
% ADMM for the L1/L2 constrained model 
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
% Date: Nov 11. 2018
%============================================================= 

[~,N] = size(A);
xg = pm.xg;
result = pm;   

rho = 100;
if isfield(pm,'rho'); rho = pm.rho; end

iter = 2*N;
if isfield(pm,'iter'); iter = pm.iter; end

reletive_error = 1e-8;
if isfield(pm,'reletive_error'); reletive_error = pm.reletive_error; end

% L1 solution
if isfield(pm,'xr')
    xr = pm.xr;
else
    xr = mL1_constrained_LP_Gurobi(A,b);
end
%% ---------------- L1/L2 ----------------
eps = 1e-9;
obj_list = [];
dis_list = [];disz_list =[]; disy_list =[];
normx_list = [];
x = xr;
y = xr;
z = xr;
v = zeros(N,1);
w = zeros(N,1);
AAt = A'*((A*A')\b);
B = eye(N) - A'*((A*A')\A);
for i = 1:iter
    x_old = x;
    % x update
    xx = 0.5*(y+z) - (1/(2*rho))*(v+w);
    x = B*xx + AAt;
    % y update
    d = x + v/rho;
    t = mfindroot(norm(z,1),norm(d),rho);
    y = t*d;
    % z update
    z = mShrink(x + w/rho, 1/(rho*norm(y)));
    % u, v, w updates
    v = v + rho*(x - y);
    w = w + rho*(x - z);
    if norm(x-x_old)/norm(x_old) < reletive_error && i>1500
        break;
    end

    obj_list = [obj_list, norm(x,1)/norm(x)];
    dis_list = [dis_list, max(norm(x-y),norm(x-z))];
    disz_list = [disz_list, norm(x-z)];
    disy_list = [disy_list, norm(x-y)];
    normx_list = [normx_list, norm(x)];
end
result.i = i;
fx = norm(x,1)/norm(x);
result.fx = fx;
fxg = norm(xg,1)/norm(xg);
result.fxg = fxg;
result.obj_list = obj_list;
result.dis_list = dis_list;
result.disz_list = disz_list;
result.disy_list = disy_list;
result.normx_list = normx_list;

end

function v = mShrink(s,lambda)

v = sign(s).*(max(abs(s) - lambda , 0));

end


function [tau] = mfindroot(c,eta,rho)

if c == 0 || eta == 0 
    tau = 1;
else
    a = 27*c/(rho*(eta^3)) + 2;
    C = ((a + (a^2 - 4)^0.5)/2)^(1/3);
    tau = (1 + C + 1/C)/3;
end

end