function [x,result] = mL1dL2_constrained_A1(A,b,pm)
%=============================================================
% the L1/L2 adaptive 1 (L1/L2-A1)
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
xg = pm.xg;

reletive_error = 1e-8;
if isfield(pm,'reletive_error'); reletive_error = pm.reletive_error; end

iterDCA = 40;
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
eps = 1e-9;
obj_list = [];
unbound = zeros(1,iterDCA);
for i = 1:iterDCA
    xold = x;
    c = alpha*x/(norm(x) + eps);
    f = ones(2*N,1) - [c; -c];
    B = [A -A];
    % Copyright 2017, Gurobi Optimization, Inc.
    clear model;
    model.A = sparse(B);
    model.obj = f;
    model.rhs = b;
    model.lb = zeros(2*N,1);
    model.sense = '=';
    clear params;
    params.outputflag = 0;
    results = gurobi(model,params);
    if strcmp(results.status,'UNBOUNDED') || strcmp(results.status,'INF_OR_UNBD')
        break;
        unbound(i) = 1;
    end
    if isfield(results,'x')
        x2 = results.x();
    else
        fprintf(results.status)
        fprintf(',No results... \n')
        break;
    end
    x = x2(1:N) - x2(N+1:2*N);
    obj_list = [obj_list, norm(x,1)-alpha*norm(x)];
    alpha = norm(x,1)/norm(x);
    alpha_list = [alpha_list,alpha];
    if norm(x-xold)/norm(x) < reletive_error
        break;
    end
end
result.i = i;
result.xnorm2_r = (norm(x)-norm(xg))/norm(xg);
result.xnorm2 = norm(x)-norm(xg);
result.norm2 = norm(x);
result.unbound = sum(unbound);
result.alpha_list = alpha_list;
result.obj_list = obj_list;
fx = norm(x,1)/norm(x);
result.fx = fx;
fxg = norm(xg,1)/norm(xg);
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
