function [x,result] = mL1dL2_constrained_BS(A,b,pm)
%=============================================================
% L1/L2-BS
%

result = pm;


[~,N] = size(A);
xg = pm.xg; % ground truth
if isfield(pm,'reletive_error'); reletive_error = pm.reletive_error; else
    reletive_error = 1e-8; end
if isfield(pm,'xr'); xr = pm.xr; else
    xr = mL1_constrained_LP_Gurobi(A,b); end
x = xr;
al = 1; ar = norm(x,1)/norm(x);
t = 0;
alpha = 1;
obj_list = [];
alpha_list = [];
xmin = x;
%% outer iteration
for i = 1:7
    x = xr;
    %% inner iteration
    for j=1:8
        xold = x;
        c = alpha*x/(norm(x));
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
            t = 1;
            break;
        end
        if isfield(results,'x')
            x2 = results.x();
        else
            fprintf(results.status)
            fprintf('  ,No results... \n')
            break;
        end
        x = x2(1:N) - x2(N+1:2*N);
        
        if norm(x,1)/norm(x) < norm(xmin,1)/norm(xmin)
            xmin = x;
        end
        if norm(x-xold)/norm(x) < reletive_error
            break;
        end
    end
    %% update alpha
    if t==1
        t=0;
        ar = alpha;
    elseif norm(x,1)-alpha*norm(x) < 0
        ar = norm(x,1)/norm(x);
    else
        al = alpha;
    end
    alpha = (ar+al)/2;
    if abs(ar-al) < 0.01
        break;
    end
    alpha_list = [alpha_list, alpha];
    obj_list = [obj_list, norm(x,1)-alpha*norm(x)];
end
result.i = i;
result.obj_list = obj_list;
result.alpha_list = alpha_list;
fx = norm(x,1)/norm(x);
result.fx = fx;
fxg = norm(xg,1)/norm(xg);
result.fxg = fxg;
% result.success = successR;
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
