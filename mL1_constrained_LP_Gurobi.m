function xr = mL1_constrained_LP_Gurobi(A,b)
%=============================================================
% gurobi for the L1 constrained model 
%
% Solves
%           min  norm(x,1)
%           s.t. Ax = b
%
% 
% Author: Chao Wang  
% Date: Nov 11. 2018
%============================================================= 

[~, N] = size(A);
f = ones(2*N,1);
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
if isfield(results,'x')
   x = results.x();
else
   fprintf(results.status);
   fprintf(',No results in L1... \n');
   x=randn(2*N,1);
end
xr = x(1:N) - x(N+1:2*N);

end

