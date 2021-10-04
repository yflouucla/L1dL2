function [xg] = data_generator_xg(pm)
%=============================================================
% DESCRIPTION: the option of x and its parameter
%             'dynamic': pm.range_factor
%             'normalized' pm.normalized
%             'min_speration': pm.min_speration
%             'sparsity': pm.sparsity



N = pm.N; % Size of xg
sparsity = pm.sparsity; % Sparsity

if isfield(pm,'min_seperation'); min_seperation = pm.min_seperation;
else;min_seperation = 1; end  % Minimum seperation, default as 1;

if isfield(pm,'normalized'); normalized = pm.normalized;
else; normalized = 0; end % Normalized xg, i.e. max(xg) = 1
if isfield(pm,'dynamic'); dynamic = pm.dynamic; else; dynamic = 0;end
% building the vector xg
xg = zeros(N,1);
supp  = randsample_separated(N,sparsity,min_seperation);

if dynamic == 1
    %--------------- Dynamic ---------------------%
    if isfield(pm,'range_factor'); range_factor = pm.range_factor;
    else;range_factor = 5; end
    xs = sign(randn(sparsity,1)).*10.^(range_factor*rand(sparsity,1));
else
    %----------------- uniform ---------------------%
    xs = randn(sparsity,1);
end

xg(supp) = xs;

if normalized == 1 % normalizing xg if wanted
    xg = xg/max(abs(xg));
end

end

function supp = randsample_separated(N,K,L)
% random sampling K integers from 1--N with spacing at least L
supp = randsample(N-L*(K),K);
supp = sort(supp);
supp = supp + (0:K-1)'*L;
end