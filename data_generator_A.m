function [A] = data_generator_A(pm)
%=============================================================
% DESCRIPTION: the type of matrix defined by mode and the parameter
%             'Gaussian': pm.r
%             'Oversampled_DCT': pm.F
%             'Topelitz': pm.rho
%             'DCT'

M = pm.M; N = pm.N;% dimention of matrix is M by N
if isfield(pm,'sen_mat'); sen_mat = pm.sen_mat; end


if strcmp(sen_mat,'Gaussian')
    %--------------- Gaussian ---------------------%
    % r is between 0 and 1
    if isfield(pm,'r'); r = pm.r; else; r = 0.2; end
    sigma = r*ones(N);
    sigma = sigma + (1-r)*eye(N);
    mu = zeros(1,N);
    A = mvnrnd(mu,sigma,M);
    
elseif strcmp(sen_mat,'Oversampled_DCT')
    %--------------- oversampled DCT ---------------------%
    if isfield(pm,'F'); F = pm.F; else; F = 10; end
    w = rand(M,1); W_col = w*ones(1,N); k_row = ones(M,1)*(1:N);
    A = cos(2*pi*W_col.*k_row/F)/sqrt(M);
    
elseif strcmp(sen_mat,'DCT')
    %----------------- DCT matrix  -------------------------%
    A = dctmtx(N);
    tmp = randperm(N-1);
    A = A([1 tmp(1:M-1)+1],:); % randomly select m rows but always include row 1
    A = A/norm(A);
else
    fprintf('This mode does not exists... \n')
end



end
