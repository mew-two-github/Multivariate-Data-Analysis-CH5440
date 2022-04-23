clear; close all;
%% load data
load('flowdata3.mat');
n_samples = length(Fmeas);
n_var = size(Fmeas,2);
std_true = std(1:5);
%% PCA
Ltrue = diag(std_true);
Lsinv = inv(Ltrue);
Zs = Lsinv*Fmeas'/sqrt(n_samples);
% Columns are samples, rows are variables. We want relationship across rows
% so we will use the u matrix.
[u s v] = svd(Zs,'econ');
n_const = 3;
Amat = zeros(n_const,n_var);
for i = 1:n_const
    Amat(i,:) = u(:,(n_var-n_const)+i)';
end
Amat_actual = Amat*Lsinv;
Ad = Amat_actual(:,[1,2,4]);
Ai = Amat_actual(:,[3,5]);
Areg = -1*Ad\Ai;
Atrue_reg = -1*Atrue(:,[1,2,4])\Atrue(:,[3,5]);
% Maximum absolute difference between the constraint matrices
maxdiff = max(abs(Areg-Atrue_reg),[],'all');
eig_values = eig(Zs*Zs');
