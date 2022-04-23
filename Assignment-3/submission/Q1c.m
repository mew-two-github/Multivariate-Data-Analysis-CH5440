%% Comments
% Same code as 1b with n_constraints alone changed from 3 to 4
clear; close all;
%% load data
load('flowdata3.mat');
Z = Fmeas'; 
n_samples = length(Fmeas);
n_var = size(Fmeas,2);
std_true = std(1:5);
%% IPCA
tol = 0.001;
val = tol;
vsmall = 1.0e-04;
covZ = cov(Fmeas);
vlb = vsmall*sqrt(diag(covZ));  % Lower bound on std
sigma_e = vlb;
vub = sqrt(diag(covZ)); % Upper bound on std
x0 = sqrt(vsmall)*vub;  % Initial estimate of error covariance taken as identify matrix
n_const = 4;
while val >= tol
    % PCA part
    for i = 1:n_var
        Zs(i,:) = Z(i,:)./sigma_e(i)/sqrt(n_samples);
    end
    % Columns are samples, rows are variables. We want relationship across rows
    % so we will use the u matrix.
    [u,s,v] = svd(Zs,'econ');
    Amat = zeros(n_const,n_var);
    for i = 1:n_const
        Amat(i,:) = u(:,(n_var-n_const)+i)';
    end
    for k = 1:n_var
        Amat(:,k) = Amat(:,k)./sigma_e(k)/sqrt(n_samples);
    end
    options = optimset('Display','iter','MaxFunEvals',50000);
    [sigma_e_est,fval,exitflag] = fmincon('cost_fn',x0,[],[],[],[],vlb,vub,[],options,Amat,Z);
    x0 = sigma_e_est;
    val = (abs(sum(sigma_e_est(n_var-n_const+1:end))-sum(sigma_e(n_var-n_const+1:end))));
    sigma_e = sigma_e_est;
end
%% Examining the obtained matrix
for i = 1:n_var
    Zs(i,:) = Z(i,:)./sigma_e(i)/sqrt(n_samples);
end
eig_values = eig(Zs*Zs')