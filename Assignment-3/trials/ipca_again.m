clear; close all;
%% load data
load('flowdata3.mat');
n_samples = length(Fmeas);
n_var = size(Fmeas,2);
std_true = std(1:5);
%% IPCA
tol = 0.001;
val = tol;
vsmall = 1.0e-04;
Zs = Fmeas';
covZ = cov(Zs');
vlb = vsmall*sqrt(diag(covZ));  % Lower bound on std
sigma_e = vlb;
vub = sqrt(diag(covZ)); % Upper bound on std
% vub = [1,1,1,1,1]';
x0 = sqrt(vsmall)*vub;  % Initial estimate of error covariance taken as identify matrix
n_const = 3;
while val >= tol
    % PCA part
    Zs = Fmeas';
    for i = 1:n_var
        Zs(i,:) = Zs(i,:)./sqrt(n_samples)./sigma_e(i);
    end
    % Columns are samples, rows are variables. We want relationship across rows
    % so we will use the u matrix.
    [u,s,v] = svd(Zs,'econ');
    Amat = zeros(n_const,n_var);
    for i = 1:n_const
        Amat(i,:) = u(:,(n_var-n_const)+i)';
    end
    options = optimset('Display','iter','MaxFunEvals',50000);
    r = Amat*Fmeas';
    [sigma_e_est,fval,exitflag] = fmincon(@(X)(cost(Amat,n_samples,r,X)),x0,[],[],[],[],vlb,vub,[],options);
    exitflag
    x0 = sigma_e_est;
    val = sum(abs(sigma_e_est(n_var-n_const+1:end)-sigma_e(n_var-n_const+1:end)));
    sigma_e = sigma_e_est;
end
%% Objective to obtain sigma_e
function J = cost(A,N,r,sigma_e)
    matrix_e = diag(sigma_e);
    J = 0;
    for i = 1:N
        J = J + r(:,i)'*(inv(A*matrix_e*A'))*r(:,i);
    end
    J = J+ N*log(det(A*matrix_e*A'));
end