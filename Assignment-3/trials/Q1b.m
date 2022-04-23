clear; close all;
%% load data
load('flowdata3.mat');
n_samples = length(Fmeas);
n_var = size(Fmeas,2);
std_true = std(1:5);
%% IPCA
sigma_e = eye(5);
tol = 0.1;
val = tol;
while val >= tol
    % PCA part
    Zs = Fmeas';
    for i = 1:n_var
        Zs(i,:) = Zs(i,:)./sigma_e(i,i)/sqrt(n_samples);
    end
    % Columns are samples, rows are variables. We want relationship across rows
    % so we will use the u matrix.
    [u,s,v] = svd(Zs,'econ');
    n_const = 3;
    Amat = zeros(n_const,n_var);
    for i = 1:n_const
        Amat(i,:) = u(:,(n_var-n_const)+i)';
    end
    % Amat_actual = Amat*Lsinv;
    % Estimating sigma^2
    % r = Amat_actual*Zs;
    % r = Amat*Zs;
    % Sr = (r*r')/n_samples;
    % Akron = kron(Amat_actual,Amat_actual);
    % Akron = kron(Amat,Amat);
    % sigma_e_est = (Akron'*Akron)\Akron'*Sr(:);
    % cost(Amat,n_samples,r,X)
    [sigma_e_est,fval,exitflag] = fmincon(@(X)(obj_val(X,Amat,Zs)),diag(sigma_e),-1*eye(5),zeros(5,1));
    exitflag
    val = norm(sigma_e_est-diag(sigma_e)');
%     sigma_e = reshape(sigma_e_est,[5,5]);
    for i = 1:n_var
        sigma_e(i,i) = sigma_e_est(i);
    end
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
% function J = cost(Akron,Sr,sigma_e)
%     mat = diag(sigma_e);
%     J = norm(Sr(:) - Akron*mat(:));
% end