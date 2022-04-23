clear; close all;
%% load data
load 'Inorfull';

Yavg = zeros(26,176);
stdavg = zeros(26,176);
concavg = zeros(26,3);

for nPC=1:26
    sindex = 5*(nPC-1)+1;
    eindex = sindex+4;
    Yavg(nPC,:) = mean(DATA(sindex:eindex,:),1);
    stdavg(nPC,:) = mean(stdDATA(sindex:eindex,:),1)/sqrt(5);
    concavg(nPC,:) = mean(CONC(sindex:eindex,:),1);
end
%% IPCR - IPCA step
tol = 0.0351;
val = tol;
vsmall = 1.0e-04;
Z = Yavg';
[n_var n_samples] = size(Z);
covZ = cov(Yavg);
vlb = vsmall*sqrt(diag(covZ));  % Lower bound on std
sigma_e = mean(stdavg)'; % Initial estimate as obtained from replicates
vub = 2*sqrt(diag(covZ)); % Upper bound on std; 2* because initial estimate
% not coming inside vub
x0 = sqrt(vsmall)*vub;  
n_const = 173; % since 176 variables and 3 concentrations
Zs = zeros(size(Z));
while val >= tol
    % PCA part
    for i = 1:n_var
        Zs(i,:) = Z(i,:)./sigma_e(i)/sqrt(n_samples);
    end
    % Columns are samples, rows are variables. We want relationship across rows
    % so we will use the u matrix.
    [u,s,~] = svd(Zs);
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
%% IPCR - OLS step
plot(1:1:176,sigma_e);
xlabel('Wavelengths');ylabel('$\hat{\sigma}$','interpreter','latex');
Yavgs = Yavg*inv(diag(sigma_e));
[~,~,v] = svd(Yavgs,'econ');
T = Yavgs*v(:,1:3);
B = (T'*T)\T'*concavg;
errors = (T*B-concavg)/size(Yavg,1);
RMSE = norm(errors(:));