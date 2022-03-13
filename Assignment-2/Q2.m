clear; close all;
%% Load data
load('vpdata.mat');
N = length(temp);
%% Part a) - Linear regression
X = 1./temp;
y = log(psat);
[alpha, beta, uhat, yhat, s] = OLS(X,y);
% Confirming answer with MATLAB in-built function
mdl_a = fitlm(X,y);
%% Part b) - non-linear regression
antoine = @(X)(X(1) - X(2)./(temp + X(3)));
% Vector valued function for lsqnonlin
f = @(X)(psat - exp(antoine(X)));
% Ensuring proper order of the initial guess
X0 = [1,1000,100];
[Params_b,resnorm,res,exitflag] = lsqnonlin(f,X0);%,[],[],opts);
figure;
plot(temp,antoine(Params_b),'o',temp,y,'x');
ylabel('log(P)'); xlabel('Temperature');
errors_b = f(Params_b);
error_var = var(errors_b);
%% Part c) - WTLS
antoine_c = @(X)(X(1) - X(2)./(X(4:end) + X(3)));
% Given measurement errors
sigma_x = 0.18;
sigma_y = 2;
% Matrix valued function for lsqnonlin - it minimises sum of squares of
% each element in the matrix
obj_c = @(X)([(psat - exp(antoine_c(X)))/sigma_y,(temp-X(4:end))/sigma_x]);
X0 = Params_b';
X0(4:length(temp)+3) = temp;
[Params_c,resnorm_c,res_c,exitflag_c] = lsqnonlin(obj_c,X0);
%% Part d)
max1 = max(abs(psat - exp(alpha*X + beta)));
max2 = max(abs(res));
max3 = max(abs(res_c(:,1))*sigma_y);
%% function to perform OLS
function [alpha, beta, uhat, yhat, s] = OLS(u,y)
    N = length(u);
    ybar = mean(y);
    ubar = mean(u);
    syy = var(y,1);
    suu = var(u,1);
    syu = 1/N*sum((y-ybar).*(u-ubar));
    uhat = u;
    alpha = syu/suu;
    beta = ybar - alpha*ubar;
    yhat = alpha*uhat + beta;
    s = struct('alpha',alpha,'beta',beta,'ybar',ybar,'ubar',ubar,'syy',syy,'syu',syu,'suu',suu);
end