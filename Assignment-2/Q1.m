clear; close all;
%% Load data
data = readmatrix('q1_data.xlsx');
X = data(:,2:5); % Time is not involved
y = data(:,end);
n = 31; p = 5;
%% Part a) - linear model
Xs = X - mean(X);
ys = y - mean(y);
beta_hat_class = (Xs'*Xs)\Xs'*y;
intercept = mean(y) - (mean(X))*beta_hat_class;
mdl_a = fitlm(X,y);
% Positively correlated with all except N2O
%% Part b) - outlier removal
error_var = sum((y-X*beta_hat_class-intercept).^2)/(n-p-1);
CI = [beta_hat_class - sqrt(error_var*inv(Xs'*Xs)), beta_hat_class + sqrt(error_var*inv(X'*X))];
res = y-X*beta_hat_class-intercept;
std_res = (res - mean(res))/sqrt(error_var);
plot(std_res,'x');hold on; yline(2,'r');yline(-2,'r');
ylim([-2.5,2.5]); title('Residual plot'); xlabel('Index'); ylabel('Residual');
% Iteration-1: remove 13th data point
indices = [1:1:31];
Xnew = X(indices(indices~=13),:);ynew = y(indices(indices~=13));
Xs = Xnew - mean(Xnew); beta_new = (Xs'*Xs)\Xs'*ynew;
beta0_new = mean(ynew) - (mean(Xnew))*beta_hat_class;
res = ynew-Xnew*beta_new-beta0_new;
std_res = (res - mean(res))/sqrt(error_var);
figure;
plot(std_res,'x');hold on; yline(2,'r');yline(-2,'r');
ylim([-2.5,2.5]); title('Residual plot'); xlabel('Index'); ylabel('Residual');
% Everything is within bounds now
%% Part c) - dropping insignificant variables
mdl_c = fitlm(Xnew,ynew)
% drop x4 (O3)
mdl_c1 = fitlm(Xnew,ynew,'y ~ 1 + x1 + x2 + x3')
% drop Intercept
mdl_c2 = fitlm(Xnew,ynew,'y ~ x1 + x2 + x3 - 1')
%% Part d)
CH4 = beta_hat_class(2)/beta_hat_class(1)*10^3;
N2O = beta_hat_class(3)/beta_hat_class(1)*10^3;