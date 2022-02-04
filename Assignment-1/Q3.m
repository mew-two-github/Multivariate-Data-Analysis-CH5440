close all; clear;
%% Data
EP = [1.98,2.31,3.29,3.56,1.23,1.57,2.05,0.66,0.31,2.82,0.13,3.15,2.72,2.31,1.92,1.56,0.94,2.27,3.17,2.36]';
CF = [1.87,2.2,3.15,3.42,1.1,1.41,1.84,0.68,0.27,2.8,0.14,3.2,2.7,2.43,1.78,1.53,0.84,2.21,3.10,2.34]';
% u - EP, y - CF
u = EP;
y = CF;
N = length(u);
%% Part a)
% OLS
[alpha_OLS, beta_OLS, uhat_OLS, yhat_OLS,s_OLS] = OLS(u,y);
sigma_e_OLS = 1/(N-2)*sum((y-alpha_OLS*u-beta_OLS).^2);
CI_OLS = [alpha_OLS-2.16*sigma_e_OLS,alpha_OLS+2.16*sigma_e_OLS];
% IOLS
[alpha_IOLS, beta_IOLS, uhat_IOLS, yhat_IOLS,s_IOLS] = IOLS(u,y);
sigma_e_IOLS = 1/(N-2)*sum((y-alpha_IOLS*u-beta_IOLS).^2);
CI_IOLS = [alpha_IOLS-2.16*sigma_e_IOLS,alpha_OLS+2.16*sigma_e_IOLS];
% TLS
[alpha_TLS, beta_TLS, uhat_TLS, yhat_TLS,s_TLS] = TLS(u,y);
sigma_e_TLS = 1/(N-2)*sum((y-alpha_TLS*u-beta_TLS).^2);
CI_TLS = [alpha_TLS-2.16*sigma_e_TLS,alpha_TLS+2.16*sigma_e_TLS];
%% Part b)
ui = 2.31;
yi = 2.20;
% OLS: u is perfect
OLS_pred = ui;
% IOLS: y is perfect
IOLS_pred = yi;
% TLS: both imperfect. doing a perpendicular projection will give us an
% estimate for EP and CF. Since we need a single estimate we will assume
% alpha = 1 and beta = 0. So yhat = uhat = (u+y)/2
TLS_pred = (ui+yi)/2;