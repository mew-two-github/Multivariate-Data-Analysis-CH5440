clear; close all;
%% Open data
A = readmatrix('CO2.csv');
time = A(:,1);
CO2 = A(:,2);
temp = A(:,3);
temp_cut = 3.6;
%% Predict maximum permissible level of CO2 - OLS
% y -> temperature deviation
% u -> CO2
u = CO2;
y = temp;
[alpha, beta, uhat, yhat, s] = OLS(u,y);
disp(s);
CO2_cut_OLS = (temp_cut-beta)/alpha;
%% Predict maximum permissible level of CO2 - TLS
[alpha_TLS, beta_TLS, uhat_TLS, yhat_TLS,s_TLS] = TLS(u,y);
CO2_cut_TLS = (temp_cut-beta_TLS)/alpha_TLS;
%% Predict year
ut = time;
yt = CO2;
[alpha_t, beta_t, uhat_t, yhat_t, s_t] = OLS(ut,yt);
t_OLS = (CO2_cut_OLS-beta_t)/alpha_t;
tpred_OLS = ceil(t_OLS);
t_TLS = (CO2_cut_TLS-beta_t)/alpha_t;
tpred_TLS = ceil(t_TLS);

