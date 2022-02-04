clear; close all;
%% Open data
A = readmatrix('defects_annotation_data.csv');
x1 = rem_NaN(A(:,1)); y1 = rem_NaN(A(:,2));
x2 = rem_NaN(A(:,4)); y2 = rem_NaN(A(:,5));
x3 = rem_NaN(A(:,7)); y3 = rem_NaN(A(:,8));
%% Defect-1
N = length(x1);
% TLS
[alpha_TLS, beta_TLS, uhat_TLS, yhat_TLS,s_TLS] = TLS(x1,y1);
sigma_e_TLS = 1/(N-2)*sum((y1-alpha_TLS*x1-beta_TLS).^2);
CI_TLS = [alpha_TLS-2.16*sigma_e_TLS,alpha_TLS+2.16*sigma_e_TLS];
% TLS - inverted
[alpha_TLS2, beta_TLS2, uhat_TLS2, yhat_TLS2,s_TLS2] = TLS(y1,x1);
sigma_e_TLS2 = 1/(N-2)*sum((y1-alpha_TLS2*x1-beta_TLS2).^2);
CI_TLS2 = [alpha_TLS2-2.16*sigma_e_TLS2,alpha_TLS2+2.16*sigma_e_TLS2];
% OLS
[alpha_OLS, beta_OLS, uhat_OLS, yhat_OLS,s_OLS] = OLS(x1,y1);
sigma_e_OLS = 1/(N-2)*sum((y1-alpha_OLS*x1-beta_OLS).^2);
CI_OLS = [alpha_OLS-2.16*sigma_e_OLS,alpha_OLS+2.16*sigma_e_OLS];
% OLS - inverted
[alpha_OLS2, beta_OLS2, uhat_OLS2, yhat_OLS2,s_OLS2] = OLS(y1,x1);
sigma_e_OLS2 = 1/(N-2)*sum((y1-alpha_OLS2*x1-beta_OLS2).^2);
CI_OLS2 = [alpha_OLS2-2.16*sigma_e_OLS2,alpha_OLS2+2.16*sigma_e_OLS2];
subplot(2,2,1);
plot(x1,y1,'x',x1,x1*alpha_OLS+beta_OLS,y1*alpha_OLS2+beta_OLS2,y1);
title('Defect-1'); xlabel('x-coordinate'); ylabel('y-coordinate');
legend('Data','OLS','IOLS');
%% Defect-2
% TLS
[alpha_TLS_def2, beta_TLS_def2, uhat_TLS_def2, yhat_TLS_def2,s_TLS_def2] = TLS(x2,y2);
sigma_e_TLS_def2 = 1/(N-2)*sum((y2-alpha_TLS_def2*x2-beta_TLS_def2).^2);
CI_TLS_def2 = [alpha_TLS_def2-2.16*sigma_e_TLS_def2,alpha_TLS_def2+2.16*sigma_e_TLS_def2];
% TLS - inverted
[alpha_TLS2_def2, beta_TLS2_def2, uhat_TLS2_def2, yhat_TLS2_def2,s_TLS2_def2] = TLS(y2,x2);
sigma_e_TLS2_def2 = 1/(N-2)*sum((y2-alpha_TLS2_def2*x2-beta_TLS2_def2).^2);
CI_TLS2_def2 = [alpha_TLS2_def2-2.16*sigma_e_TLS2_def2,alpha_TLS2_def2+2.16*sigma_e_TLS2_def2];
% OLS
[alpha_OLS_def2, beta_OLS_def2, uhat_OLS_def2, yhat_OLS_def2,s_OLS_def2] = OLS(x2,y2);
sigma_e_OLS_def2 = 1/(N-2)*sum((y2-alpha_OLS_def2*x2-beta_OLS_def2).^2);
CI_OLS_def2 = [alpha_OLS_def2-2.16*sigma_e_OLS_def2,alpha_OLS_def2+2.16*sigma_e_OLS_def2];
% OLS - inverted
[alpha_OLS2_def2, beta_OLS2_def2, uhat_OLS2_def2, yhat_OLS2_def2,s_OLS2_def2] = OLS(y2,x2);
sigma_e_OLS2_def2 = 1/(N-2)*sum((y2-alpha_OLS2_def2*x2-beta_OLS2_def2).^2);
CI_OLS2_def2 = [alpha_OLS2_def2-2.16*sigma_e_OLS2_def2,alpha_OLS2_def2+2.16*sigma_e_OLS2_def2];
subplot(2,2,2);
plot(x2,y2,'x',x2,x2*alpha_OLS_def2+beta_OLS_def2,y2*alpha_OLS2_def2+beta_OLS2_def2,y2);
title('Defect-2'); xlabel('x-coordinate'); ylabel('y-coordinate');
legend('Data','OLS','IOLS');
%% Defect-3
% TLS
[alpha_TLS_def3, beta_TLS_def3, uhat_TLS_def3, yhat_TLS_def3,s_TLS_def3] = TLS(x3,y3);
sigma_e_TLS_def3 = 1/(N-2)*sum((y3-alpha_TLS_def3*x3-beta_TLS_def3).^2);
CI_TLS_def3 = [alpha_TLS_def3-2.16*sigma_e_TLS_def3,alpha_TLS_def3+2.16*sigma_e_TLS_def3];
% TLS - inverted
[alpha_TLS2_def3, beta_TLS2_def3, uhat_TLS2_def3, yhat_TLS2_def3,s_TLS2_def3] = TLS(y3,x3);
sigma_e_TLS2_def3 = 1/(N-2)*sum((y3-alpha_TLS2_def3*x3-beta_TLS2_def3).^2);
CI_TLS2_def3 = [alpha_TLS2_def3-2.16*sigma_e_TLS2_def3,alpha_TLS2_def3+2.16*sigma_e_TLS2_def3];
% OLS
[alpha_OLS_def3, beta_OLS_def3, uhat_OLS_def3, yhat_OLS_def3,s_OLS_def3] = OLS(x3,y3);
sigma_e_OLS_def3 = 1/(N-2)*sum((y3-alpha_OLS_def3*x3-beta_OLS_def3).^2);
CI_OLS_def3 = [alpha_OLS_def3-2.16*sigma_e_OLS_def3,alpha_OLS_def3+2.16*sigma_e_OLS_def3];
% OLS - inverted
[alpha_OLS2_def3, beta_OLS2_def3, uhat_OLS2_def3, yhat_OLS2_def3,s_OLS2_def3] = OLS(y3,x3);
sigma_e_OLS2_def3 = 1/(N-2)*sum((y2-alpha_OLS2_def3*x2-beta_OLS2_def3).^2);
CI_OLS2_def3 = [alpha_OLS2_def3-2.16*sigma_e_OLS2_def3,alpha_OLS2_def3+2.16*sigma_e_OLS2_def3];
subplot(2,2,3);
plot(x3,y3,'x',x3,x3*alpha_OLS_def3+beta_OLS_def3,y3*alpha_OLS2_def3+beta_OLS2_def3,y3);
title('Defect-3'); xlabel('x-coordinate'); ylabel('y-coordinate');
legend('Data','OLS','IOLS');
%% function to remove NaN values
function vnew = rem_NaN(v)
    vnew = v(~isnan(v));
end