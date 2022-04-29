clc
clear all
X = xlsread('ghg-concentrations_1984-2014.xlsx','B6:E36');
y = xlsread('ghg-concentrations_1984-2014.xlsx','G6:G36');

% Apply Mutlilinear regression. Mean shift the data or add a column of ones to X.  Scaling the independent variables will not
% change the results

[N p] = size(X);
Xmean = mean(X);
ymean = mean(y);
Xs = X - ones(size(X))*diag(Xmean);
ys = y - ones(size(y))*ymean;
B = inv(Xs'*Xs)*Xs'*ys;
b0 = ymean - Xmean*B;

%  Part b
% Obtain prediction errors and estimate output error variance
err = ys - Xs*B;
dof = N - p -1;
sigmaerrest = sqrt(sum(err.^2)/dof);

% Compute confidence intervals for all regression slope coefficients
C = inv(X'*X);
tval = tinv(0.975,dof);
for i = 1:p
    CIL(i) = B(i) - sigmaerrest*tval*sqrt(C(i,i));
    CIU(i) = B(i) + sigmaerrest*tval*sqrt(C(i,i));
end

%  Compute standardized residuals to identify outliers
normerr = err/sigmaerrest;

%  Part c
Xssub = Xs(:,1:3);
Bsub = inv(Xssub'*Xssub)*Xssub'*ys;
b0sub = ymean - mean(Xssub)*Bsub;
