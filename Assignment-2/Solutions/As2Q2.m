clc
clear all
load ../vpdata

% Part a

x  = 1./temp;
y = log(psat);
Smat = cov(x,y);
alpha = Smat(1,2)/Smat(1,1);
beta = mean(y) - alpha*mean(x);
psatest = exp(alpha*x + beta);
err = psatest - psat;
max(abs(err))

% Part b

%  Initial values od decision variables based on results of part a
%
x0(1) = beta;
x0(2) = -alpha;
x0(3) = 0;
options = optimoptions('lsqnonlin','MaxIterations',1000,'MaxfunctionEvaluations',10000,'Display','iter');
xstar = lsqnonlin(@(x) myfun(x,temp,psat),x0,[],[],options);
max(abs(myfun(xstar,temp,psat)))

% Part c
%  Initial values of A, B and C from part b
x0(1) = beta; % xstar(1);
x0(2) = -alpha; % xstar(2);
x0(3) = 0; % xstar(3);
stderr(1) = 0.18; % Std of error in T
stderr(2) = 2; % Std of error in Psat
options = optimoptions('lsqnonlin','MaxIterations',1000,'MaxfunctionEvaluations',10000,'Display','iter');
xstar = lsqnonlin(@(x) myfunTLS(x,temp,psat,stderr),x0,[],[],options);

% Estimate temp for each sample based on optimal values of parameters A, B,
% and C and use it to compute Psat estimate

n = length(temp);
options = optimoptions('fminunc','Display','off');
Test = zeros(n,1);
psatest = zeros(n,1);
for i = 1:n
    t0 = temp(i);
    [Test(i), fval] = fminunc(@(test) myDR(test, temp(i), psat(i), stderr(1), stderr(2), xstar(1), xstar(2), xstar(3)), t0, options);
    psatest(i) = exp(xstar(1) - xstar(2)/(Test(i) + xstar(3)));
end

max(abs(psatest - psat));


