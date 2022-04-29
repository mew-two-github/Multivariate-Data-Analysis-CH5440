function err = myfun(x,T,psat)
%
% x is a decisionvector with first elemnt A, second is B and third is C.  y
% is the actual measured outputs
%
n = length(T);
psatest = zeros(n,1);
for i = 1:n
    psatest(i) = exp(x(1) - x(2)/(T(i)+x(3)));
end
err = psatest - psat;