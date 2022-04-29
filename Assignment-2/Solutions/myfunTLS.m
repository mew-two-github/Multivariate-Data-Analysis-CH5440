function err = myfunTLS(x,T,Psat,stderr)
%
% x is a decision vector with first element A, second is B and third is C.  y
% is the actual measured outputs
%
n = length(T);
% 
%  Find the TLS estimates of T and Psat for each sample for the given value
%  of parameters
%
err = zeros(n,1);
options = optimoptions('fminunc','Display','off');
for i = 1:n
    x0 = T(i);
    [Test, fval] = fminunc(@(y) myDR(y, T(i), Psat(i), stderr(1), stderr(2), x(1), x(2), x(3)), x0, options);
    err(i) = sqrt(fval);
end
