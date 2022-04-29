function f = myDRfun(x, T, Psat, stdT, stdP, A, B, C)
%
%  Computes the WTLS function value for given estimate of x and measured
%  values of T and Psat and parameters A, B, C
%  x is estimated value of T 
%
Pest = exp(A - B/(x + C)); % Estimate P using the model and given estimate of T
f = ((T - x)/stdT)^2 + ((Psat - Pest)/stdP)^2;
