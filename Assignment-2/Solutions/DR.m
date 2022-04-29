function [Test, Pest, err] = DR(T, Psat, stdT, stdP, A, B, C);
%
%  Given parameters A, B, C and measured values of T and Psat determine the
%  reconciled values of Temp and Psat in the WTLS sense subject to the
%  model constraint 
%
