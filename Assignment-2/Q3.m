close all; clear;
%% Data
xbar = [9;68;129];
S = [7 21 34;21 64 102;34 102 186];
%% Part a
[v,d] = eig(S);
ev = diag(d);
%% Part b
val = ev(3)/sum(ev)
% one eigen value
%% Part c
% corresponding to smaller two evs
relations = v(:,1:2);
%% Part d
x = [10.1;73;135.5];
scores = (x-xbar)'*v;
%% Part e 
z_hats = [0.233 -0.5135;0.9589 -0.02]\[-7.99-0.8258*73;-13.2+0.2831*73] ;
mass_e = z_hats(1);
%% Part f

