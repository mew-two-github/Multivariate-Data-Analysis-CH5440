clc
clear all

% Part a

zmean = [9 68 129];
S = [7 21 34; 21 64 102; 34 102 186];
[V D] = eig(S);

% Part b

lambda = diag(D);
lambda(3)/sum(lambda)

% Part c

b1 = V(:,1)'*zmean';
b2 = V(:,2)'*zmean';

% Part d

score = ([10.1 73 135.5] - zmean)*V(:,3);
