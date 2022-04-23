% Test PCA
clc
clear all

nvar = 5;
nsamples = 1000;
% True constraint matrix
Atrue = [1 1 -1 0 0; 0 0 1 -1 0; 0 -1 0 1 -1];
% True values of variables
Fbase(1,1) = 10;
Fbase(1,2) = 10;
Fbase(1,3) = Fbase(1,1)+Fbase(1,2);
Fbase(1,4) = Fbase(1,3);
Fbase(1,5) = Fbase(1,4) - Fbase(1,2);

% Parameters for random walk model
lamda(1) = 1.0;
lamda(2) = 2.0;
lamda(3) = 0.8;
lamda(4) = 1.2;


% Standard deviation of measurement errors (No errors)
% std(1) = 0.0;
% std(2) = 0.0;
% std(3) = 0.0;
% std(4) = 0.0;
% std(5) = 0.0;

% Standard deviation of measurement errors (equal variances)
% std(1) = 0.1;
% std(2) = 0.1;
% std(3) = 0.1;
% std(4) = 0.1;
% std(5) = 0.1;
% 
% Standard deviation of measurement errors (unequal variances - high SNR)
std(1) = 0.1;
std(2) = 0.08;
std(3) = 0.15;
std(4) = 0.2;
std(5) = 0.18;
% 
% Standard deviation of measurement errors (unequal variances - low SNR)
% std(1) = 1;
% std(2) = 0.8;
% std(3) = 1.5;
% std(4) = 2.0;
% std(5) = 1.8;
% 
Ltrue = diag(std);

% Tolerance for convergence
tol = 1e-06;
% Reset random number generator to initial state
randn('state',5);

for j = 1:nsamples
    for k = 1:2
        Ftrue(j,k) = Fbase(k) + lamda(k)*randn; 
    end
    Ftrue(j,3) = Ftrue(j,1) + Ftrue(j,2);
    Ftrue(j,4) = Ftrue(j,3);
    Ftrue(j,5) = Ftrue(j,4) - Ftrue(j,2);
end

for j = 1:nsamples
    error = randn(nvar,1);
    Fmeas(j,:) = Ftrue(j,:) + (Ltrue*error)';
end

SY = zeros(nvar);
for j = 1:nsamples
    SY = SY + Fmeas(j,:)'*Fmeas(j,:);
end
meanY = mean(Fmeas);
covY = (SY - nsamples*meanY'*meanY)/(nsamples-1);

% Different scaling strategies
% Lsinv = eye(nvar);  % No scaling
% Lsinv = diag(ones(nvar,1)./sqrt(diag(covY)));  % Scaling by measurement variances
Lsinv = inv(Ltrue);  % Scaling by cholesky factor of true error covariance matrix

nfact = 2;
%  Determine how good the constraint matrix has been estimated by PCA
Ys =Lsinv*Fmeas'/sqrt(nsamples);

[u s v] = svds(Ys,nvar);
for k = nfact+1:nvar
    Amat(k-nfact,:) = u(:,k)';
end
Amat = Amat*Lsinv;
sval = diag(s);
% Amat
spca = diag(s);
theta_pca = 180*subspace(Atrue', Amat')/pi;

% Determine how well the model matches with the true constraint matrix.
% For this determine the minimum distance of each true constraint vector from the
% row space of model constraints
for i = 1:3
    bcol = Atrue(i,:)';
    dist_pca(i) = norm(bcol - Amat'*inv(Amat*Amat')*Amat*bcol);
end

dist_pca
sum(dist_pca)
theta_pca
sval
