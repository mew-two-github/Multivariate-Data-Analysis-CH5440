function [RMSE] = LOOCV_PCR(X,Y,nfact)
%
%  Function for performing leave one sample out cross validation using PCR
%
%  X : N x n matrix of inputs where N is number of samples and n number of
%  variables
%  Y : N x p output vector (assumed to be only one output)
[nsamples nvar] = size(Y);
%
% Estimate OLS regression matrix leaving one sample out in turn 
%
RMSE = zeros(1,nvar);
%
%  Build PCR model dropping each sample out in turn
for i = 1:nsamples
    Xsub = [X(1:i-1,:); X(i+1:end,:)];
    Ysub = [Y(1:i-1,:); Y(i+1:end,:)];
    [u s v] = svd(Xsub,'econ');
    Tsub = Xsub*v(:,1:nfact);
    B = inv(Tsub'*Tsub)*Tsub'*Ysub;
    prederr = Y(i,:) - X(i,:)*v(:,1:nfact)*B;
    RMSE = RMSE + prederr.*prederr;
end
RMSE = sqrt(RMSE/nsamples);
    
    
    
