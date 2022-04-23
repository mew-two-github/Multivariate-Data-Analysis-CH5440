function [RMSE] = LOOCV_PCR(X,Y,nfact)
[nsamples,nvar] = size(Y);

RMSE = zeros(1,nvar);

for i = 1:nsamples
    Xsub = [X(1:i-1,:); X(i+1:end,:)];
    Ysub = [Y(1:i-1,:); Y(i+1:end,:)];
    [~,~,v] = svd(Xsub,'econ');
    Tsub = Xsub*v(:,1:nfact);
    B = (Tsub'*Tsub)\Tsub'*Ysub;
    prederr = Y(i,:) - X(i,:)*v(:,1:nfact)*B;
    RMSE = RMSE + prederr.*prederr;
end
RMSE = sqrt(RMSE/nsamples);
    
    
    
