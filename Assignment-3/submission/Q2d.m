clear; close all;
%% load data
load 'Inorfull';

Yavg = zeros(26,176);
stdavg = zeros(26,176);
concavg = zeros(26,3);

for nPC=1:26
    sindex = 5*(nPC-1)+1;
    eindex = sindex+4;
    Yavg(nPC,:) = mean(DATA(sindex:eindex,:),1);
    stdavg(nPC,:) = mean(stdDATA(sindex:eindex,:),1)/sqrt(5);
    concavg(nPC,:) = mean(CONC(sindex:eindex,:),1);
end
[n_samples, n_var] = size(Yavg);
%% Applying MLPCR
% adding constant to have better covariance matrix
% stdavg = stdavg + 1;
[vecs,scores] = MLPCA(Yavg,stdavg,3);
B = (scores'*scores)\scores'*concavg;
errors = (scores*B-concavg)/size(Yavg,1);
RMSE = sqrt(sum((errors.*errors),'all')/n_samples);
%% LOOCV for MLPCR
avg_RMSE = zeros(10,3);
for nPC = 1:10
    avg_RMSE(nPC,:) = LOOCV_MLPCR(Yavg,concavg,stdavg,nPC);
end
PRESS = sum(avg_RMSE,2);
plot(PRESS);
%% function to find LOOCV for MLPCR
function RMSE = LOOCV_MLPCR(Yavg,concavg,sigma,p)
    [nsamples,nvar] = size(concavg);
    RMSE = zeros(1,nvar);
    for i = 1:nsamples
        Ysub = [Yavg(1:i-1,:); Yavg(i+1:end,:)];
        concsub = [concavg(1:i-1,:); concavg(i+1:end,:)];
        sigmasub = [sigma(1:i-1,:); sigma(i+1:end,:)];
        [vecs,Tsub] = MLPCA(Ysub,sigmasub,p);
        B = (Tsub'*Tsub)\Tsub'*concsub;
        covZ = diag(sigma(i,:).^2);
        prederr = concavg(i,:) - Yavg(i,:)*inv(covZ)*vecs*inv(vecs'*inv(covZ)*vecs)*B;
        RMSE = RMSE + prederr.*prederr;
    end
    RMSE = sqrt(RMSE/nsamples);
end
    