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
%% Perform PCR with LOOCV for different number of PCs
avg_RMSE = zeros(10,3);
% PCR on averaged values
for nPC = 1:10
    avg_RMSE(nPC,:) = LOOCV_PCR(Yavg,concavg,nPC);
end
PRESS = sum(avg_RMSE,2);
plot(PRESS)
% take 3 components according to PRESS
[~,~,v] = svd(Yavg,'econ');
T = Yavg*v(:,1:3);
B = (T'*T)\T'*concavg;
errors = (T*B-concavg)/size(Yavg,1);
RMSE = norm(errors(:));
%% Part b) - scaled PCR
% Plot the average standard deviations
figure;
plot(stdavg');
title('Plot of average standard deviations');
ylabel('\sigma');xlabel('index');
% Assume variation only along wavelength and not along sample
stds = mean(stdavg);
% Scale accordingly and perform PCR
L = diag(stds);
Yavgs = Yavg*inv(L);
avg_RMSE_scaled = zeros(10,3);
for nPC = 1:10
    avg_RMSE_scaled(nPC,:) = LOOCV_PCR(Yavgs,concavg,nPC);
end
PRESS_scaled = sum(avg_RMSE_scaled,2);
figure;
plot(PRESS_scaled);
% Scaled PCR giving us better results in terms of PRESS. And same as in
% previous case, nPC = 3 is the best pick according to PRESS
% take 3 components according to PRESS
[~,~,v2] = svd(Yavgs,'econ');
T2 = Yavgs*v2(:,1:3);
B2 = (T2'*T2)\T2'*concavg;
% Need to preprocess with same sigma^2 to use this B
errors_scaled = (T2*B2-concavg)/size(Yavg,1);
RMSE_scaled = norm(errors_scaled(:));