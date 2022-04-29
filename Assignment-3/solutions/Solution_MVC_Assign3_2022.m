%  Test the expt. data sets downloaded from Wentzell's site both with MLPCA
%  and IPCA (assuming covariance matrix of errors for all samples are
%  same).  Here the data at different wavelengths is assumed to be the
%  different samples.  Absorbance for 26 different mixture concentrations at 176
%  different wavelengths are taken and each is repeated 5 times to give a
%  data matrix of size 130 x 176 (n x N). Corresponding to each conc. since
%  5 readings are taken at all wavelengths, we compute the average of every
%  5 rows to get a data matrix of size 26 x 176.
clear all
clc

load 'Inorfull'

%  Compute average of every 5 rows of DATA matrix and store.  To be used
%  with lab data generated experimentally by Wentzell and loaded in his web site.

Yavg = [];
Ysample = [];
stdavg = [];
concavg = [];
for i=1:26
    istart = 5*(i-1)+1;
    Ysample = [Ysample; DATA(istart,:)];
    iend = istart+4;
    avg = mean(DATA(istart:iend,:),1);
    stda = mean(stdDATA(istart:iend,:),1)/sqrt(5);
    Yavg = [Yavg; avg];
    stdavg = [stdavg; stda];
    avg = mean(CONC(istart:iend,:),1);
    concavg = [concavg;avg];
end


% -------------------Multilinear OLS Regression--------------------
% Data for multilinear regression
% Maximum absorbances for Ni, Cr, Co corresponds to array indices 12, 54 and 106 
% YsampleMLR = [Ysample(:,48) Ysample(:,54) Ysample(:,106)];
% YavgMLR = [Yavg(:,48) Yavg(:,54) Yavg(:,106)];
% 
% RMSEsampleMLR = LOOCV_OLS(YsampleMLR,concavg);


% Multilinear regression on individual samples
% BsampleMLR = inv(YsampleMLR'*YsampleMLR)*YsampleMLR'*concavg;
% RMSEsampleMLR = norm(concavg - YsampleMLR*BsampleMLR);

% Multilinear regression on averaged values
% BavgMLR = inv(YavgMLR'*YavgMLR)*YavgMLR'*concavg;
% RMSEavgMLR = norm(concavg - YavgMLR*BavgMLR);
% 
nPC = 5;
% -------------------Principal Component Regression--------------------
% PCR on individual samples
% RMSEPCRsample = [];
% for i = 1:nPC
%       RMSEPCRsample = [RMSEPCRsample; LOOCV_PCR(Ysample,concavg,i)];
% end

% PCR on averaged values
% RMSEPCRavg = [];
% for i = 1:nPC
%     RMSEPCRavg = [RMSEPCRavg; LOOCV_PCR(Yavg,concavg,i)];
% end
% 
% -------------------Scaled Principal Component Regression--------------------
% Scaled PCR on individual samples
% RMSESPCRsample = [];
% stda = mean(stdavg);
% for i = 1:nPC
%     Ysamplescale = Ysample*inv(diag(stda)); % Does not matter if we scale use std of averaged values, since scaling factor does not make a difference to solution
%     RMSESPCRsample = [RMSESPCRsample; LOOCV_PCR(Ysamplescale, concavg, i)];
% end

% % Scaled PCR on averaged samples
% RMSESPCRavg = [];
% for i = 1:nPC
%     Yavgscale = Yavg*inv(diag(stda)); % Does not matter if we scale use std of averaged values, since scaling factor does not make a difference to solution
%     RMSESPCRavg = [RMSESPCRavg; LOOCV_PCR(Yavgscale, concavg, i)];
% end;

% -------------------Iterative Principal Component Regression--------------------
% % IPCR on individual samples
% RMSEIPCRsample = [];
% % Estimate std of errors by dividing the absorbance data into blocks of
% % manageable size 26 x 25 (wavelength range divided into blocks of 25
% nblocks = 7;
% %  For each choice of number of PCs
% for i = 1:nPC
%     sdall = [];
%     for j = 1:7
%         Yblock = Ysample(:,25*(j-1)+1:25*j);
%         if (j == 7)
%             Yblock = Ysample(:,25*(j-1)+1:end);
%         end
%         flag = 1;
%         nwave = size(Yblock,2);
%         sdold = ones(nwave,1);
%         Linv = eye(nwave);  % INitial guess of std is identity matrix
%         %  estimate std of errors for different wavelengths and iterate
%         %  till convergence or maximum iterations are exceeded
%         iter = 0;
%         while (flag)
%             iter = iter + 1;
%             Yblocks = Yblock*Linv;
%             [U S V] = svd(Yblocks,'econ');
%             Amat = V(:,i+1:end)'*Linv;
%             sdnew = stdest(Amat,Yblock');
%             if ( ( norm(sdold - sdnew) < 1.0E-04 ) || ( iter > 10 ) )
%                 flag = 0;
%             else
%                 sdold = sdnew;
%                 Linv = diag(1./sdnew);
%             end
%         end
%         sdall = [sdall; sdnew];  % Combine estimated error variances of all blocks into one vector
%     end
%     %
%     Yscale = Ysample*inv(diag(sdall)); % Scale using estimated error std for all wavelengths
%     %     [U S V] = svd(Yavgscale,0); % Economical svd
%     %     T = Yavgscale*V(:,1:nPC);  % Compute scores for chosen number of PCs
%     %     BPCR = inv(T'*T)*T'*concavg;
%     RMSEIPCRsample = [RMSEIPCRsample; LOOCV_PCR(Yscale, concavg, i)];
% end

% -------------Maximum Likelihood Principal Component Regression (Uncorrelated errors)--------------
% % MLPCR on individual samples
RMSEMLPCRsample = [];
for k = 1:nPC
    %
    % Perform leave one cross validation based on MLPCA for each choice of # PCs
    %
    [nsamples nspecies] = size(concavg);
    %
    % Estimate OLS regression matrix leaving one sample out in turn
    %
    RMSE = zeros(1,nspecies);
    %
    %  Build MLPCR model dropping each sample out in turn
    for i = 1:nsamples
        Xsub = [Ysample(1:i-1,:); Ysample(i+1:end,:)];
        Ysub = [concavg(1:i-1,:); concavg(i+1:end,:)];
        stdsub = [stdavg(1:i-1,:); stdavg(i+1:end,:)];
        [u s v,sobj,errflag] = MLPCA(Xsub,stdsub,k);
%         sobj
%         errflag
        V1 = v(:,1:k);
        Tsub = Xsub*V1;
        B = inv(Tsub'*Tsub)*Tsub'*Ysub;
        Sinv = inv(diag(stdavg(i,:)));  % Diagonal matrix containing variance of errors for sample i
        tpred = Ysample(i,:)*Sinv*V1*inv(V1'*Sinv*V1);
        prederr = concavg(i,:) - tpred*B;
        RMSE = RMSE + prederr.*prederr;
    end
    RMSE = sqrt(RMSE/nsamples);
    RMSEMLPCRsample = [RMSEMLPCRsample; RMSE];
end

