%% ======================================================
% This file is the channel estimation algorithm WITH noise shown in Section V 
% of the paper "Channel Estimation for Intelligent Reflecting Surface Assisted Multiuser Communications: Framework, 
% Algorithms, and Analysis" by Zhaorui Wang, Liang Liu, and Shuguang Cui.
% Channel Model: Correlated Racian Fading. It can be replaced by others as long as the channels are linearly independent
% Date: 01/12/2019
% Modified on:10/04/2020
% Modified on:19/06/2020. real correlation matrix--> complex correlation matrix; only NLOS channel
% All Rights Reserved
% For any code questions, please send email to zrwang2009@gmail.com.
%% ==================================================================
clc
clear
iteration=20000; % 20000 times for a stable and precise result
M=32; % number of receive antennas
K=8; % number of users
N=32; % number of IRS elements
T0_min=K; % minimum number of training symbols in the first stage
T1_min=N; % minimum number of training symbols in the second stage
if M>=N
    T2_min=K-1; %  training symbols in the third stage
else
    T2_min=(K-1)*ceil(N/M); %  training symbols in the third stage
end
T_ini=T0_min+T1_min+T2_min;
T=[50:10:100];
%% Generate channels for Monte Caro Simulation
% channel parameters
d0=0.01;
d_bs=d0*100^(-2.2); % pathloss between IRS and BS
d_ue_all=[9.92929414915431e-05;0.000933140557886915;8.93686991222412e-05;0.000895057527534615;7.92296698429811e-05;9.80306203829406e-05;8.92266090070821e-05;0.000905361709970152;8.92929414915431e-05;0.000933140557886915;9.83686991222412e-05;0.000245057527534615;3.62296698429811e-05;9.10306203829406e-05;8.22266090070821e-05;0.000105361709970152];
d_ue=d_ue_all(1:K,1); % pathloss  between IRS and users
d_drl_all=[9.93910190825605e-11;9.97952844902355e-11;8.94221515470294e-11;9.93355022335517e-11;8.97219456983751e-11;9.94251035209681e-11;9.82479583303659e-11;9.30203404175918e-11;9.03910190825605e-11;8.97952844902355e-11;9.84221515470294e-11;2.83355022335517e-11;3.97219456983751e-11;2.94251035209681e-11;2.82479583303659e-11;3.30203404175918e-11];
d_drl=d_drl_all(1:K,1); % pathloss of direct channels
SNR=140;
% generate channels and the channel statistics for Monte Caro Simulation
N_para=iteration; % number of iterations
[G1_cov,G1_pwr,G1_all,lambda_cov,lambda_pwr,coeff_all,H_drl_mtx_all,C_bs_ur]=channel_parameters(M,N,K,d_drl,d_bs,d_ue,N_para); 
G1_cov_inv=inv(G1_cov);
for ii=1:iteration
    ii 
    tic
    % channels for ii-th iteration
    coeff=coeff_all(:,ii); % reflecting channel ratio
    G1=G1_all((ii-1)*M+1:ii*M,:); % user1 's reflecting channel 
    H_drl_mtx=H_drl_mtx_all((ii-1)*M+1:ii*M,:); % all users direct channels
    for i=1:length(T)
        % Benchmark scheme, i.e., estimated channels with traditional method
        mse_tr(ii,i)=MMSE_estimator2(M,K,N,T0_min,T(i),SNR,d_drl,coeff,G1,H_drl_mtx); 
        % Our proposed scheme. Allocating All Extra training symbols to Phase II 
        mse_pr2(ii,i)=Proposed_MMSE2_6(T(i),T_ini,T0_min,T1_min,T2_min,M,K,N,SNR,d_drl,...
            C_bs_ur,G1_cov_inv,lambda_cov,lambda_pwr,coeff,G1,H_drl_mtx); 
    end
    toc
    MSE_tr=mean(mse_tr,1);
    MSE_pr2=mean(mse_pr2,1);
end
%%
semilogy(T,MSE_tr,'r-','LineWidth',1);
xlabel('Overall Pilot Sequence Length');
ylabel('Normalized MSE');
hold on
grid on
set(gca,'XTick',T)
semilogy(T,MSE_pr2,'g-','LineWidth',1);
legend('Benchmark Scheme','Allocating All Extra Symbols to Phase II')