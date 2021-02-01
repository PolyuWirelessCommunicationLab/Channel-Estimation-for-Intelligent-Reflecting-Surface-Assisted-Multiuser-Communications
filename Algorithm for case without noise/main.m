%% ======================================================
% This file is the channel estimation algorithm WITHOUT noise shown in Section IV 
% of the paper "Channel Estimation for Intelligent Reflecting Surface Assisted Multiuser Communications: Framework, 
% Algorithms, and Analysis" by Zhaorui Wang, Liang Liu, and Shuguang Cui.
% Channel Model: Correlated Racian Fading. It can be replaced by others as long as the channels are linearly independent
% Date: 01/12/2019
% Modified on:10/04/2020
% Modified on:19/06/2020. real correlation matrix--> complex correlation matrix; only NLOS channel
% All Rights Reserved
% For any code questions, please send email to zrwang2009@gmail.com.
%% ========================================================
clc
clear
M=32; % number of receive antennas
K=8; % number of users
N=32; % number of IRS elements
T0=K; % minimum number of training symbols in the first stage
T1=N; % minimum number of training symbols in the second stage
if M>=N
    T2=K-1; % minimum number of training symbols in the third stage
else
    T2=ceil(N*(K-1)/M); % minimum number of training symbols in the third stage
end
T=T0+T1+T2;
%% channel parameters
d0=0.01;
d_bs=d0*100^(-2.2); % pathloss between IRS and BS
d_ue=d0*(10*rand(K,1)+5).^(-2.1); % pathloss  between IRS and users
d_drl=d0*(10*rand(K,1)+100).^(-4.2); % pathloss of direct channels
c_irs=0.3+0.4*1j; % channel correlation
c_bs=0.5+0.3*1j;  % channel correlation
C_irs=C_gen(N,c_irs); % generate correlation matrix. Please refer to the file 'algorithm for the case with noise' for C_irss are different among different users. 
C_bs=C_gen(M,c_bs); % generate of correlation matrix Please refer to the file 'algorithm for the case with noise' for C_bss are different among different users. 
SNR=10000; % i.e., no noise. 
%% ===============Phase I===========================
F0_temp=dftmtx(T0);
A0=F0_temp(1:K,:); % pilot matrix

H_drl_mtx=[];
for i_h=1:K
   H_drl_temp=sqrtm(C_bs)*1/sqrt(2)*(randn(M,1)+1j*randn(M,1))*sqrt(d_drl(i_h));
   H_drl_mtx=[H_drl_mtx H_drl_temp];
end
Y0=H_drl_mtx*A0;
H_drl_mtx_es=Y0*A0'/T0;
%% =========Phase II============================== 
P1_bit=1;
F1=dftmtx(T1);
Phi1=F1(1:N,:); % phi in phase II;  none nomalized DFT matrix. MMSE stage:T=N
a1=diag(ones(T1,1))*sqrt(P1_bit); % pilot; a=diag(1,1,...,1)*sqrt(P1/N)

H_bs_hat=1/sqrt(2)*(randn(M,N)+1j*randn(M,N));
H_bs=sqrt(d_bs)*sqrtm(C_bs)*H_bs_hat*sqrtm(C_irs);
H_ue1_hat=1/sqrt(2)*(randn(N,1)+1j*randn(N,1));
H_ue1=sqrt(d_ue(1))*sqrtm(C_irs)*H_ue1_hat;
H1_mtx=H_bs*diag(H_ue1);
a1_drl=ones(1,T1);
Y1=H1_mtx*Phi1*a1+H_drl_mtx(:,1)*a1_drl; 
Y1=Y1-H_drl_mtx_es(:,1)*a1_drl;
H1_mtx_es=Y1*Phi1'/T1;
%% ========Phase III====================================
coeff=zeros(N,K); % note: the user 1' coeff is set to zero. 
for k=2:K
 % generating channels from user 2 to K.
 H_uek_hat=1/sqrt(2)*(randn(N,1)+1j*randn(N,1));
 H_uek=sqrt(d_ue(k))*sqrtm(C_irs)*H_uek_hat;
 coeff(:,k)=H_uek./H_ue1;
end
% channel estimation
if M>=N
  coeff_es=coeffes2_MgeN(H1_mtx_es,K,coeff,H1_mtx,H_drl_mtx_es,H_drl_mtx); % channel estimation algorithm for M>=N
else
  coeff_es=coeffes2_MlN(M,H1_mtx_es,N,K,coeff,H1_mtx,H_drl_mtx_es,H_drl_mtx); % channel estimation algorithm for M<N
end
%% recover all real channles and computing normalized MSE
coeff_es(:,1)=[];
coeff(:,1)=[];
coeff_es=reshape(coeff_es,(K-1)*N,1);
coeff=reshape(coeff,(K-1)*N,1);
H_drl_es=reshape(H_drl_mtx_es,M*K,1);
H_drl=reshape(H_drl_mtx,M*K,1);
H1_es=reshape(H1_mtx_es,M*N,1);
H1=reshape(H1_mtx,M*N,1);
H_es=[H_drl_es;H1_es];
H_real=[H_drl;H1];
for iii=1:K-1
  H_es_temp=H1_mtx_es*diag(coeff_es((iii-1)*N+1:iii*N));
  H_real_temp=H1_mtx*diag(coeff((iii-1)*N+1:iii*N));
  H_es_temp2=reshape(H_es_temp,N*M,1);
  H_real_temp2=reshape(H_real_temp,N*M,1);
  H_es=[H_es;H_es_temp2];
  H_real=[H_real;H_real_temp2];       
end  
p_s=H_real'*H_real; % total channel power approximation
MSE=(H_es-H_real)'*(H_es-H_real)/p_s; 
fprintf('Normalized Mean Square Error: '); disp(MSE);
if MSE<10^(-10)
   fprintf('Perfect Channel Estimation');
end