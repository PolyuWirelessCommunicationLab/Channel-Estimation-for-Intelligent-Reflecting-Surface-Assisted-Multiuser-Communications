function [MSE]=Proposed_MMSE2_6(T,T_ini,T0_min,T1_min,T2_min,M,K,N,SNR,d_drl,C_bs_ur,G1_cov_inv,lambda_cov,lambda_pwr,coeff,G1,H_drl_mtx)
%% ======================================
% This file is the channel estimation algorithm WITH noise shown in Section V 
% of the paper "Channel Estimation for Intelligent Reflecting Surface Assisted Multiuser Communications: Framework, 
% Algorithms, and Analysis" by Zhaorui Wang, Liang Liu, and Shuguang Cui.
% Zhaorui WANG, 10/04/2020
% G1_mean: mean of the user 1's refecting channel
% G1_cov : covariance matrix of the user 1's refecting channel
% G1 : user 1's overall reflecting channel
% lambda_mean: mean of reflecting channel ratios between user k and user 1
% lambda_cov : covariance matrix of reflecting channel ratios between user k and user 1
% lambda_pwr : power of reflecting channel ratios between user k and user 1
% coeff  : reflecting channel ratios between user k and user 1
% H_drl_mtx: direct channels between BS and users
% C_bs_ur      : channel correlation matrix at BS
%% =========================================
% time allocation scheme
T_res=T-T_ini;
% T0=T0_min+T_res; T1=T1_min;       T2=T2_min; %scheme 1;
T0=T0_min;       T1=T1_min+T_res; T2=T2_min; %scheme 2;
% T0=T0_min;       T1=T1_min;       T2=T2_min+T_res; %scheme 3;
% T0=T0_min+floor(T_res/3); T1=T1_min+floor(T_res/3); T2=T-T0-T1; %scheme 4;
% T0=T0_min+floor(T_res*T0_min/T_ini); T1=T1_min+floor(T_res*T1_min/T_ini); T2=T-T0-T1; %scheme 5;
%% ===============Phase I===========================
F0_temp=dftmtx(T0);
A0=F0_temp(1:K,:); % pilots
P0_bit=1; % transmit power 

noise0=1/sqrt(2)*(randn(M,T0)+1j*randn(M,T0)); 
sigma2_0=P0_bit*10^(-SNR/10);
Y0=H_drl_mtx*A0+sqrt(sigma2_0)*noise0;
H_drl_mtx_es=[];
for i_h=1:K
   H_drl_mtx_es_temp=d_drl(i_h)/(d_drl(i_h)*T0+sigma2_0)*Y0*A0(i_h,:)';
   H_drl_mtx_es=[H_drl_mtx_es H_drl_mtx_es_temp];
end 
%% =========Phase II============================== 
F1=dftmtx(T1);
Phi1=F1(1:N,:); % phi; none nomalized DFT matrix. MMSE stage:T=N
P1_bit=1; % transmit power 
a1=diag(ones(T1,1))*sqrt(P1_bit); % pilot; a=diag(1,1,...,1)*sqrt(P1/N);

noise1=1/sqrt(2)*(randn(M,T1)+1j*randn(M,T1));     
sigma2_1=P1_bit*10^(-SNR/10);
H1_mtx=G1;
a1_drl=ones(1,T1);
Y1=H1_mtx*Phi1*a1+H_drl_mtx(:,1)*a1_drl+sqrt(sigma2_1)*noise1; 
Y1=Y1-H_drl_mtx_es(:,1)*a1_drl;
% channel estimation
err_0_1=M*d_drl(1)*sigma2_1/(d_drl(1)*T0+sigma2_1);
Psi=err_0_1*ones(T1)+M*sigma2_1*eye(T1);
Psi_inv=inv(Psi);
temp3=Psi_inv*Phi1'*inv(Phi1*Psi_inv*Phi1'+G1_cov_inv);
H1_mtx_es=(Y1)*temp3;
%% ========Phase III====================================  
% Phi2: active elements index range on IRS
% A2: allocated time slots for each user.
[Phi2,A2]=pdc_mtx2(M,K,N,T2); 
P2_bit=1; % transmit power

s=1;
sigma2_2=P2_bit*10^(-SNR/10);
coeff_es=[];
mse_thy=0;
for k=1:K-1
    % channel estimation for user k
    % produce received signals at each slot for user k;
    C_bs=C_bs_ur(k*M+1:(k+1)*M,:);
    lambda_cov_k=lambda_cov((k-1)*N+1:(k)*N,:); % cov of user k's coefficients
    tk=A2(k,1); % allocated time slots for user k
    H_drl_k=H_drl_mtx(:,k+1); % k-th direct link
    H_drl_es_k=H_drl_mtx_es(:,k+1); % k-th estimated-direct link
    cov=sigma2_2*eye(M*ceil(N/M)); % noise covariance matrix within M*ceil(N/M) 
    y2=zeros(tk*M,1);   % recieved signals within tk*M   
    coeff_k=coeff((k-1)*N+1:k*N,1); %k-th coefficients
    Phi2_k=Phi2(s:s+tk-1,:);  % allocated IRS-elements for user k
    phi2=diag(ones(1,N)); % values of IRS-elements for user k
    a2_k=ones(1,tk);  % values of pilots for user k
    for t=1:tk
      Phi2_kt=Phi2_k(t,1):Phi2_k(t,2); % active irs elements at time t for user k
      phi2_tk=phi2(Phi2_kt,Phi2_kt); % values of active irs elements at time t for user k
      a2_tk=a2_k(t); % values of pilots for user k at time t
      coeff_kt=coeff_k(Phi2_kt); % values of coefficients at time t for user k
      H1_mtx_t=H1_mtx(:,Phi2_kt); % active usr 1's reflecting channels
      noise2=1/sqrt(2)*(randn(M,1)+1j*randn(M,1));
      y2((t-1)*M+1:t*M,1)=H_drl_k*a2_tk+a2_tk*H1_mtx_t*phi2_tk*coeff_kt+sqrt(P2_bit*10^(-SNR/10))*noise2;
      Phi2_kt=[];
      coeff_kt=[];
      H1_mtx_t=[];
    end
    s=s+tk;
    % received signal combination for the estimaion of a same channel.
    num_g=floor(tk/ceil(N/M));
    num_r=tk-num_g*ceil(N/M);
    y2_temp=reshape(y2(1:num_g*ceil(N/M)*M),ceil(N/M)*M,num_g);
    Y2=sum(y2_temp,2);
    for i_y=1:num_r
         Y2((i_y-1)*M+1:i_y*M,1)=Y2((i_y-1)*M+1:i_y*M,1)+y2(num_g*ceil(N/M)*M+(i_y-1)*M+1:num_g*ceil(N/M)*M+i_y*M,1);
         Y2((i_y-1)*M+1:i_y*M,1)=Y2((i_y-1)*M+1:i_y*M,1)/(num_g+1);
         cov((i_y-1)*M+1:i_y*M,(i_y-1)*M+1:i_y*M)=cov((i_y-1)*M+1:i_y*M,(i_y-1)*M+1:i_y*M)/(num_g+1);
    end 
    for i_y=num_r+1:ceil(N/M)
          Y2((i_y-1)*M+1:i_y*M,1)=Y2((i_y-1)*M+1:i_y*M,1)/(num_g);
          cov((i_y-1)*M+1:i_y*M,(i_y-1)*M+1:i_y*M)=cov((i_y-1)*M+1:i_y*M,(i_y-1)*M+1:i_y*M)/(num_g);
    end
    % estimate channels
    Phi2_kk=Phi2_k(1:ceil(N/M),:); % IRS index
    a2_k_up=a2_k;  % pilots
    for i_t2=1:ceil(N/M)
         a2_tk=a2_k_up(i_t2); % pilots
         Phi2_ik=Phi2_kk(i_t2,1):Phi2_kk(i_t2,2); % IRS index
         phi2_tk=phi2(Phi2_ik,Phi2_ik); %phase
         lambda_cov_kt=lambda_cov_k(Phi2_ik,Phi2_ik); % cov of user k's coefficients at time i_t2
         cov_i=cov((i_t2-1)*M+1:i_t2*M,(i_t2-1)*M+1:i_t2*M); % covariance of noise
         cov_err_k=(sigma2_2/(d_drl(k+1)*T0+sigma2_2))^2*C_bs*d_drl(k+1)...
             +(d_drl(k+1)/(d_drl(k+1)*T0+sigma2_2))^2*T0*sigma2_2*eye(M);
         cov_i_allnoise=cov_i+cov_err_k;
         cov_i_allnoise_inv=inv(cov_i_allnoise);
         Y2_i=Y2((i_t2-1)*M+1:i_t2*M)-H_drl_es_k*a2_tk; 
         V=H1_mtx_es(:,Phi2_ik);
         temp1=inv(inv(lambda_cov_kt)+phi2_tk'*V'*cov_i_allnoise_inv*V*phi2_tk)*a2_tk'*phi2_tk'*V'*cov_i_allnoise_inv;  
         coeff_temp=temp1*(Y2_i);
         coeff_es=[coeff_es;coeff_temp];
         temp2=inv(inv(lambda_cov_kt)+phi2_tk'*V'*cov_i_allnoise_inv*V*phi2_tk);
         mse_thy=mse_thy+abs(trace(temp2));
         lambda_cov_kt=[];
         Y2_i=[];
         V=[];
         Phi2_ik=[];
         cov_i=[];
         cov_i_allnoise=[];
         coeff_temp=[];
         phi2_tk=[];
    end
    tk=[];
    H_drl_k=[];
    H_drl_k=[];
    y2=[];
    coeff_k=[];
    Phi2_k=[];
    cov=[];
end
    % recover all real channles
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
cc=H_real'*H_real;
MSE=(H_es-H_real)'*(H_es-H_real)/cc; %!!