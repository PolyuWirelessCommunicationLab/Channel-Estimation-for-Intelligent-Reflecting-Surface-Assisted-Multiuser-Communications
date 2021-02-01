function [G1_cov,G1_pwr,G1_all,lambda_cov,lambda_pwr,coeff_out,H_drl_mtx_all,C_bs_ur]=channel_parameters(M,N,K,d_drl,d_bs,d_ue,N_para)
%% ==========================
% generating channels
% Zhaorui WANG, 10/04/2020
% modified on: 06/19/2020
% G1_cov : covariance matrix of the user 1's refecting channel
% G1_pwr : user 1's refelcting channel power
% G1_all : user 1's overall reflecting channel
% lambda_cov : covariance matrix of reflecting channel ratios between user k and user 1
% lambda_pwr : power of reflecting channel ratios between user k and user 1
% coeff_out  : reflecting channel ratios between user k and user 1
% H_drl_mtx_all: direct channels between BS and users
% C_bs_ur      : channel correlation matrix at BS
%% ===========================
rand_c=[0.1+1j*0.2 0.2+1j*0.1 0.1+1j*0.1 0.2+1j*0.4 0.2+1j*0.3 0.2+1j*0.3 0.3+1j*0.3 0.2+1j*0.1 0.1+1j*0.3 0.1+1j*0.2 0.1]; % correlation value
c_irs=0.3+1j*0.2; % correlation 
c_bs=0.1+1j*0.2; % correlation
C_irs=C_gen(N,c_irs); % generating of square root correlation matrix
C_bs=C_gen(M,c_bs); % generating of square root correlation matrix
c_irs_ur=rand_c(1:K); % correlation
c_bs_ur=rand_c(1:K); % correlation
C_irs_ur=[];
C_bs_ur=[];
for k=1:K
    temp1=C_gen(M,c_bs_ur(k));
    temp2=C_gen(N,c_irs_ur(k));
    C_bs_ur=[temp1;C_bs_ur];
    C_irs_ur=[temp2;C_irs_ur];
end
%
G1_cov=zeros(N,N);
G1_pwr=zeros(1,1);
G1_all=zeros(M*N_para,N);
lambda_cov=zeros((K-1)*N,N);
lambda_pwr=0;
H_drl_mtx_all=[];
for i=1:N_para
    % direct channels   
    for i_h=1:K
       C_bs_k=C_bs_ur((i_h-1)*M+1:i_h*M,:);
       H_drl_mtx(:,i_h)=C_bs_k*1/sqrt(2)*(randn(M,1)+1j*randn(M,1))*sqrt(d_drl(i_h));
    end
    H_drl_mtx_all=[H_drl_mtx;H_drl_mtx_all];
    % reflecting chanel of user 1
    H_bs_hat=1/sqrt(2)*(randn(M,N)+1j*randn(M,N));
    H_bs=sqrt(d_bs)*C_bs*H_bs_hat*C_irs;
    C_irs_1=C_irs_ur(1:N,:);
    H_ue1_hat=1/sqrt(2)*(randn(N,1)+1j*randn(N,1));
    H_ue1=sqrt(d_ue(1))*C_irs_1*H_ue1_hat;
    H1_mtx=H_bs*diag(H_ue1);
    G1_all((i-1)*M+1:i*M,:)=H1_mtx;
    G1_cov=G1_cov+H1_mtx'*H1_mtx;
    H1_clm=reshape(H1_mtx,M*N,1);
    G1_pwr=G1_pwr+H1_clm'*H1_clm;
    H_ue1_hat=1/sqrt(2)*(randn(N,1)+1j*randn(N,1));
    H_ue1=sqrt(d_ue(1))*C_irs_1*H_ue1_hat;
    % reflecting channel of users 2-K; coefficients
      for k=2:K
         % generating matrix from user 2 to K.
         C_irs_k=C_irs_ur((k-1)*N+1:k*N,:);
         H_uek_hat=1/sqrt(2)*(randn(N,1)+1j*randn(N,1));
         H_uek=sqrt(d_ue(k))*C_irs_k*H_uek_hat;
         coeff((k-2)*N+1:(k-1)*N,i)=H_uek./H_ue1;
         lambda_cov((k-2)*N+1:(k-1)*N,:)=lambda_cov((k-2)*N+1:(k-1)*N,:)+coeff((k-2)*N+1:(k-1)*N,i)*coeff((k-2)*N+1:(k-1)*N,i)';
      end
      lambda_pwr=lambda_pwr+coeff(:,i)'*coeff(:,i);
end
G1_cov=G1_cov/N_para;
G1_pwr=G1_pwr/N_para;
coeff_out=coeff;
for k=2:K
    lambda_cov((k-2)*N+1:(k-1)*N,:)=lambda_cov((k-2)*N+1:(k-1)*N,:)/N_para;
end
lambda_pwr=lambda_pwr/N_para;
