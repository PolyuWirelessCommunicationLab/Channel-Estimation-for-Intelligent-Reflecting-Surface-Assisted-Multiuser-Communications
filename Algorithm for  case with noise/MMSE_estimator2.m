function MSE=MMSE_estimator2(M,K,N,T0,T,SNR,d_drl,coeff,G1,H_drl_mtx)
% Benchmark scheme
%% ============Phase I============================
F0_temp=dftmtx(T0);
A0=F0_temp(1:K,:); % pilots
P0_bit=1;

noise0=1/sqrt(2)*(randn(M,T0)+1j*randn(M,T0)); 
sigma2_0=P0_bit*10^(-SNR/10);
Y0=H_drl_mtx*A0+sqrt(sigma2_0)*noise0;
H_drl_mtx_es=[];
for i_h=1:K
   H_drl_mtx_es_temp=d_drl(i_h)/(d_drl(i_h)*T0+sigma2_0)*Y0*A0(i_h,:)';
   H_drl_mtx_es=[H_drl_mtx_es H_drl_mtx_es_temp];
end 
H_drl(1:M*K,1)=reshape(H_drl_mtx,M*K,1);
H_drl_es(1:M*K,1)=reshape(H_drl_mtx_es,M*K,1);
%% ===========Phase II============================
T1=T-T0;
P_bit=1;

H1_es=zeros(N*K*M,1);
H1_1=G1;
H1(1:M*N,1)=reshape(H1_1,M*N,1);
for i_ch=2:K
   H1_i=H1_1*diag(coeff((i_ch-2)*N+1:(i_ch-1)*N,1));
   H1((i_ch-1)*M*N+1:i_ch*M*N,1)=reshape(H1_i,M*N,1);         
end
sigma2=P_bit*10^(-SNR/10);
for u=1:T1
    noise=1/sqrt(2)*(randn(M,1)+1j*randn(M,1)); 
    idx_u=ceil(u/N);
    H1_i=H1((u-1)*M+1:u*M,1);
    H_drl_u=H_drl((idx_u-1)*M+1:idx_u*M,1);
    H_drl_es_u=H_drl_es((idx_u-1)*M+1:idx_u*M,1);
    y=H1_i+H_drl_u+sqrt(sigma2)*noise;
    y=y-H_drl_es_u;     
    H1_es((u-1)*M+1:u*M,1)=y;
end
H_es=[H_drl_es;H1_es];
H=[H_drl;H1];
cc=H'*H;
MSE=(H_es-H)'*(H_es-H)/cc;
