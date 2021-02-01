function coeff_es=coeffes2_MgeN(H1_mtx_es,K,coeff,H1_mtx,H_drl_es,H_drl)
for k=2:K
    y=H1_mtx*coeff(:,k)+H_drl(:,k);
    y=y-H_drl_es(:,k);
    coeff_es(:,k)=pinv(H1_mtx_es)*y;
end