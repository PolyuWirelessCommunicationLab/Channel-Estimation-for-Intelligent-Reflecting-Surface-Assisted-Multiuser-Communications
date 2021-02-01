function coeff_es=coeffes2_MlN(M,H1_mtx_es,N,K,coeff,H1_mtx,H_drl_es,H_drl)
% define sets and generate pilots and phi.
      rho=floor(N/M);
      v=N-M*rho;
      PhaI=zeros(N-v,K);
      PhaII=zeros(v,K);  
      for k=2:K
          u_id_set=(k-2)*v+1:(k-1)*v;
          for i_v=1:v
              u_id_ele=u_id_set(i_v);
              PhaII(i_v,k)=u_id_ele-(ceil(u_id_ele/N)-1)*N;
          end        
      end
      for k=2:K
          N_set=1:N;
          N_set(PhaII(:,k))=[];
          PhaI(:,k)=N_set;
      end
% channel estimation
      coeff_es=zeros(N,K);
      % channel estimation for i=1:(K-1)*rho
      for i=1:(K-1)*rho
          ki=(i-(ceil(i/rho)-1)*rho-1)*M;
          u_id=ceil(i/rho)+1;
          Omega=PhaI(ki+1:ki+M,u_id);          
          y=H1_mtx(:,Omega)*coeff(Omega,u_id)+H_drl(:,u_id);
          y=y-H_drl_es(:,u_id);
          coeff_es(Omega,u_id)=inv(H1_mtx_es(:,Omega))*y;
      end
      % channel estimation for i=(K-1)*rho+1:ceil((K-1)*N/M)
      for i=(K-1)*rho+1:ceil((K-1)*N/M)
          J_set=(i-(K-1)*rho-1)*M+1:min((i-(K-1)*rho)*M,(K-1)*N-(K-1)*M*rho);
          u_set=[];
          irs_set=[];
          for j=1:size(J_set,2)
              j_ele=J_set(j);
              u_idx=ceil(j_ele/v)+1;
              u_set=[u_set,u_idx];
              irs_set=[irs_set PhaII(j_ele-(ceil(j_ele/v)-1)*v,u_idx)];
          end
          u_set_unq=unique(u_set);
          ref_lik=sum(H1_mtx(:,irs_set)*coeff(irs_set,u_set_unq),2);
          dir_lik=sum(H_drl(:,u_set_unq),2);
          dir_lik_es=sum(H_drl_es(:,u_set_unq),2);
          y=ref_lik+dir_lik;
          y=y-dir_lik_es;
          s_temp=0;
          for kk=1:length(u_set_unq)
              k=u_set_unq(kk);
              idx=find(u_set==k);
              S=irs_set;
              S(idx)=[];
              for nn=1:length(S)
                  n=S(nn);
                  s_temp=s_temp+H1_mtx_es(:,n)*coeff_es(n,k);
              end
              idx=[];
              S=[];
          end
          y=y-s_temp;
          coeff_es_temp=pinv(H1_mtx_es(:,irs_set))*y;
          for j=1:size(J_set,2)
              j_ele=J_set(j);
              u_idx=ceil(j_ele/v)+1;
              irs_idx=PhaII(j_ele-(ceil(j_ele/v)-1)*v,u_idx);
              coeff_es(irs_idx,u_idx)=coeff_es_temp(j);
          end
          J_set=[];
      end
end