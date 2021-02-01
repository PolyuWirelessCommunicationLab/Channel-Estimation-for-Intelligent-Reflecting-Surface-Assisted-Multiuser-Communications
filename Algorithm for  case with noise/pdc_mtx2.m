function [Phi2,A2]=pdc_mtx2(M,K,N,T2)
% Phi2: active elements index range on IRS
% allocated time slots for each user.
% Wang zhaorui, 01/12/2019
g1=floor(T2/((K-1)*ceil(N/M)));
r1=T2-g1*(K-1)*ceil(N/M);
g2=floor(r1/ceil(N/M));
r2=r1-g2*ceil(N/M);
A2=ones(K-1,1)*ceil(N/M)*g1;
A2(1:g2,1)=A2(1:g2,1)+ceil(N/M);
A2(1+g2,1)=A2(1+g2,1)+r2;

c_temp=ceil(N/M);
Phi_temp=zeros(c_temp,2);
Phi2=[];
for i=1:c_temp
    Phi_temp(i,:)=[(i-1)*M+1 i*M];
    if i==c_temp
      Phi_temp(i,:)=[(i-1)*M+1 N];  
    end
end
for ii=1:K-1
    a_i=A2(ii,1);
    g1a=floor(a_i/c_temp);
    r1a=a_i-g1a*c_temp;
    for iii=1:g1a
        Phi2=[Phi2;Phi_temp];
    end
    Phi2=[Phi2;Phi_temp(1:r1a,:)];
end
end