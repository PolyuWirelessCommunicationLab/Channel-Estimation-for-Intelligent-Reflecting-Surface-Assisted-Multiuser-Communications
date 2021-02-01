function C=C_gen(D,c_bs)
for i=1:D
    for j=1:D
        if j>=i
           C_temp(i,j)=c_bs^(j-i);
        else
           C_temp(i,j)=C_temp(j,i)';
        end
    end
end
C=sqrtm(C_temp);