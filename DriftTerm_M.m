function driftTerm_M = DriftTerm_M(C,C_new,DivDT,DerQ,DerM,F_Q,F_M,Pix,h_vox)

N=size(C,2);

coefDrift=zeros(2,N);
driftTerm_M=zeros(4,N);

for i=1:N
    x=C(1,i);
    y=C(2,i);
    [col]=find_voxel(x,y,Pix,h_vox);
    
    if col==0
       DT=zeros(2);
    else
       DT=[C_new(3,col) C_new(4,col); C_new(4,col) C_new(5,col)];
    end
    
    coefDrift(:,i)=F_Q(i).*DerQ(:,i)-F_M(i).*DerM(:,i);
    driftTerm_M(1:2,i)=C(1:2,i);
    driftTerm_M(3:4,i)=DT*coefDrift(:,i)-DivDT(3:4,i);

end

end