function [CON_M,b_M] = matrices_M(C,L,C_new_ODF,X,Pix,F_Q,F_M,Q,P_M,h_vox)

%MATRICES_M: Assembles the convection matrix and right-hand side vector for the tumor model.

N=size(C,2);   %number od vertices
s=size(L);     %number of triangles

load ('Q1_derQ1_match.mat','DivDT_mesh')
DerQ=gradSol2D(C,L,Q);

U=X(1:N);
DerM=gradSol2D(C,L,U);

driftTerm_N=DriftTerm_M(C,C_new_ODF,DivDT_mesh,DerQ,DerM,F_Q,F_M,Pix,h_vox);

CON_M=sparse(N,N);
b_M=zeros(N,1); 

for k=1:s(1)
  K=C(1:2,L(k,:));
  Con=Conv_M(K,driftTerm_N(3:4,L(k,:)));
  RT_M=RightTerm_M(K,U(L(k,:)),P_M(L(k,:))); 
     
  for i=1:3
       for j=1:3 
          CON_M(L(k,i),L(k,j))=CON_M(L(k,i),L(k,j)) + Con(i,j);
       end
        b_M(L(k,i))=b_M(L(k,i)) + RT_M(i);
  end

end