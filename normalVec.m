function nor_vec=normalVec(Bordo_N)

  N=zeros(2,size(Bordo_N,1));
  x_vec1=Bordo_N(end,1)-Bordo_N(1,1);
  y_vec1=Bordo_N(end,2)-Bordo_N(1,2);
  x_vec2=Bordo_N(2,1)-Bordo_N(1,1);
  y_vec2=Bordo_N(2,2)-Bordo_N(1,2);
  
  v1=[x_vec1 y_vec1]';
  v1=v1./norm(v1);
  v2=[x_vec2 y_vec2]';
  v2=v2./norm(v2);
 
  Rot_CW=[0 1;-1 0];
  Rot_ACW=[0 -1;1 0];
  n1=Rot_ACW*v1;
  n2=Rot_CW*v2;
  
  N(:,1)=n1+n2;
  N(:,1)=N(:,1)./norm(N(:,1));

  for i=2:size(Bordo_N,1)-1
     x_vec1=Bordo_N(i-1,1)-Bordo_N(i,1);
     y_vec1=Bordo_N(i-1,2)-Bordo_N(i,2);
     x_vec2=Bordo_N(i+1,1)-Bordo_N(i,1);
     y_vec2=Bordo_N(i+1,2)-Bordo_N(i,2); 
      
     v1=[x_vec1 y_vec1]';
     v1=v1./norm(v1);
     v2=[x_vec2 y_vec2]';
     v2=v2./norm(v2);
     
     n1=Rot_ACW*v1;
     n2=Rot_CW*v2;
     N(:,i)=n1+n2;
     N(:,i)=N(:,i)./norm(N(:,i));
  end
   
  x_vec1=Bordo_N(end-1,1)-Bordo_N(end,1);
  y_vec1=Bordo_N(end-1,2)-Bordo_N(end,2);
  x_vec2=Bordo_N(1,1)-Bordo_N(end,1);
  y_vec2=Bordo_N(1,2)-Bordo_N(end,2);
  
  v1=[x_vec1 y_vec1]';
  v1=v1./norm(v1);
  v2=[x_vec2 y_vec2]';
  v2=v2./norm(v2);
 
  n1=Rot_ACW*v1;
  n2=Rot_CW*v2;
  
  N(:,end)=n1+n2;
  N(:,end)=N(:,end)./norm(N(:,end));
  
  nor_vec=N;
end