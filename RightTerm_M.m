function F_N=RightTerm_M(K,U,P_M)


% P1=K(:,1);
% P2=K(:,2);
% P3=K(:,3);
% 
B=[K(:,2)-K(:,1) K(:,3)-K(:,1)]; 
AreaT=abs(det(B))./2;  


TP1=U(1).*P_M(1);

TP2=U(1).*P_M(2);

TP3=U(3).*P_M(3);

%%%%%%%
TP12=0.5*(TP1+TP2);
TP13=0.5*(TP1+TP3);
TP23=0.5*(TP3+TP2);
TP123=1/3 *(TP1+TP2+TP3);

F_N=AreaT./60*[3*TP1+4*(TP12+TP13)+9*TP123   3*TP2+4*(TP12+TP23)+9*TP123   3*TP3+4*(TP13+TP23)+9*TP123]';
%%%%%%%
%F_N=AreaT./3.*[TP1 TP2 TP3]';

end





