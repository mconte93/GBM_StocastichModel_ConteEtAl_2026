function [Con_term]=Conv_M(K,Drift_term)

%   Takes a 2x3 matrix with the nodes of the element and computes
%   the convection matrix corresponding to that element. 
%   It uses a III order Gaussian quadrature rule and can implement an 
%   anisotropic diffusion tensor. 

B=[K(:,2)-K(:,1) K(:,3)-K(:,1)];

AreaT=abs(det(B))./2;

gradFIhat=[-1 1 0; -1 0 1];

gradFI=(inv(B))'*gradFIhat;

bconP1=Drift_term(:,1);
bconP2=Drift_term(:,2);
bconP3=Drift_term(:,3);

bconP12=(bconP1+bconP2)./2;
bconP13=(bconP1+bconP3)./2;
bconP23=(bconP2+bconP3)./2;
bconP123=(bconP1+bconP2+bconP3)./3;

Dtensor1=AreaT./60.*(3.*(bconP1) + 4.*(bconP12+bconP13) + 9.*bconP123);
Dtensor2=AreaT./60.*(3.*(bconP2) + 4.*(bconP12+bconP23) + 9.*bconP123);
Dtensor3=AreaT./60.*(3.*(bconP3) + 4.*(bconP23+bconP13) + 9.*bconP123);

Con_term(:,1)=Dtensor1'*gradFI;
Con_term(:,2)=Dtensor2'*gradFI;
Con_term(:,3)=Dtensor3'*gradFI;

end



