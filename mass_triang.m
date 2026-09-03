function M=mass_triang(K,U)

B=[K(:,2)-K(:,1) K(:,3)-K(:,1)];
AreaT=abs(det(B))./2;

fac1=U(1);
fac2=U(2);
fac3=U(3);

fac12=(fac1+fac2)/2;
fac13=(fac1+fac3)/2;
fac23=(fac3+fac2)/2;

fac123=(fac1+fac2+fac3)/3;

M=AreaT./60.*(3.*(fac1+fac2+fac3) + 8.*(fac12+fac13+fac23) + 27.*fac123);

end