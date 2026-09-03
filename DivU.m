function [DivU] = DivU(K,U)

%Function for the gradient of the solution U at the element K

P1=K(:,1);
P2=K(:,2);
P3=K(:,3);

V=P2-P1;
W=P3-P1;
T=P3-P2;

M1=[V';W'];

b1=[(U(2)-U(1)); U(3)-U(1)];
GradU1=M1\b1; 

M2=[-V';T'];

b2=[(U(1)-U(2)); U(3)-U(2)];
GradU2=M2\b2;


M3=[-W';-T'];

b3=[(U(1)-U(3)); U(2)-U(3)];
GradU3=M3\b3; 

DivU=[GradU1 GradU2 GradU3];

end