%% Illustrative figure for initial setting - Entire Brain

clear all
close all

load ('Mesh_2E_mod.mat')
C=C_mod; 
L=L_mod;
h_vox=0.0875;

load ('Evo_noTer_th16_05margin')
U_end=U(:,end);

GTV_T2=find(U_end>=1/100*16);
x_GTV_T2=C(1,GTV_T2)';
y_GTV_T2=C(2,GTV_T2)';
bordo_GTV_T2=boundary(x_GTV_T2,y_GTV_T2,0.1);
Bordo_GTV_T2(:,1)=x_GTV_T2(bordo_GTV_T2);
Bordo_GTV_T2(:,2)=y_GTV_T2(bordo_GTV_T2);

Bordo_GTV_T2_mod=Bordo_GTV_T2(1:end-1,:);   
nor_vec=normalVec(Bordo_GTV_T2_mod);
nor_vec=nor_vec';
Bordo_PTV_T2=zeros(size(Bordo_GTV_T2_mod,1),2);
   
for i=1:size(Bordo_GTV_T2_mod,1)
   xB=Bordo_GTV_T2_mod(i,1);
   yB=Bordo_GTV_T2_mod(i,2);
  
   c_temp=11.43.*h_vox.*nor_vec(i,:);  %5.71=ratio between margin (0.5cm) and h_vox (0.0875cm)
   
   Bordo_PTV_T2(i,1)=xB+c_temp(1);
   Bordo_PTV_T2(i,2)=yB+c_temp(2);
end

Bordo_PTV_T2=[Bordo_PTV_T2;[Bordo_PTV_T2(1,1) Bordo_PTV_T2(1,2)]];

subplot(1,4,1)
trisurf(L,C(1,:),C(2,:),U(:,1),'FaceColor','interp','Edgecolor','none','Facelighting','phong')
view(2)
axis equal
axis([-70 70 -100 80])
xticks([-60 -20 20 60 ])
yticks([-100 -60 -20 20 60])
colormap jet
hold on
title(['$M(t=0)$'],'Interpreter','latex')
clim([0 1])

subplot(1,4,2)
trisurf(L,C(1,:),C(2,:),U(:,end),'FaceColor','interp','Edgecolor','none','Facelighting','phong')
view(2)
axis equal
axis([-70 70 -100 80])
xticks([-60 -20 20 60 ])
yticks([-100 -60 -20 20 60])
colormap jet
hold on
title(['$M(t=84d)$'],'Interpreter','latex')
clim([0 1])

subplot(1,4,3)
trisurf(L,C(1,:),C(2,:),U(:,end),'FaceColor','interp','Edgecolor','none','Facelighting','phong')
view(2)
axis equal
axis([-70 70 -100 80])
xticks([-60 -20 20 60 ])
yticks([-100 -60 -20 20 60])
colormap jet
hold on
clim([0 1])
hold on
title(['GTV'],'Interpreter','latex')
plot3(Bordo_GTV_T2(:,1),Bordo_GTV_T2(:,2),1*ones(size(Bordo_GTV_T2(:,1),1)),'r',LineWidth=0.7,LineStyle='--')

subplot(1,4,4)
trisurf(L,C(1,:),C(2,:),U(:,end),'FaceColor','interp','Edgecolor','none','Facelighting','phong')
view(2)
axis equal
axis([-70 70 -100 80])
xticks([-60 -20 20 60 ])
yticks([-100 -60 -20 20 60])
colormap jet
hold on
clim([0 1])
hold on
title(['PTV'],'Interpreter','latex')
plot3(Bordo_PTV_T2(:,1),Bordo_PTV_T2(:,2),1*ones(size(Bordo_PTV_T2(:,1),1)),'r',LineWidth=0.7,LineStyle='--')

 
%% Additional countour plot for T2-no th comparison

clear all
close all

load ('Mesh_2E_mod.mat')
C=C_mod; 
L=L_mod;

load ('Evo_noTer_th16_05margin')
U_end=U(:,end);

GTV_T2=find(U_end>=1/100*16);
x_GTV_T2=C(1,GTV_T2)';
y_GTV_T2=C(2,GTV_T2)';
bordo_GTV_T2=boundary(x_GTV_T2,y_GTV_T2,0.1);
Bordo_GTV_T2(:,1)=x_GTV_T2(bordo_GTV_T2);
Bordo_GTV_T2(:,2)=y_GTV_T2(bordo_GTV_T2);
  
GTV_noTh=find(U_end>0.001);
x_GTV_noTh=C(1,GTV_noTh)';
y_GTV_noTh=C(2,GTV_noTh)';
bordo_GTV_noTh=boundary(x_GTV_noTh,y_GTV_noTh,0.1);
Bordo_GTV_noTh(:,1)=x_GTV_noTh(bordo_GTV_noTh);
Bordo_GTV_noTh(:,2)=y_GTV_noTh(bordo_GTV_noTh);
   
subplot(1,2,1)
plot3(Bordo_GTV_T2(:,1),Bordo_GTV_T2(:,2),1*ones(size(Bordo_GTV_T2(:,1),1)),'r',LineWidth=0.7)
hold on
plot3(Bordo_GTV_noTh(:,1),Bordo_GTV_noTh(:,2),1*ones(size(Bordo_GTV_noTh(:,2),1)),'r--',LineWidth=0.7)
trisurf(L,C(1,:),C(2,:),U(:,end),'FaceColor','interp','Edgecolor','none','Facelighting','phong')
view(2)
axis equal
axis([-70 70 -100 80])
xticks([-60 -20 20 60 ])
yticks([-100 -60 -20 20 60])
colormap turbo
hold on
clim([0 1])

subplot(1,2,2)
plot3(Bordo_GTV_T2(:,1),Bordo_GTV_T2(:,2),1*ones(size(Bordo_GTV_T2(:,1),1)),'r',LineWidth=1.5)
hold on
plot3(Bordo_GTV_noTh(:,1),Bordo_GTV_noTh(:,2),1*ones(size(Bordo_GTV_noTh(:,2),1)),'r--',LineWidth=1.5)
trisurf(L,C(1,:),C(2,:),U(:,end),'FaceColor','interp','Edgecolor','none','Facelighting','phong')
view(2)
axis equal
axis([-35 5 -15 25])
xticks([-35 -25 -15 -5 5])
yticks([-15 -5 5 15 25])
colormap turbo
hold on
clim([0 1])
