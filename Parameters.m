T=86400; %d
h_vox=0.875;   %mm 

%%%% Interins equation %%%%
k1plus=0.034*T;              %1/d     
k2plus=0.034*T;              %1/d     
kminus=0.01*T;              %1/d

%%%% Tumor equation %%%%%
lambda0=0.0001*T;                   %1/d  
lambda1=0.0001*T;                  %1/d
KM=10^5;                          %cell/mm^3
speed_M=10^(-3)*0.0084*T;           %mm/s                 
mu_M=8.44*10^(-7)*T;               %1/d
dr=2;
alfaM=0.0906;   %Tra 0.025 e 0.033 
betaM=0.006;

%%%% Tissue equation %%%
d_Q=0.058*10^(-8)*T;
alfaQ=0.0025/10;  %0.025
betaQ=0.00005; %0.0025
r_Q=0; %mu_M/5;

%%%% Additional parameters %%%% 
acc=(speed_M^2)/lambda0;   %nodim                               
var_eps=.001;
