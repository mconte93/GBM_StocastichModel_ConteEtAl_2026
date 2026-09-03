clear all
close all
Tmax1=7*14;

%% Calculate total TCP,NTCT,UTCP 
Num_sim=195;
TCP_2perc=zeros(Tmax1,1);
NTCP=zeros(Tmax1,1);
UTCP_2perc=zeros(Tmax1,1);

for i=1:Tmax1
    numM_1_5perc=0;
    numM_2perc=0;
    numQ=0;
    for j=1:Num_sim
         load (['Results_DoseEx_15day/Evo15day_Num_',num2str(j),'.mat'])
         if tau_M_2perc>0 && tau_M_2perc<=i
             numM_2perc=numM_2perc+1;
         end
         if tau_Q>0 && tau_Q<=i
             numQ=numQ+1;
         end
    end

    TCP_2perc(i)=numM_2perc/(Num_sim);
    NTCP(i)=numQ/(Num_sim);
    UTCP_2perc(i)=TCP_2perc(i)*(1-NTCP(i));
end
% 
%% Calculate TCP,NTCT,UTCP mean and std

Tmax1=7*14;
TCP_2vec=zeros(Tmax1,8);
NTCPvec=zeros(Tmax1,8);
UTCP_2vec=zeros(Tmax1,8);

Num_sim=200;
for k=1:8
Set_sim=k;
Block_sim=25;
NS=25;
for i=1:Tmax1
    numM_2perc=0;
    numQ=0;
    for j=(Set_sim*Block_sim-NS)+1:Set_sim*Block_sim
        if j~=162
         load (['Results_DoseEx_45day/Evo45day_Num_',num2str(j),'.mat'])
         if tau_M_2perc>0 && tau_M_2perc<=i
             numM_2perc=numM_2perc+1;
         end
         if tau_Q>0 && tau_Q<=i
             numQ=numQ+1;
         end
        end
    end

    TCP_2vec(i,k)=numM_2perc/(Block_sim);
    NTCPvec(i,k)=numQ/(Block_sim);
    UTCP_2vec(i,k)=TCP_2vec(i,k)*(1-NTCPvec(i,k));
end
end

TCP_2vec_ave=mean(TCP_2vec,2);
TCP_2vec_std=(std(TCP_2vec'))';
UTCP_2vec_ave=mean(UTCP_2vec,2);
UTCP_2vec_std=((UTCP_2vec'))';
NTCPvec_ave=mean(NTCPvec,2);
NTCPvec_std=((NTCPvec'))';

save (['Evo45day_AveTCP.mat'],'TCP_2vec_ave','UTCP_2vec_ave','NTCPvec_ave','TCP_2vec_std','UTCP_2vec_std','NTCPvec_std');


%% Calculate volume of tumor and tissue   
Tvol_vec=zeros(Tmax1,1);
Tred_vec=zeros(Tmax1,1);

Tvol_vec_T2=zeros(Tmax1,1);
Tred_vec_T2=zeros(Tmax1,1);

Qvol_vec=zeros(Tmax1,1);
Qred_vec=zeros(Tmax1,1);

Qvol_vec_T2=zeros(Tmax1,1);
Qred_vec_T2=zeros(Tmax1,1);

Tvol_in=[97.5659 85.2739];
n=1;
for j=1:Num_sim
    if j~=162
    load (['Results_DoseEx_45day/Evo45day_Num_',num2str(j),'.mat'])
    Tvol(1)=Tvol_in(1);
    Tvol_T2(1)=Tvol_in(2);

    Tvol_vec=[Tvol_vec Tvol'];
    Tvol_vec_T2=[Tvol_vec_T2 Tvol_T2'];
    
    Tred_vec=[Tred_vec -((Tvol-Tvol_in(1))./Tvol_in(1))'];
    Tred_vec_T2=[Tred_vec_T2 -((Tvol_T2-Tvol_in(2))./Tvol_in(2))'];
    
    Qvol_vec=[Qvol_vec Qvol'];
    Qvol_vec_T2=[Qvol_vec_T2 Qvol_T2'];

    Qred_vec=[Qred_vec -((Qvol-Qvol(1))./Qvol(1))'];
    Qred_vec_T2=[Qred_vec_T2 -((Qvol_T2-Qvol(1))./Qvol(1))'];
    
    timeRec_vec(j)=timeRec_M;

end

Tvol_ave=mean(Tvol_vec(:,2:end)');
Tvol_std=std(Tvol_vec(:,2:end)');
Tvol_max=max(Tvol_vec(:,2:end)');
Tvol_min=min(Tvol_vec(:,2:end)');

Tvol_ave_T2=mean(Tvol_vec_T2(:,2:end)');
Tvol_max_T2=max(Tvol_vec_T2(:,2:end)');
Tvol_min_T2=min(Tvol_vec_T2(:,2:end)');
Tvol_std_T2=std(Tvol_vec_T2(:,2:end)');

Tred_ave=mean(Tred_vec(:,2:end)').*100;
Tred_max=max(Tred_vec(:,2:end)').*100;
Tred_min=min(Tred_vec(:,2:end)').*100;
Tred_std=std(Tred_vec(:,2:end)').*100;

Tred_ave_T2=mean(Tred_vec_T2(:,2:end)').*100;
Tred_max_T2=max(Tred_vec_T2(:,2:end)').*100;
Tred_min_T2=min(Tred_vec_T2(:,2:end)').*100;
Tred_std_T2=std(Tred_vec_T2(:,2:end)').*100;

Qvol_ave=mean(Qvol_vec(:,2:end)');
Qvol_max=max(Qvol_vec(:,2:end)');
Qvol_min=min(Qvol_vec(:,2:end)');
Qvol_std=std(Qvol_vec(:,2:end)');

Qred_ave=mean(Qred_vec(:,2:end)').*100;
Qred_max=max(Qred_vec(:,2:end)').*100;
Qred_min=min(Qred_vec(:,2:end)').*100;
Qred_std=std(Qred_vec(:,2:end)').*100;

Qvol_ave_T2=mean(Qvol_vec_T2(:,2:end)');
Qvol_max_T2=max(Qvol_vec_T2(:,2:end)');
Qvol_min_T2=min(Qvol_vec_T2(:,2:end)');
Qvol_std_T2=std(Qvol_vec_T2(:,2:end)');

Qred_ave_T2=mean(Qred_vec_T2(:,2:end)').*100;
Qred_max_T2=max(Qred_vec_T2(:,2:end)').*100;
Qred_min_T2=min(Qred_vec_T2(:,2:end)').*100;
Qred_std_T2=std(Qred_vec_T2(:,2:end)').*100;

save (['Evo45day_Final.mat'],'timeRec_vec','Tvol_ave','Tvol_min','Tvol_max','Tvol_std','Tvol_min_T2','Tvol_max_T2','Tvol_std_T2','Tvol_ave_T2','Tred_ave','Tred_min','Tred_max','Tred_std','Tred_ave_T2','Tred_min_T2','Tred_max_T2','Tred_std_T2','Qred_min','Qred_max','Qvol_min','Qvol_max','Qvol_ave','Qvol_std','Qred_ave','Qred_std','Qvol_ave_T2','Qvol_min_T2','Qvol_std_T2','Qvol_max_T2','Qred_ave_T2','Qred_min_T2','Qred_max_T2','Qred_std_T2','Tmax','Tmax1');





