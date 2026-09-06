%% Tumor Volume - Test 1
clear all
close all
grayColor = [.7 .7 .7];

load('Evo15day_Final.mat')

time=1:1:Tmax1;
Time_to_Rec_ave_vec=[];
Time_to_Rec_ave_vec_min=[];
Time_to_Rec_ave_vec_max=[];

Time_to_Rec_ave_T2_vec=[];
Time_to_Rec_ave_vec_T2_min=[];
Time_to_Rec_ave_vec_T2_max=[];

Tred_endT_ave_vec=[];
Tred_endT_min_vec=[];
Tred_endT_max_vec=[];
Tred_FU_ave_vec=[];
Tred_FU_min_vec=[];
Tred_FU_max_vec=[];

Tred_endT_ave_T2_vec=[];
Tred_endT_min_T2_vec=[];
Tred_endT_max_T2_vec=[];
Tred_FU_ave_T2_vec=[];
Tred_FU_min_T2_vec=[];
Tred_FU_max_T2_vec=[];

Qred_FU_ave_vec=[];
Qred_FU_min_vec=[];
Qred_FU_max_vec=[];

Terapies=[15 20 25 30 35 40 45];
for j=1:size(Terapies,2)

    load (['Evo',num2str(Terapies(j)),'day_Final.mat'])

figure(1)
subplot(1,2,1)
plot(time,Tvol_ave,LineWidth=1)
ylabel('Tumor Volume',Interpreter='latex',FontSize=16)
axis([1 Tmax1 0 140])
xticks([0 20 40 60 80])
yticks([0 20 40 60 80 100 120 140])
xlabel('Time',Interpreter='latex',FontSize=16)
title('No detection threshold',Interpreter='latex',FontSize=14)
axis square 
hold on
subplot(1,2,2)
plot(time,Tvol_ave_T2,LineWidth=1)
axis([1 Tmax1 0 140])
xticks([0 20 40 60 80])
yticks([0 20 40 60 80 100 120 140])
xlabel('Time',Interpreter='latex',FontSize=16)
title('$T_2$ detection threshold',Interpreter='latex',FontSize=14)
axis square 
hold on

if j<5
figure(20)
subplot(2,4,j)
plot(time,Tvol_ave,'k',LineWidth=1.5)
hold on
plot(time, Tvol_max,'--b',LineWidth=1.5)
plot(time,Tvol_min,'--b',LineWidth=1.5)
if j==1 
ylabel('$V_T$',Interpreter='latex',FontSize=16)
end
axis([1 Tmax1 0 140])
xticks([0 40 80])
yticks([0 40 80 120])
title(['$D_T$=',num2str(Terapies(j))],Interpreter='latex',FontSize=14)
axis square 
subplot(2,4,j+4)
plot(time,Tvol_ave_T2,'k',LineWidth=1.5)
hold on
plot(time, Tvol_max_T2,'--b',LineWidth=1.5)
plot(time,Tvol_min_T2,'--b',LineWidth=1.5)
if j==1 
ylabel('$V_T$',Interpreter='latex',FontSize=16)
end
xlabel('Time [d]',Interpreter='latex',FontSize=16)
axis([1 Tmax1 0 140])
axis square 
xticks([0 40 80])
yticks([0 40 80 120])

else
figure(21)
subplot(2,3,j-4)
plot(time,Tvol_ave,'k',LineWidth=1.5)
hold on
plot(time, Tvol_max,'--b',LineWidth=1.5)
plot(time,Tvol_min,'--b',LineWidth=1.5)
if j==5 
ylabel('$V_T$',Interpreter='latex',FontSize=16)
end
axis([1 Tmax1 0 140])
xticks([0 40 80])
yticks([0 40 80 120])
title(['$D_T$=',num2str(Terapies(j))],Interpreter='latex',FontSize=14)
axis square 
subplot(2,3,j-1)
plot(time,Tvol_ave_T2,'k',LineWidth=1.5)
hold on
plot(time, Tvol_max_T2,'--b',LineWidth=1.5)
plot(time,Tvol_min_T2,'--b',LineWidth=1.5)
if j==5 
ylabel('$V_T$',Interpreter='latex',FontSize=16)
end
xlabel('Time [d]',Interpreter='latex',FontSize=16)
axis([1 Tmax1 0 140])
axis square 
xticks([0 40 80])
yticks([0 40 80 120])  

end
Therapyend=Terapies(j)/5*7+1;
Followup=Therapyend+7*4;

Tred_endT_ave_vec=[Tred_endT_ave_vec Tred_ave(Therapyend)];
Tred_endT_min_vec=[Tred_endT_min_vec Tred_min(Therapyend)];
Tred_endT_max_vec=[Tred_endT_max_vec Tred_max(Therapyend)];
Tred_endT_ave_T2_vec=[Tred_endT_ave_T2_vec Tred_ave_T2(Therapyend)];
Tred_endT_min_T2_vec=[Tred_endT_min_T2_vec Tred_min_T2(Therapyend)];
Tred_endT_max_T2_vec=[Tred_endT_max_T2_vec Tred_max_T2(Therapyend)];

Tred_FU_ave_vec=[Tred_FU_ave_vec Tred_ave(Followup)];
Tred_FU_min_vec=[Tred_FU_min_vec Tred_min(Followup)];
Tred_FU_max_vec=[Tred_FU_max_vec Tred_max(Followup)];
Tred_FU_ave_T2_vec=[Tred_FU_ave_T2_vec Tred_ave_T2(Followup)];
Tred_FU_min_T2_vec=[Tred_FU_min_T2_vec Tred_min_T2(Followup)];
Tred_FU_max_T2_vec=[Tred_FU_max_T2_vec Tred_max_T2(Followup)];

Qred_FU_ave_vec=[Qred_FU_ave_vec Qred_ave(Followup)];
Qred_FU_min_vec=[Qred_FU_min_vec Qred_min(Followup)];
Qred_FU_max_vec=[Qred_FU_max_vec Qred_max(Followup)];

end

figure(1)
hold on
legend('15 day','20 day','25 day','30 day','35 day','40 day','45 day')

%% TCP/NTP/UTCP and R-score - Test 2

Terapies=[15 20 25 30 35 40 45];
load('Evo15day_Final.mat')

time=1:1:Tmax1;

for j=1:size(Terapies,2)

    load (['Evo',num2str(Terapies(j)),'day_AveTCP.mat'])

figure(7)
subplot(1,3,1)
plot(time,TCP_2vec_ave,LineWidth=1.5)
hold on
xlabel('Time [d]',Interpreter='latex',FontSize=20)
ylabel('TCP',Interpreter='latex',FontSize=20)
axis([1 45 0 1])
xticks([0 15 30 45])
yticks([0 0.2 .4 .6 .8 1])
axis square 

subplot(1,3,2)
plot(time,NTCPvec_ave,LineWidth=1.5)
hold on
xlabel('Time [d]',Interpreter='latex',FontSize=20)
ylabel('NTCP',Interpreter='latex',FontSize=20)
axis([1 45 0 1])
xticks([0 15 30 45])
yticks([0 0.2 .4 .6 .8 1])
axis square 

subplot(1,3,3)
plot(time,UTCP_2vec_ave,LineWidth=1.5)
hold on
xlabel('Time [d]',Interpreter='latex',FontSize=20)
ylabel('UTCP',Interpreter='latex',FontSize=20)
axis([1 45 0 1])
xticks([0 15 30 45])
yticks([0 0.2 .4 .6 .8 1])
axis square 

figure(18)

for l=1:size(NTCPvec_ave,1)
    if NTCPvec_ave(l)<10^(-3)
        NTCP_Rscore(l)=10^(-3);
    else
        NTCP_Rscore(l)=NTCPvec_ave(l,1);
    end
end
R_score=TCP_2vec_ave(:,1)./((NTCP_Rscore(:,1)).^1);
Rscore_Cum_45day(j)=sum(R_score);

end

max_Rscore=max(Rscore_Cum_45day);
Rscore_Cum_45day_mod=Rscore_Cum_45day./max_Rscore;

figure(18)
bar(Terapies,Rscore_Cum_45day_mod,'FaceColor',[.9 .9 .9],'EdgeColor',[0 0 0],'LineWidth',1.5)
hold on
plot(Terapies,Rscore_Cum_45day_mod,"k*-",MarkerSize=12,LineWidth=1)
xlabel('$D_T$',Interpreter='latex',FontSize=20)
ylabel('R-score',Interpreter='latex',FontSize=20)
axis square 
axis([10 50 0.94 1.0005])
xticks([15 20 25 30 35 40 45])
yticks([.94 .96 .98 1])

figure(7)
hold on
legend('15 day','20 day','25 day','30 day','35 day','40 day','45 day')

%% Recit criteria - Test 3
Tred_endT_NS_vec=[];
Tred_FU_NS_vec=[];
Tred_endT_NS_T2_vec=[];
Tred_FU_NS_T2_vec=[];

Terapies=[15 20 25 30 35 40 45];
for j=1:size(Terapies,2) 
    load (['Evo',num2str(Terapies(j)),'dayNoStoc_Final.mat'])

    Therapyend=Terapies(j)/5*7+1;
    Followup=Therapyend+7*4;

    Tred_endT_NS_vec=[Tred_endT_NS_vec Tred_vec(Therapyend)];
    Tred_endT_NS_T2_vec=[Tred_endT_NS_T2_vec Tred_vec_T2(Therapyend)];
    Tred_FU_NS_vec=[Tred_FU_NS_vec Tred_vec(Followup)];
    Tred_FU_NS_T2_vec=[Tred_FU_NS_T2_vec Tred_vec_T2(Followup)];
end
 
figure(5)
subplot(2,2,1)
plot(Terapies,Tred_endT_ave_vec,'k-*',LineWidth=1)
hold on
plot(Terapies,Tred_endT_NS_vec,'g:*',LineWidth=2)
plot(Terapies,Tred_endT_min_vec,'r--*',LineWidth=1)
plot(Terapies,Tred_endT_max_vec,'r--*',LineWidth=1)
plot([15 45], [100 100],'--','Color', grayColor,LineWidth=1)
plot([15 45], [65 65],'-','Color', grayColor,LineWidth=1)
plot([15 45], [-73 -73],':','Color', grayColor,LineWidth=2)
ylabel('$RV_T(\tau_1)$',Interpreter='latex',FontSize=16)
title('No detection threshold',Interpreter='latex',FontSize=14)
axis([15 45 -150 150])
xticks([15 30 45])
yticks([-150 -100 -50 0 50 100 150])
 axis square 

subplot(2,2,2)
plot(Terapies,Tred_endT_ave_T2_vec,'k-*',LineWidth=1)
hold on
plot(Terapies,Tred_endT_NS_T2_vec,'g:*',LineWidth=2)
plot(Terapies,Tred_endT_min_T2_vec,'r--*',LineWidth=1)
plot(Terapies,Tred_endT_max_T2_vec,'r--*',LineWidth=1)
plot([15 45], [100 100],'--','Color', grayColor,LineWidth=1)
plot([15 45], [65 65],'-','Color', grayColor,LineWidth=1)
plot([15 45], [-73 -73],':','Color', grayColor,LineWidth=2)
title('$T_2$ detection threshold',Interpreter='latex',FontSize=14)
axis([15 45 -150 150])
xticks([15 30 45])
yticks([-150 -100 -50 0 50 100 150])
 axis square 

subplot(2,2,3)
plot(Terapies,Tred_FU_ave_vec,'k-*',LineWidth=1)
hold on
plot(Terapies,Tred_FU_NS_vec,'g:*',LineWidth=2)
plot(Terapies,Tred_FU_min_vec,'r--*',LineWidth=1)
plot(Terapies,Tred_FU_max_vec,'r--*',LineWidth=1)
plot([15 45], [100 100],'--','Color', grayColor,LineWidth=1)
plot([15 45], [65 65],'-','Color', grayColor,LineWidth=1)
plot([15 45], [-73 -73],':','Color', grayColor,LineWidth=2)
xlabel('$D_T$',Interpreter='latex',FontSize=16)
ylabel('$RV_T(\tau_2)$',Interpreter='latex',FontSize=16)
axis([15 45 -150 150])
xticks([15 30 45])
yticks([-150 -100 -50 0 50 100 150])
 axis square 

subplot(2,2,4)
plot(Terapies,Tred_FU_ave_T2_vec,'k-*',LineWidth=1)
hold on
plot(Terapies,Tred_FU_NS_T2_vec,'g:*',LineWidth=2)
plot(Terapies,Tred_FU_min_T2_vec,'r--*',LineWidth=1)
plot(Terapies,Tred_FU_max_T2_vec,'r--*',LineWidth=1)
plot([15 45], [100 100],'--','Color', grayColor,LineWidth=1)
plot([15 45], [65 65],'-','Color', grayColor,LineWidth=1)
plot([15 45], [-73 -73],':','Color', grayColor,LineWidth=2)
legend('Mean','Min','Max','CR','PR','PD','Location','southwest')
xlabel('$D_T$',Interpreter='latex',FontSize=16)
axis([15 45 -150 150])
xticks([15 30 45])
yticks([-150 -100 -50 0 50 100 150])
axis square 


%% Kaplan Meier - Test 4
Terapies = [15 20 25 30 35 40 45];
num_protocols = length(Terapies);
T_max = 7 * 14; 

protocols_tau_matrix = zeros(195, num_protocols);
final_success_rate = zeros(num_protocols, 1);
rmst_values = zeros(num_protocols, 1);

for j = 1:num_protocols
    load(['Evo', num2str(Terapies(j)), 'day_Final.mat'], 'timeRec_vec');
    protocols_tau_matrix(:, j) = timeRec_vec(1:195);
end

figure;
hold on;

for p = 1:num_protocols
    tau_rec = protocols_tau_matrix(:, p);
    
    durations = tau_rec;
    is_censored = (tau_rec == 0);
    durations(is_censored) = T_max; 
    
    [f, x] = ecdf(durations, 'Censoring', is_censored, 'Function', 'survivor');
    final_success_rate(p) = f(end) * 100;
    if x(1) > 0
        x_rmst = [0; x];
        f_rmst = [1; f];
    else
        x_rmst = x;
        f_rmst = f;
    end
    rmst_values(p) = trapz(x_rmst, f_rmst);
    stairs(x, f, 'LineWidth', 1.5,'DisplayName', sprintf('%d day', Terapies(p)));
end

yline(0.5, 'k:', 'LineWidth', 2.0, 'HandleVisibility', 'off');
xlabel('Time [d]', Interpreter='latex', FontSize=14);
ylabel('Relapse-free probability', Interpreter='latex', FontSize=14);
legend('Location', 'northeast', Interpreter='latex', FontSize=14);

xlim([60 T_max]); 
ylim([0 1]);
hold off;
axis square

%% Kaplar-Meier in SM

subplot_positions = [1, 2, 3, 4, 5, 6, 8]; 

for p = 1:num_protocols
    tau_rec = protocols_tau_matrix(:, p);
    
    durations = tau_rec;
    is_censored = (tau_rec == 0);
    durations(is_censored) = T_max; 

    [f, x, flo, fup] = ecdf(durations, 'Censoring', is_censored, 'Function', 'survivor');
    
    if x(1) > 0
        x = [0; x];
        f = [1; f];
        flo = [1; flo];
        fup = [1; fup];
    end

    subplot(3, 3, subplot_positions(p));
    hold on;
   
    valid_idx = ~isnan(flo) & ~isnan(fup);
    x_clean = x(valid_idx);
    f_clean = f(valid_idx);
    flo_clean = flo(valid_idx);
    fup_clean = fup(valid_idx);
    
   
    x_patch = [x_clean; flipud(x_clean)];
    y_patch = [fup_clean; flipud(flo_clean)];
    fill(x_patch, y_patch, [0.75 0.75 0.75], 'EdgeColor', 'none', 'FaceAlpha', 0.35, 'HandleVisibility', 'off');
    

    stairs(x_clean, f_clean, 'LineWidth', 1.8, 'Color', [0.1 0.1 0.1],'DisplayName', sprintf('$D_T = %d$ d', Terapies(p)));
    yline(0.5, 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    
    title(sprintf('$D_T = %d$ days', Terapies(p)), Interpreter='latex', FontSize=20);
    xlabel('Time [d]', Interpreter='latex', FontSize=18);
    ylabel('Relapse-free prob.', Interpreter='latex', FontSize=18);
    
    xlim([60 T_max]);
    ylim([0 1.05]);
    axis square;
    
    hold off;
end

%% TCP/NTCP/UTCP with mean and std - Test 2 SM

time=1:1:Tmax1;

figure(10)
Terapies=[15 20 25 30 35 40 45];
for j=1:size(Terapies,2)

    load (['Evo',num2str(Terapies(j)),'day_AveTCP.mat'])
   if j<7
   subplot(3,3,j)
   errorbar(time,TCP_2vec_ave,TCP_2vec_std,'k')
   ylabel('TCP',Interpreter='latex',FontSize=20)
   xlabel('Time [d]',Interpreter='latex',FontSize=20)
   xticks([0 15 30 45])
   yticks([0 0.2 .4 .6 .8 1])
   title(['$D_T$=',num2str(Terapies(j))],Interpreter='latex',FontSize=20)
   axis([1 45 0 1])
   axis square 
   else 
       subplot(3,3,8)
   errorbar(time,TCP_2vec_ave,TCP_2vec_std,'k')
   ylabel('TCP',Interpreter='latex',FontSize=20)
   xlabel('Time [d]',Interpreter='latex',FontSize=20)
   xticks([0 15 30 45])
   yticks([0 0.2 .4 .6 .8 1])
   title(['$D_T=$',num2str(Terapies(j))],Interpreter='latex',FontSize=20)
   axis([1 45 0 1])
   axis square
   end
end

