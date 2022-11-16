%%% Fitting models to the first recurrence time


%%% Overview over the code:
% 1. Fitting to a part of the 1st recurrence data, PMQ+ vs blood-stage
% 2. Fitting to data grouped by drug and study
% 3. CIs for fitting to data grouped by drug and study


% This script uses the following functions:
% modelfit.m
% modelfit_1rec.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 1. Fitting to a part of the 1st recurrence data, PMQ+ vs blood-stage %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Only use day x[i-1] to x[i] data and ignore drug washout time
t_days = [0;60;120;180;240;300;1000]; % every 60 days & longer last interval

nllh = zeros(length(t_days)-1,1); % negative loglikelihood
par = table('Size',[2,6],'VariableTypes',{'double','double','double','double','double','double'},...
    'VariableNames',{'t60','t120','t180','t240','t300','t1000'},...
    'rowNames',{'r_PMQ','r_blood'});
% model fitting:
n_par = 2; % number of estimated parameters for each model
% parameters: r_pmq, r_blood
k = 100; % number of random initial conditions tried for the fit
opts = optimoptions('fmincon','Display','iter','UseParallel',true);


for l=1:(length(t_days)-1)
    time_min = t_days(l); % only use data from day time_min to time_max
    time_max = t_days(l+1);
    
    % load data:
    combined_time_data = cell2mat(struct2cell(load('Combined_Time_Data.mat')));
    
    % grouping of drugs as "blood" and "PMQ":
    drugs2 = string(combined_time_data.arm_num);
    drugs2(ismember(combined_time_data.arm_num,'CHQ/PMQ')) = repmat('PMQ',sum(ismember(combined_time_data.arm_num,'CHQ/PMQ')),1);
    drugs2(ismember(combined_time_data.arm_num,'AS')) = repmat('blood',sum(ismember(combined_time_data.arm_num,'AS')),1);
    drugs2(ismember(combined_time_data.arm_num,'CHQ')) = repmat('blood',sum(ismember(combined_time_data.arm_num,'CHQ')),1);
    
    % extract data of time of first recurrence for each individual:
    id = unique(string(combined_time_data.patientid));
    time1 = nan(length(id),1);
    event1 = nan(length(id),1);
    drug = strings(length(id),1);
    
    for i=2:length(id)
        time1(i) = combined_time_data.Time_to_event(ismember(combined_time_data.patientid,id{i}) & combined_time_data.episode==2);
        event1(i) = 1-combined_time_data.Censored(ismember(combined_time_data.patientid,id{i}) & combined_time_data.episode==2);
        drug(i) = unique(drugs2(ismember(combined_time_data.patientid,id{i})));
    end
    
    % exclude individual censored at time 0:
    ind = find(time1==0);
    ind = [1:(ind-1),(ind+1):length(time1)];
    data = table(id(ind),time1(ind),event1(ind),drug(ind),'VariableNames',{'id','time','event','drug'});
    
    % exclude individuals with first recurrence before time_min
    data = data(data.time>time_min,:);
    
    % censor individuals with the first recurrence after time_max
    for i=1:size(data,1)
        if(data.time(i)>time_max)
            data.time(i) = time_max;
            data.event(i) = 0;
        end
    end
    
    % shift times s.t. time_min is at 0:
    data.time = data.time-time_min;
    
    times_cen_pmq = data.time(data.event==0 & data.drug=="PMQ");
    times_cen_blood = data.time(data.event==0 & data.drug=="blood");
    times_event_pmq = data.time(data.event==1 & data.drug=="PMQ");
    times_event_blood = data.time(data.event==1 & data.drug=="blood");
    
    % model fitting using the function "modelfit.m" (with numerical ODE
    % solution):
    % this function returns the negative loglikelihood
    fun = @(par) (-1)*(sum(log(modelfit(times_cen_pmq,1,[0,0,par(1)])))+...
        sum(log(modelfit(times_event_pmq-1,1,[0,0,par(1)])-modelfit(times_event_pmq,1,[0,0,par(1)])))+...
        sum(log(modelfit(times_cen_blood,1,[0,0,par(2)])))+...
        sum(log(modelfit(times_event_blood-1,1,[0,0,par(2)])-modelfit(times_event_blood,1,[0,0,par(2)]))));
    
    parameters = zeros(k,n_par);
    nllhs = zeros(k,1);
    for i=1:k
        while nllhs(i,1)==0
            try
                [par_tmp,nllh_tmp] = fmincon(fun,[rand(1),rand(1)],[],[],[],[],[0,0],[Inf,Inf],[],opts);
                parameters(i,:) = par_tmp;
                nllhs(i,1) = nllh_tmp;
                disp(i)
            catch
                disp(".")
            end
        end
    end
    [nllh_min,ind_min] = min(nllhs);
    parameters_min = parameters(ind_min,:);
    
    nllh(l) = nllh_min;
    par{:,l} = parameters_min';
end

% Estimated contribution of relapses to the overall number of recurrences:
contr = [(par{2,1}-par{1,1})/par{2,1},(par{2,2}-par{1,2})/par{2,2},(par{2,3}-par{1,3})/par{2,3},...
    (par{2,4}-par{1,4})/par{2,4},(par{2,5}-par{1,5})/par{2,5},(par{2,6}-par{1,6})/par{2,6}];

% best fit after 100 optimizations with random initial values:
% % 0-60 days
% nllh(1) = 1.5553e+03;
% par{:,1} = [0.00043377;0.013901];
% contr(1) = 0.9688;
% % 60-120 days
% nllh(2) = 486.1646;
% par{:,2} = [0.0004603;0.0074761];
% contr(2) = 0.9384;
% % 120-180 days
% nllh(3) = 276.6285;
% par{:,3} = [0.00047973;0.0036491];
% contr(3) = 0.8685;
% % 180-240 days
% nllh(4) = 162.6269;
% par{:,4} = [0.00042352;0.0011361];
% contr(4) = 0.6272;
% % 240-300 days
% nllh(5) = 182.6010;
% par{:,5} = [0.00057726;0.0013703];
% contr(5) = 0.5787;
% % >300 days
% nllh(6) = 218.8106;
% par{:,6} = [0.00073971;0.0015044];
% contr(6) = 0.5083;


% Visualize fit (Fig. S9):
t1 = 0:60;
U_bl_1 =  modelfit(t1,1,[0,0,par{2,1}]);
U_pmq_1 =  modelfit(t1,1,[0,0,par{1,1}]);
U_bl_2 =  modelfit(t1,1,[0,0,par{2,2}]);
U_pmq_2 =  modelfit(t1,1,[0,0,par{1,2}]);
U_bl_3 =  modelfit(t1,1,[0,0,par{2,3}]);
U_pmq_3 =  modelfit(t1,1,[0,0,par{1,3}]);
U_bl_4 =  modelfit(t1,1,[0,0,par{2,4}]);
U_pmq_4 =  modelfit(t1,1,[0,0,par{1,4}]);
U_bl_5 =  modelfit(t1,1,[0,0,par{2,5}]);
U_pmq_5 =  modelfit(t1,1,[0,0,par{1,5}]);
t6 = 0:100;
U_bl_6 =  modelfit(t6,1,[0,0,par{2,6}]);
U_pmq_6 =  modelfit(t6,1,[0,0,par{1,6}]);

subplot(3,3,1)
plot(t1,U_bl_1,'r','LineWidth',3)
ylim([0 1.05])
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t1,U_pmq_1,'b','LineWidth',3)
hold off

subplot(3,3,2)
plot(t1,U_bl_2,'r','LineWidth',3)
ylim([0 1.05])
xticks([0 20 40 60])
xticklabels({'60','80','100','120'})
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t1,U_pmq_2,'b','LineWidth',3)
hold off

subplot(3,3,3)
plot(t1,U_bl_3,'r','LineWidth',3)
ylim([0 1.05])
xticks([0 20 40 60])
xticklabels({'120','140','160','180'})
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t1,U_pmq_3,'b','LineWidth',3)
hold off

subplot(3,3,4)
plot(t1,U_bl_4,'r','LineWidth',3)
ylim([0 1.05])
xticks([0 20 40 60])
xticklabels({'180','200','220','240'})
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t1,U_pmq_4,'b','LineWidth',3)
hold off

subplot(3,3,5)
plot(t1,U_bl_5,'r','LineWidth',3)
ylim([0 1.05])
xticks([0 20 40 60])
xticklabels({'240','260','280','300'})
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t1,U_pmq_5,'b','LineWidth',3)
hold off

subplot(3,3,6)
plot(t6,U_bl_6,'r','LineWidth',3)
ylim([0 1.05])
xticks([0 20 40 60 80 100])
xticklabels({'300','320','340','360','380','400'})
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t6,U_pmq_6,'b','LineWidth',3)
hold off

subplot(3,3,7:9)
plot([0,60,60,120,120,180,180,240,240,300,300,400],[par{1,1},par{1,1},par{1,2},par{1,2},...
    par{1,3},par{1,3},par{1,4},par{1,4},par{1,5},par{1,5},par{1,6},par{1,6}],'b','LineWidth',3)
xlabel('Time [days]')
ylabel('Recurrence rate [per day]')
hold on
plot([0,60,60,120,120,180,180,240,240,300,300,400],[par{2,1},par{2,1},par{2,2},par{2,2},...
    par{2,3},par{2,3},par{2,4},par{2,4},par{2,5},par{2,5},par{2,6},par{2,6}],'r','LineWidth',3)
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% 2. Fitting to data grouped by drug and study %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% load data:
combined_time_data = cell2mat(struct2cell(load('Combined_Time_Data.mat')));

% grouping of drugs as AS, CHQ, CHQ/PMQ, DP/PMQ
drugs2 = string(combined_time_data.arm_num);
drugs2(ismember(combined_time_data.PMQ_partner,'DP')) = repmat('DP/PMQ',sum(ismember(combined_time_data.PMQ_partner,'DP')),1);

% extract data of time of first & second recurrence for each individual:
id = unique(string(combined_time_data.patientid));
time1 = nan(length(id),1);
event1 = nan(length(id),1);
time2 = nan(length(id),1);
event2 = nan(length(id),1);
drug = strings(length(id),1);
study = strings(length(id),1);
drug_study = strings(length(id),1);
followup = nan(length(id),1);

for i=1:length(id)
    time1(i) = combined_time_data.Time_to_event(ismember(combined_time_data.patientid,id{i}) & combined_time_data.episode==2);
    event1(i) = 1-combined_time_data.Censored(ismember(combined_time_data.patientid,id{i}) & combined_time_data.episode==2);
    if max(combined_time_data.episode(ismember(combined_time_data.patientid,id{i})))>=3
        time2(i) = combined_time_data.Time_to_event(ismember(combined_time_data.patientid,id{i}) & combined_time_data.episode==3);
        event2(i) = 1-combined_time_data.Censored(ismember(combined_time_data.patientid,id{i}) & combined_time_data.episode==3);
    end
    drug(i) = unique(drugs2(ismember(combined_time_data.patientid,id{i})));
    study(i) = string(extractBetween(id{i},1,3));
    drug_study(i) = strcat(drug(i),'_',study(i));
    followup(i) = combined_time_data.FU_time(ismember(combined_time_data.patientid,id{i}) & combined_time_data.episode==2);
end

% exclude individual censored at time 0:
ind = find(time1==0);
ind = [1:(ind-1),(ind+1):length(time1)];
data = table(id(ind),time1(ind),event1(ind),time2(ind),event2(ind),drug(ind),study(ind),drug_study(ind),followup(ind),...
    'VariableNames',{'id','time1','event1','time2','event2','drug','study','drug_study','followup'});


% Model fit:
% fit models 1 to 4 with 2 reinfection rates, daily follow-up, and 10 risk
% groups for models 3 and 4

n = 10; % number of relapse risk groups for model 3
f = 1; % use daily follow-up scheme
reinf = 2; % 2 reinfection rates

n_par = [9,10,10,11]; % number of estimated parameters for each model

k = 100; % number of optimisations
opts = optimoptions('fmincon','Display','iter','UseParallel',true);

for m=1:4
    if m==1
        nllh_1 = zeros(k,1);
        par_1 = zeros(k,n_par(m));
        
        fun =@(par) modelfit_1rec(data,m,f,par,reinf);
        
        for i=1:k
            [par_tmp,nllh_tmp] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10-5,rand(1)*10,...
                rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10],...
                [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,0],...
                [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
            nllh_1(i) = nllh_tmp;
            par_1(i,:) = par_tmp;
            disp(i)
        end
        [nllh_min1,ind_min1] = min(nllh_1);
        par_min1 = par_1(ind_min1,:);
        
    elseif m==2
        nllh_2 = zeros(k,1);
        par_2 = zeros(k,n_par(m));
        
        fun =@(par) modelfit_1rec(data,m,f,par,reinf);
        
        warning off
        for i=1:k
            [par_tmp,nllh_tmp] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10-5,rand(1)*10,...
                rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10],...
                [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,0,0],...
                [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
        end
        [nllh_min2,ind_min2] = min(nllh_2);
        par_min2 = par_2(ind_min2,:);
        warning on
    
    elseif m==3
        nllh_3 = zeros(k,1);
        par_3 = zeros(k,n_par(m));
        
        fun =@(par) modelfit_1rec(data,m,f,[par,n],reinf);
        
        warning off
        for i=1:k
            [par_tmp,nllh_tmp] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10-5,rand(1)*10,...
                rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10],...
                [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,-Inf,0],...
                [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
            nllh_3(i) = nllh_tmp;
            par_3(i,:) = par_tmp;
            disp(i)
        end
        [nllh_min3,ind_min3] = min(nllh_3);
        par_min3 = par_3(ind_min3,:);
        warning on
        
    elseif m==4
        nllh_4 = zeros(k,1);
        par_4 = zeros(k,n_par(m));
        
        fun =@(par) modelfit_1rec(data,m,f,[par,n],reinf);
        
        warning off
        for i=1:k
            [par_tmp,nllh_tmp] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10-5,rand(1)*10,...
                rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10,rand(1)*10],...
                [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,-Inf,0,0],...
                [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
            nllh_4(i) = nllh_tmp;
            par_4(i,:) = par_tmp;
            disp(i)
        end
        [nllh_min4,ind_min4] = min(nllh_4);
        par_min4 = par_4(ind_min4,:);
        warning on
    end
        
end

%%% Visualize model fits with no relapses for PMQ: (Fig. S10, Fig. 2)
 
% Model fit results:
% nllh_min1 = 2.8340e+03;
% par_min1 = [3.2637;1.6906;4.0960;1.3505;3.8848;0.0571;0.0721;0.001063726;0.000611389];
% nllh_min2 = 2.7380e+03;
% par_min2 = [2.8145;0.2562;3.4201;0.3801;3.8855;0.0572;0.000885083;0.000529934;0.0879;0.0287];
% nllh_min3 = 2.7315e+03;
% par_min3 = [3.107285344;0.283339868;3.681248508;0.417569833;3.901824196;0.063138469;0.000784488;...
%     0.00054339;-1.875040576;5.046851789;10];
% nllh_min4 = 2.7313e+03;
% par_min4 = [3.092853108080442;0.270203211961760;3.679732616373312;0.423061058175370;3.874721472437503;...
%     0.062464924763542;0.000830851635662;0.000545916745933;-1.597612718997407;3.934106528974070;0.008126360695056;10];

% Plot model 1:
t = 0:400;
U_as = modelfit(t,1,[par_min1(1),par_min1(2),par_min1(7)]);
U_chq = modelfit(t,1,[par_min1(3),par_min1(4),par_min1(7)]);
U_chq_pmq_vhx = modelfit(t,1,[par_min1(3),par_min1(4),par_min1(8)]);
U_chq_pmq_bpd = modelfit(t,1,[par_min1(3),par_min1(4),par_min1(9)]);
U_dp_pmq_bpd = modelfit(t,1,[par_min1(5),par_min1(6),par_min1(9)]);

plot(t,U_as,'r','LineWidth',3)
ylim([0 1.05])
title('Model 1: constant relapse rate')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_chq,'g','LineWidth',3)
plot(t,U_chq_pmq_vhx,'Color',[0.3010,0.7450,0.9330],'LineWidth',3)
plot(t,U_chq_pmq_bpd,'b','LineWidth',3)
plot(t,U_dp_pmq_bpd,'Color',[0.9290,0.6940,0.1250],'LineWidth',3)
text(250,0.3,['AIC = ',num2str(round(2*nllh_min1+2*9))])
hold off

% Plot models 2 and 3:
t = 0:400;
U_as_2 = modelfit(t,2,[par_min2(1),par_min2(2),par_min2(7),par_min2(9),par_min2(10)]);
U_chq_2 = modelfit(t,2,[par_min2(3),par_min2(4),par_min2(7),par_min2(9),par_min2(10)]);
U_chq_pmq_vhx_2 = modelfit(t,2,[par_min2(3),par_min2(4),par_min2(7),0,par_min2(10)]);
U_chq_pmq_bpd_2 = modelfit(t,2,[par_min2(3),par_min2(4),par_min2(8),0,par_min2(10)]);
U_dp_pmq_bpd_2 = modelfit(t,2,[par_min2(5),par_min2(6),par_min2(8),0,par_min2(10)]);
U_as_3 = sum(modelfit(t,3,[par_min3(1),par_min3(2),par_min3(7),par_min3(11),par_min3(9),par_min3(10)]),1);
U_chq_3 = sum(modelfit(t,3,[par_min3(3),par_min3(4),par_min3(7),par_min3(11),par_min3(9),par_min3(10)]),1);
U_chq_pmq_vhx_3 = sum(modelfit(t,3,[par_min3(3),par_min3(4),par_min3(7),par_min3(11),0,0]),1);
U_chq_pmq_bpd_3 = sum(modelfit(t,3,[par_min3(3),par_min3(4),par_min3(8),par_min3(11),0,0]),1);
U_dp_pmq_bpd_3 = sum(modelfit(t,3,[par_min3(5),par_min3(6),par_min3(8),par_min3(11),0,0]),1);

subplot(1,2,1)
plot(t,U_as_2,'r','LineWidth',3)
ylim([0 1.05])
title('Model 2: temporal heterogeneity')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_chq_2,'g','LineWidth',3)
plot(t,U_chq_pmq_vhx_2,'Color',[0.3010,0.7450,0.9330],'LineWidth',3)
plot(t,U_chq_pmq_bpd_2,'b','LineWidth',3)
plot(t,U_dp_pmq_bpd_2,'Color',[0.9290,0.6940,0.1250],'LineWidth',3)
text(250,0.3,['AIC = ',num2str(round(2*nllh_min2+2*10))])
hold off

subplot(1,2,2)
plot(t,U_as_3,'r','LineWidth',3)
ylim([0 1.05])
title('Model 3: population heterogeneity')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_chq_3,'g','LineWidth',3)
plot(t,U_chq_pmq_vhx_3,'Color',[0.3010,0.7450,0.9330],'LineWidth',3)
plot(t,U_chq_pmq_bpd_3,'b','LineWidth',3)
plot(t,U_dp_pmq_bpd_3,'Color',[0.9290,0.6940,0.1250],'LineWidth',3)
text(250,0.3,['AIC = ',num2str(round(2*nllh_min3+2*10))])
hold off

% Plot model 4:
t = 0:400;
U_as_4 = sum(modelfit(t,4,[par_min4(1),par_min4(2),par_min4(7),par_min4(12),par_min4(9),par_min4(10),par_min4(11)]),1);
U_chq_4 = sum(modelfit(t,4,[par_min4(3),par_min4(4),par_min4(7),par_min4(12),par_min4(9),par_min4(10),par_min4(11)]),1);
U_chq_pmq_vhx_4 = sum(modelfit(t,4,[par_min4(3),par_min4(4),par_min4(7),par_min4(12),0,0,0]),1);
U_chq_pmq_bpd_4 = sum(modelfit(t,4,[par_min4(3),par_min4(4),par_min4(8),par_min4(12),0,0,0]),1);
U_dp_pmq_bpd_4 = sum(modelfit(t,4,[par_min4(5),par_min4(6),par_min4(8),par_min4(12),0,0,0]),1);

plot(t,U_as_4,'r','LineWidth',3)
ylim([0 1.05])
title('Model 4: temporal & population heterogeneity')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_chq_4,'g','LineWidth',3)
plot(t,U_chq_pmq_vhx_4,'Color',[0.3010,0.7450,0.9330],'LineWidth',3)
plot(t,U_chq_pmq_bpd_4,'b','LineWidth',3)
plot(t,U_dp_pmq_bpd_4,'Color',[0.9290,0.6940,0.1250],'LineWidth',3)
text(250,0.3,['AIC = ',num2str(round(2*nllh_min4+2*11))])
hold off




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% 3. CIs for fitting to data grouped by drug and study %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Confidence intervals for parameters fit to first recurrence data
% (Tables S14-17)

% load data:
combined_time_data = cell2mat(struct2cell(load('Combined_Time_Data.mat')));

% grouping of drugs as AS, CHQ, CHQ/PMQ, DP/PMQ
drugs2 = string(combined_time_data.arm_num);
drugs2(ismember(combined_time_data.PMQ_partner,'DP')) = repmat('DP/PMQ',sum(ismember(combined_time_data.PMQ_partner,'DP')),1);

% extract data of time of first & second recurrence for each individual:
id = unique(string(combined_time_data.patientid));
time1 = nan(length(id),1);
event1 = nan(length(id),1);
time2 = nan(length(id),1);
event2 = nan(length(id),1);
drug = strings(length(id),1);
study = strings(length(id),1);
drug_study = strings(length(id),1);
followup = nan(length(id),1);

for i=1:length(id)
    time1(i) = combined_time_data.Time_to_event(ismember(combined_time_data.patientid,id{i}) & combined_time_data.episode==2);
    event1(i) = 1-combined_time_data.Censored(ismember(combined_time_data.patientid,id{i}) & combined_time_data.episode==2);
    if max(combined_time_data.episode(ismember(combined_time_data.patientid,id{i})))>=3
        time2(i) = combined_time_data.Time_to_event(ismember(combined_time_data.patientid,id{i}) & combined_time_data.episode==3);
        event2(i) = 1-combined_time_data.Censored(ismember(combined_time_data.patientid,id{i}) & combined_time_data.episode==3);
    end
    drug(i) = unique(drugs2(ismember(combined_time_data.patientid,id{i})));
    study(i) = string(extractBetween(id{i},1,3));
    drug_study(i) = strcat(drug(i),'_',study(i));
    followup(i) = combined_time_data.FU_time(ismember(combined_time_data.patientid,id{i}) & combined_time_data.episode==2);
end

% exclude individual censored at time 0:
ind = find(time1==0);
ind = [1:(ind-1),(ind+1):length(time1)];
data = table(id(ind),time1(ind),event1(ind),time2(ind),event2(ind),drug(ind),study(ind),drug_study(ind),followup(ind),...
    'VariableNames',{'id','time1','event1','time2','event2','drug','study','drug_study','followup'});

m = 4; % model (1, 2, 3, or 4)
n_par = [9,10,10,11]; % number of estimated parameters for each model
opts = optimoptions('fmincon','Display','iter','UseParallel',true);

n = 10; % number of relapse risk groups for models 3 and 4
f = 1; % use daily follow-up scheme
reinf = 2; % 2 reinfection rates

n_fit = 1000; % number of fits
n_rand = 10; % number of random initial parameter values for each fit

nllhs_all = zeros(n_fit,1);
param_all = zeros(n_fit,n_par(m));
for i=1:n_fit
    ind_rnd = randi(size(data,1),size(data,1),1);
    data_tmp = data(ind_rnd,:); % bootstrapped data
    
    nllh_tmp = zeros(n_rand+1,1);
    par_tmp = zeros(n_rand+1,n_par(m));
    
    if m==1
        fun =@(par) modelfit_1rec(data_tmp,m,f,par,reinf);
        for j=1:(n_rand+1)
            if j==1 % use best fit as initial parameter values
                warning off
                [par_tmp2,nllh_tmp2] = fmincon(fun,[3.2637;1.6906;4.0960;1.3505;3.8848;0.0571;0.0721;0.001063726;0.000611389],...
                    [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,0],...
                    [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            else
                warning off
                [par_tmp2,nllh_tmp2] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10-5,rand(1)*10,...
                    rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10],...
                    [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,0],...
                    [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            end
        end
    elseif m==2
        fun =@(par) modelfit_1rec(data_tmp,m,f,par,reinf);
        for j=1:(n_rand+1)
            if j==1 % use best fit as initial parameter values
                warning off
                [par_tmp2,nllh_tmp2] = fmincon(fun,[2.8145;0.2562;3.4201;0.3801;3.8855;0.0572;0.000885083;0.000529934;0.0879;0.0287],...
                    [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,0,0],...
                    [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            else
                warning off
                [par_tmp2,nllh_tmp2] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10-5,rand(1)*10,...
                    rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10],...
                    [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,0,0],...
                    [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            end
        end
    elseif m==3
        fun =@(par) modelfit_1rec(data_tmp,m,f,[par,n],reinf);
        for j=1:(n_rand+1)
            if j==1 % use best fit as initial parameter values
                %                 warning off
                %                 pctRunOnAll warning off
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[3.107285344,0.283339868,3.681248508,0.417569833,3.901824196,0.063138469,...
                    0.000784488,0.00054339,-1.875040576,5.046851789],...
                    [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,-Inf,0],...
                    [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            else
                %                 warning off
                pctRunOnAll warning off
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10-5,rand(1)*10,...
                    rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10],...
                    [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,-Inf,0],...
                    [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            end
        end
    elseif m==4
        fun =@(par) modelfit_1rec(data,m,f,[par,n],reinf);
        for j=1:(n_rand+1)
            if j==1 % use best fit as initial parameter values
                %                 warning off
                %                 pctRunOnAll warning off
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[3.092853108080442;0.270203211961760;3.679732616373312;...
                    0.423061058175370;3.874721472437503;0.062464924763542;0.000830851635662;0.000545916745933;...
                    -1.597612718997407;3.934106528974070;0.008126360695056]',...
                    [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,-Inf,0,0],...
                    [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            else
                %                 warning off
                pctRunOnAll warning off
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10-5,rand(1)*10,...
                    rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10,rand(1)*10],...
                    [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,-Inf,0,0],...
                    [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            end
        end
    end
    
    [nllh_min,ind_min] = min(nllh_tmp);
    par_min = par_tmp(ind_min,:);
    nllhs_all(i,1) = nllh_min;
    param_all(i,:) = par_min;
    
    disp("Model fits done:")
    disp(i)
end
if m==1
    save('bootstrap_CI_data(11opt)m1(1st-rec).mat','nllhs_all','param_all')
elseif m==2
    save('bootstrap_CI_data(11opt)m2(1st-rec).mat','nllhs_all','param_all')
elseif m==3
    save('bootstrap_CI_data(11opt)m3(1st-rec).mat','nllhs_all','param_all')
elseif m==4
    save('bootstrap_CI_data(11opt)m4(1st-rec).mat','nllhs_all','param_all')
end

% 95% CIs for the different models:
alpha = 95;

for m=1:4 % for each model
    if m==1
        load('bootstrap_CI_data(11opt)m1(1st-rec).mat')
        ci_1 = zeros(size(param_all,2),2);
        for i=1:size(param_all,2)
            ci_1(i,1) = prctile(param_all(:,i),(100-alpha)/2);
            ci_1(i,2) = prctile(param_all(:,i),alpha+(100-alpha)/2);
        end
    elseif m==2
        load('bootstrap_CI_data(11opt)m2(1st-rec).mat')
        ci_2 = zeros(size(param_all,2),2);
        for i=1:size(param_all,2)
            ci_2(i,1) = prctile(param_all(:,i),(100-alpha)/2);
            ci_2(i,2) = prctile(param_all(:,i),alpha+(100-alpha)/2);
        end
    elseif m==3
        load('bootstrap_CI_data(11opt)m3(1st-rec).mat')
        ci_3 = zeros(size(param_all,2),2);
        for i=1:size(param_all,2)
            ci_3(i,1) = prctile(param_all(:,i),(100-alpha)/2);
            ci_3(i,2) = prctile(param_all(:,i),alpha+(100-alpha)/2);
        end
    elseif m==4
        load('bootstrap_CI_data(11opt)m4(1st-rec).mat')
        ci_4 = zeros(size(param_all,2),2);
        for i=1:size(param_all,2)
            ci_4(i,1) = prctile(param_all(:,i),(100-alpha)/2);
            ci_4(i,2) = prctile(param_all(:,i),alpha+(100-alpha)/2);
        end
    end
end


