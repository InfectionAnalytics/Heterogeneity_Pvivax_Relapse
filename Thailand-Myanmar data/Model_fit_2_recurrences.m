%%% Fitting models to the first and second recurrence time


%%% Overview over the code:
% 1. Fitting to the 1st and 2nd recurrence data
% 2. Confidence intervals via bootstrapping
% 3. Model fit results
% 4. Different plots


% This script uses the following functions:
% modelfit.m, modelfit_2rec.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 1. Fitting to the 1st and 2nd recurrence data %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% parameters for fitting the models:
reinf = 2; % number of different reinfection rates: 1 -> same reinfection rates for both studies,
% 2 -> different reinfection rate for the two different studies

% model fits
nllh = zeros(4,1); % nllh for each model
if reinf==1
    par = table('Size',[14,4],'VariableTypes',{'double','double','double','double'},...
        'VariableNames',{'model1','model2','model3','model4'},...
        'rowNames',{'washout1_AS','washout2_AS','washout1_CHQ','washout2_CHQ','washout1_DP',...
        'washout2_DP','recurrence_AS/CHQ','recurrence_PMQ',...
        'reinfection','relapse_risk1_AS/CHQ','relapse_risk2_AS/CHQ','number_relapse_groups',...
        'I_AS/CHQ','d'});
    n_par = [8,9,9,10]; % number of estimated parameters for each model
elseif reinf==2
    par = table('Size',[16,4],'VariableTypes',{'double','double','double','double'},...
        'VariableNames',{'model1','model2','model3','model4'},...
        'rowNames',{'washout1_AS','washout2_AS','washout1_CHQ','washout2_CHQ','washout1_DP',...
        'washout2_DP','recurrence_AS/CHQ_VHX','recurrence_PMQ_VHX','recurrence_PMQ_BPD',...
        'reinfection_VHX','reinfection_BPD','relapse_risk1_AS/CHQ','relapse_risk2_AS/CHQ',...
        'number_relapse_groups','I_AS/CHQ','d'});
    n_par = [9,10,10,11]; % number of estimated parameters for each model
end

k = 100; % number of random initial conditions tried for the fit
n = 10; % number of relapse risk groups in models 3 & 4
f = 1; % follow-up scheme used for the fit
opts = optimoptions('fmincon','Display','iter','UseParallel',true);

for m=1:4
    
    if m==1 % constant relapse rate
        fun =@(par) modelfit_2rec(data,m,f,par,reinf);
        parameters1 = zeros(k,n_par(m));
        nllhs1 = zeros(k,1);
        for i=1:k
            if reinf==1
                [par_tmp,nllh_tmp] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10-5,rand(1)*10,...
                    rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10],[],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0],...
                    [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
            elseif reinf==2
                [par_tmp,nllh_tmp] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10-5,rand(1)*10,...
                    rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10],[],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,0],...
                    [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
            end
            parameters1(i,:) = par_tmp;
            nllhs1(i,1) = nllh_tmp;
            disp(i)
        end
        [nllh_min1,ind_min1] = min(nllhs1);
        parameters_min1 = parameters1(ind_min1,:);
        
        % try best fit from time to first recurrence fit:
        if reinf==1
            [par_tmp,nllh_tmp] = fmincon(fun,[3.2637;1.6906;4.0960;1.3505;3.8848;0.0571;0.0721;mean([0.001063726;0.000611389])],[],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,0],...
                [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
        elseif reinf==2
            [par_tmp,nllh_tmp] = fmincon(fun,[3.2637;1.6906;4.0960;1.3505;3.8848;0.0571;0.0721;0.001063726;0.000611389],[],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,0],...
                [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
        end
        if nllh_tmp<nllh_min1
            nllh(m,1) = nllh_tmp;
            if reinf==1
                par{1:8,m} = par_tmp';
            elseif reinf==2
                par{1:9,m} = par_tmp';
            end
        else
            nllh(m,1) = nllh_min1;
            if reinf==1
                par{1:8,m} = parameters_min1';
            elseif reinf==2
                par{1:9,m} = parameters_min1';
            end
        end
        
    elseif m==2 % temporal heterogeneity
        fun =@(par) modelfit_2rec(data,m,f,par,reinf);
        parameters2 = zeros(k,n_par(m));
        nllhs2 = zeros(k,1);
        for i=1:k
            if reinf==1
                [par_tmp,nllh_tmp] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10-5,rand(1)*10,...
                    rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10],[],[],[],[],...
                    [-Inf,0,-Inf,0,-Inf,0,0,0,0],[Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
            elseif reinf==2
                [par_tmp,nllh_tmp] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10-5,rand(1)*10,...
                    rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10],[],[],[],[],...
                    [-Inf,0,-Inf,0,-Inf,0,0,0,0,0],[Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
            end
            parameters2(i,:) = par_tmp;
            nllhs2(i,1) = nllh_tmp;
        end
        [nllh_min2,ind_min2] = min(nllhs2);
        parameters_min2 = parameters3(ind_min2,:);
        
        % try best fit from time to first recurrence fit:
        if reinf==1
            [par_tmp,nllh_tmp] = fmincon(fun,[2.8145;0.2562;3.4201;0.3801;3.8855;0.0572;mean([0.000885083;0.000529934]);...
                0.0879;0.0287],[],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,0],...
                [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
        elseif reinf==2
            [par_tmp,nllh_tmp] = fmincon(fun,[2.8145;0.2562;3.4201;0.3801;3.8855;0.0572;...
                0.000885083;0.000529934;0.0879;0.0287],[],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,0,0],...
                [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
        end
        if nllh_tmp<nllh_min2
            nllh(m,1) = nllh_tmp;
            if reinf==1
                par{[1:6,9,13:14],m} = par_tmp;
            elseif reinf==2
                par{[1:6,10,11,15:16],m} = par_tmp;
            end
        else
            nllh(m,1) = nllh_min2;
            if reinf==1
                par{[1:6,9,13:14],m} = parameters_min2';
            elseif reinf==2
                par{[1:6,10,11,15:16],m} = parameters_min2';
            end
        end
        
    elseif m==3 % population heterogeneity
        rel_type = "lognormal"; % relapse risk distribution: "lognormal","gamma" or "exp"
        if reinf==1
            fun =@(par) modelfit_2rec(data,m,f,[par(1:7),n,par(8:9)],reinf,rel_type);
        elseif reinf==2
            fun =@(par) modelfit_2rec(data,m,f,[par(1:8),n,par(9:10)],reinf,rel_type);
        end
        parameters3 = zeros(k,n_par(m));
        nllhs3 = zeros(k,1);
        for i=1:k
            warning off
            while nllhs3(i,1)==0
                try
                    if reinf==1
                        [par_tmp,nllh_tmp] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10-5,rand(1)*10,...
                            rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10],...
                            [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,-Inf,0],...
                            [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                    elseif reinf==2
                        [par_tmp,nllh_tmp] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10-5,rand(1)*10,...
                            rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10],...
                            [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,-Inf,0],...
                            [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                    end
                    parameters3(i,:) = par_tmp;
                    nllhs3(i,1) = nllh_tmp;
                catch
                    disp(".")
                end
            end
            warning on
        end
        [nllh_min3,ind_min3] = min(nllhs3);
        parameters_min3 = parameters3(ind_min3,:);
        
        % try best fit from time to first recurrence fit:
        if reinf==1
            [par_tmp,nllh_tmp] = fmincon(fun,[3.107285344;0.283339868;3.681248508;0.417569833;3.901824196;0.063138469;...
                mean([0.000784488;0.00054339]);-1.875040576;5.046851789]',...
                [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,-Inf,0],...
                [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
        elseif reinf==2
            [par_tmp,nllh_tmp] = fmincon(fun,[3.107285344;0.283339868;3.681248508;0.417569833;3.901824196;0.063138469;...
                0.000784488;0.00054339;-1.875040576;5.046851789]',...
                [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,-Inf,0],...
                [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
        end
        if nllh_tmp<nllh_min3
            nllh(m,1) = nllh_tmp;
            if reinf==1
                par{[1:6,9:12],m} = [par_tmp';n];
            elseif reinf==2
                par{[1:6,10:14],m} = [par_tmp';n];
            end
        else
            nllh(m,1) = nllh_min3;
            if reinf==1
                par{[1:6,9:12],m} = [parameters_min3';n];
            elseif reinf==2
                par{[1:6,10:14],m} = [parameters_min3';n];
            end
        end
        
        
    elseif m==4 % temporal & population heterogeneity
        if reinf==1
            fun =@(par) modelfit_2rec(data,m,f,[par(1:7),n,par(8:10)],reinf);
        elseif reinf==2
            fun =@(par) modelfit_2rec(data,m,f,[par(1:8),n,par(9:11)],reinf);
        end
        parameters4 = zeros(k,n_par(m));
        nllhs4 = zeros(k,1);
        for i=1:k
            warning off
            while nllhs4(i,1)==0
                try
                    if reinf==1
                        [par_tmp,nllh_tmp] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10-5,rand(1)*10,...
                            rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10,rand(1)*10],...
                            [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,-Inf,0,0],...
                            [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                    elseif reinf==2
                        [par_tmp,nllh_tmp] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10-5,rand(1)*10,...
                            rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10,rand(1)*10],...
                            [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,-Inf,0,0],...
                            [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                    end
                    parameters4(i,:) = par_tmp;
                    nllhs4(i,1) = nllh_tmp;
                catch
                    disp(".")
                end
            end
            warning on
        end
        [nllh_min4,ind_min4] = min(nllhs4);
        parameters_min4 = parameters4(ind_min4,:);
        
        % try best fit from time to first recurrence fit:
        if reinf==1
            [par_tmp,nllh_tmp] = fmincon(fun,[3.092853108080442;0.270203211961760;3.679732616373312;0.423061058175370;3.874721472437503;...
                0.062464924763542;mean([0.000830851635662,0.000545916745933]);-1.597612718997407;3.934106528974070;0.008126360695056]',...
                [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,-Inf,0,0],...
                [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
        elseif reinf==2
            [par_tmp,nllh_tmp] = fmincon(fun,[3.092853108080442;0.270203211961760;3.679732616373312;0.423061058175370;3.874721472437503;...
                0.062464924763542;0.000830851635662;0.000545916745933;-1.597612718997407;3.934106528974070;0.008126360695056]',...
                [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,-Inf,0,0],...
                [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
        end
        if nllh_tmp<nllh_min4
            nllh(m,1) = nllh_tmp;
            if reinf==1
                par{[1:6,9:12,14],m} = [par_tmp';n];
            elseif reinf==2
                par{[1:6,10:14,16],m} = [par_tmp';n];
            end
        else
            nllh(m,1) = nllh_min4;
            if reinf==1
                par{[1:6,9:12,14],m} = [parameters_min4';n];
            elseif reinf==2
                par{[1:6,10:13,16,14],m} = [parameters_min4';n];
            end
        end
    end
end

% Model fit results: (see also below)
% nllh(1,1) = 4.3881e+03;
% par{1:9,1} = [3.44260377;1.61060976;4.051098296;1.15040291;3.88587443;0.05723774;...
%     0.08141406;0.001056392;0.00064541822];
% nllh(2,1) = 4.2362e+03;
% par{[1:6,10,11,15:16],2} = [2.7914;0.23686;3.3847;0.3071;3.882;0.056692;0.00097171;0.00056813;0.070277;0.024911];
% nllh(3,1) = 4.2030e+03;
% par{[1:6,10:14],3} = [2.95051812;0.24074922;3.53368025;0.32837685;3.88237672;...
%     0.05992895;0.00082751;0.00057377;-3.72783179;2.80439599;10];
% nllh(4,1) = 4.191225315140548e+03;
% par{[1:6,10:13,16,14],4} = [2.926849386154156;0.238603871853848;3.536214055784965;0.334375576933973;...
%     -10.226838072357639;26.619525370717724;0.000892369002219;0.000618753470553;-3.007393884117966;...
%     1.855962389046954;0.01213258809181310;10];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% 2. Confidence intervals via bootstrapping %%%%%%%%%%%%%%%%%
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


n = 10; % number of relapse risk groups for model 3
f = 1; % use daily follow-up scheme
reinf = 2; % 2 reinfection rates

m = 4; % model
n_par = [9,10,10,11]; % number of estimated parameters for each model
opts = optimoptions('fmincon','Display','iter','UseParallel',true);

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
        fun =@(par) modelfit_2rec(data_tmp,m,f,par,reinf);
        for j=1:(n_rand+1)
            if j==1 % use best fit as initial parameter values
                warning off
                [par_tmp2,nllh_tmp2] = fmincon(fun,[3.44260377;1.61060976;4.051098296;1.15040291;3.88587443;...
                    0.05723774;0.08141406;0.001056392;0.00064541822],...
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
        fun =@(par) modelfit_2rec(data_tmp,m,f,par,reinf);
        for j=1:(n_rand+1)
            if j==1 % use best fit as initial parameter values
                warning off
                [par_tmp2,nllh_tmp2] = fmincon(fun,[2.7914;0.23686;3.3847;0.3071;3.882;0.056692;...
                    0.00097171;0.00056813;0.070277;0.024911],...
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
        fun =@(par) modelfit_2rec(data,m,f,[par(1:8),n,par(9:10)],reinf);
        for j=1:(n_rand+1)
            if j==1 % use best fit as initial parameter values
                %                 pctRunOnAll warning off
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[2.95051812;0.24074922;3.53368025;0.32837685;...
                    3.88237672;0.05992895;0.00082751;0.00057377;-3.72783179;2.80439599]',...
                    [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,-Inf,0],...
                    [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            else
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
        fun =@(par) modelfit_2rec(data,m,f,[par(1:8),n,par(9:11)],reinf);
        for j=1:(n_rand+1)
            if j==1 % use best fit as initial parameter values
                %                 pctRunOnAll warning off
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[2.926849386154156;0.238603871853848;...
                    3.536214055784965;0.334375576933973;-10.226838072357639;26.619525370717724;...
                    0.000892369002219;0.000618753470553;-3.007393884117966;1.855962389046954;...
                    0.01213258809181310]',...
                    [],[],[],[],[-Inf,0,-Inf,0,-Inf,0,0,0,-Inf,0,0],...
                    [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            else
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
    save('bootstrap_CI_data(11opt)m1.mat','nllhs_all','param_all')
elseif m==2
    save('bootstrap_CI_data(11opt)m2.mat','nllhs_all','param_all')
elseif m==3
    save('bootstrap_CI_data(11opt)m3.mat','nllhs_all','param_all')
elseif m==4
    save('bootstrap_CI_data(11opt)m4.mat','nllhs_all','param_all')
end

% 95% CIs for the different models:
alpha = 95;

for m=1:4 % for each model
    if m==1
        load('bootstrap_CI_data(11opt)m1.mat')
        ci_1 = zeros(size(param_all,2),2);
        for i=1:size(param_all,2)
            ci_1(i,1) = prctile(param_all(:,i),(100-alpha)/2);
            ci_1(i,2) = prctile(param_all(:,i),alpha+(100-alpha)/2);
        end
    elseif m==2
        load('bootstrap_CI_data(11opt)m2.mat')
        ci_2 = zeros(size(param_all,2),2);
        for i=1:size(param_all,2)
            ci_2(i,1) = prctile(param_all(:,i),(100-alpha)/2);
            ci_2(i,2) = prctile(param_all(:,i),alpha+(100-alpha)/2);
        end
    elseif m==3
        load('bootstrap_CI_data(11opt)m3.mat')
        ci_3 = zeros(size(param_all,2),2);
        for i=1:size(param_all,2)
            ci_3(i,1) = prctile(param_all(:,i),(100-alpha)/2);
            ci_3(i,2) = prctile(param_all(:,i),alpha+(100-alpha)/2);
        end
    elseif m==4
        load('bootstrap_CI_data(11opt)m4.mat')
        ci_4 = zeros(size(param_all,2),2);
        for i=1:size(param_all,2)
            ci_4(i,1) = prctile(param_all(:,i),(100-alpha)/2);
            ci_4(i,2) = prctile(param_all(:,i),alpha+(100-alpha)/2);
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% 3. Model fit results %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Model fit results using the above code (1.) for different follow-up
% schemes, 1 or 2 reinfection rates, and different numbers of relapse risk
% groups for model 3
% All model fits with k=100 (100 random initial conditions)

%%% Model 1: constant relapse rate
% reinf=1; f=1;
nllh(1,1) = 4.3913e+03;
par{1:8,1} = [3.4391297;1.61612896;4.04232847;1.14716182;3.89350259;0.05831517;0.0809168;0.00073168];

% reinf=2; f=1;
nllh(1,1) = 4.3881e+03;
par{1:9,1} = [3.44260377;1.61060976;4.051098296;1.15040291;3.88587443;0.05723774;0.08141406;0.001056392;0.00064541822];

% reinf=2; f=2;
nllh(1,1) = 2.9236e+03;
par{1:9,1} = [3.25528984;1.83093938;3.977038737;1.22871697;3.82570173;0.04870433;0.082020754;0.001056794;0.000640362];

% reinf=2; f=3;
nllh(1,1) = 2.3843e+03;
par{1:9,1} = [2.88256452;2.30270699;3.87114951;1.34928499;3.7143667;0.03271831;0.08110879;0.001046854;0.000630766];

% reinf=2; f=4;
nllh(1,1) = 1.8472e+03;
par{1:9,1} = [1.54875928;4.10050329;3.64224885;1.64977028;3.52952104;0.02078084;0.0942147;0.001036436;0.000619133];

% reinf=2; f=5;
nllh(1,1) = 2.1356e+03;
par{1:9,1} = [2.7744;2.3836;3.7592;1.4554;3.5353;0.070517;0.083922;0.0010353;0.00061911];

% reinf=2; f=6;
nllh(1,1) = 2.2920e+03;
par{1:9,1} = [2.94936003;2.16422512;3.810031807;1.392310901;3.540521788;0.02519966;0.082028755;0.001038405;0.000619971];

% reinf=2; f=7;
nllh(1,1) = 2.5024e+03;
par{1:9,1} = [3.07594585;2.01133962;3.853021409;1.34118233;3.53641982;0.069025375;0.07991755;0.001041173;0.000620325];


%%% Model 2: temporal heterogeneity
% reinf=1; f=1;
nllh(2,1) = 4.2405e+03;
par{[1:6,9,13:14],2} = [2.7814;0.23586;3.3785;0.31353;3.8893;0.057723;0.00066443;0.06531;0.022986];

% reinf=2; f=1;
nllh(2,1) = 4.2362e+03;
par{[1:6,10,11,15:16],2} = [2.7914;0.23686;3.3847;0.3071;3.882;0.056692;0.00097171;0.00056813;0.070277;0.024911];

% reinf=2; f=2;
nllh(2,1) = 2.7949e+03;
par{[1:6,10,11,15:16],2} = [2.5860;0.2290;3.2696;0.3272;3.8220;0.0482;0.000956467;0.00056354;0.0645;0.0249];

% reinf=2; f=3;
nllh(2,1) = 2.2879e+03;
par{[1:6,10,11,15:16],2} = [2.2770;0.2223;3.1155;0.3194;3.704856984;0.03132212197;0.000946449;0.000553096;0.0584;0.0250];

% reinf=2; f=4;
nllh(2,1) = 1.8164e+03;
par{[1:6,10,11,15:16],2} = [-28.8705;22.6540;2.7168;0.2372;3.4457;0.0577;0.000888602;0.000536003;0.0438;0.0213];

% reinf=2; f=5;
nllh(2,1) = 2.0100e+03;
par{[1:6,10,11,15:16],2} = [2.1242;0.0959;2.9930;0.2305;3.4428;0.1889;0.000973178;0.000541365;0.0639;0.0285];

% reinf=2; f=6;
nllh(2,1) = 2.1493e+03;
par{[1:6,10,11,15:16],2} = [2.4012;0.1367;3.0927;0.2538;3.4861;0.1482;0.001015644;0.000545369;0.0751;0.0304];

% reinf=2; f=7;
nllh(2,1) = 2.3505e+03;
par{[1:6,10,11,15:16],2} = [2.5713;0.1688;3.1755;0.2739;3.4600;0.0482;0.001025447;0.000545646;0.0829;0.0310];


%%% Model 3: population heterogeneity
% reinf=1; n=10; f=1;
nllh(3,1) = 4.2048e+03;
par{[1:6,9:12],3} = [2.9447;0.23984;3.5277;0.329;3.8791;0.059409;0.00063404;-3.7175;2.7055;10];

% reinf=2; n=10; f=1;
nllh(3,1) = 4.2030e+03;
par{[1:6,10:14],3} = [2.95051812;0.24074922;3.53368025;0.32837685;3.88237672;0.05992895;0.00082751;0.00057377;-3.72783179;2.80439599;10];

% reinf=2; n=10; f=2;
nllh(3,1) = 2.7620e+03;
par{[1:6,10:14],3} = [2.7681;0.2368;3.4288;0.3417;3.8176;0.0503;0.000811145;0.000568073;-3.7833;2.7445;10];

% reinf=2; n=10; f=3;
nllh(3,1) = 2.2562e+03;
par{[1:6,10:14],3} = [2.5477;0.2227;3.3021;0.3394;3.7045;0.0331;0.0008028323;0.000559069;-3.8062;2.7774;10];

% reinf=2; n=10; f=4;
nllh(3,1) = 1.7929e+03;
par{[1:6,10:14],3} = [2.009593696;0.16045247;2.993471699;0.2963044078;3.445498662;0.02797017;0.000785771;0.000540561;-3.835265919;2.856245202;10];

% reinf=2; n=10; f=5;
nllh(3,1) = 1.9850e+03;
par{[1:6,10:14],3} = [2.2650;0.1284;3.1318;0.2226;3.4564;0.0633;0.000778126;0.000544682;-3.8855;2.6541;10];

% reinf=2; n=10; f=6;
nllh(3,1) = 2.1219e+03;
par{[1:6,10:14],3} = [2.5167;0.1360;3.2231;0.2511;-20.9493;39.7260;0.000792029;0.000597037;-3.8185;2.7628;10];

% reinf=2; n=10; f=7;
nllh(3,1) = 2.3217e+03;
par{[1:6,10:14],3} = [2.6949;0.1585;3.3016;0.2753;-10.2251;24.5367;0.0008042049;0.0006039908;-3.7529;2.8629;10];

% reinf=2; n=5; f=1;
nllh(3,1) = 4.2028e+03;
par{[1:6,10:14],3} = [2.9639;0.2429;3.5502;0.3361;3.8775;0.0673;0.000838419;0.000575624;-3.6376;3.0019;5];

% reinf=2; n=15; f=1;
nllh(3,1) = 4.2030e+03;
par{[1:6,10:14],3} = [2.9484;0.2388;3.5306;0.3267;3.8794;0.0607;0.000822403;0.000574456;-3.7413;2.7767;15];

% reinf=2; n=20; f=1;
nllh(3,1) = 4.2030e+03;
par{[1:6,10:14],3} = [2.9482;0.2390;3.5315;0.3271;3.8830;0.0611;0.00082275;0.00057491;-3.7423;2.7718;20];

% Gamma distribution of relapse risks:
% reinf=2; n=10; f=1;
nllh(3,1) = 4.2041e+03;
par{[1:6,10:14],3} = [2.903686703943764,0.241085890062690,3.513597074580805,0.338705364820449,3.875813718206521,...
    0.058915640688840,0.000903125587265,0.000578465253674,0.303098845785234,0.384639708871805,10];
ex_min2_gam = 2;

% Exponential distribution of relapse risks:
% reinf=2; n=10; f=1;
nllh(3,1) = 4.2567e+03;
par{[1:6,10:12,14],3} = [2.718041005749698,0.220700541966448,3.338564771709158,0.311841142475142,3.886725287979828,...
    0.060571307252208,0.000642750709578,0.000567349056870,0.035443379089622,10];
ex_min2_exp = 2;


%%% Model 4: temporal & population heterogeneity
% reinf=2; n=10; f=1;
nllh(4,1) = 4.191225315140548e+03;
par{[1:6,10:14,16],4} = [2.926849386154156;0.238603871853848;3.536214055784965;0.334375576933973;...
    -10.226838072357639;26.619525370717724;0.000892369002219;0.000618753470553;-3.007393884117966;...
    1.855962389046954;10;0.01213258809181310];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4. Different plots %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot model fit to time to first recurrence and time to second recurrence
% for reinf = 2, f=1, n=10 (Fig. S11)

% Parameter values:
nllh_min1 = 4.3881e+03;
par_min1 = [3.44260377;1.61060976;4.051098296;1.15040291;3.88587443;0.05723774;0.08141406;0.001056392;0.00064541822];
nllh_min2 = 4.2362e+03;
par_min2 = [2.7914;0.23686;3.3847;0.3071;3.882;0.056692;0.00097171;0.00056813;0.070277;0.024911];
nllh_min3 = 4.2030e+03;
par_min3 = [2.95051812;0.24074922;3.53368025;0.32837685;3.88237672;0.05992895;0.00082751;0.00057377;-3.72783179;2.80439599;10];
nllh_min4 = 4.191225315140548e+03;
par_min4 = [2.926849386154156;0.238603871853848;3.536214055784965;0.334375576933973;...
    -10.226838072357639;26.619525370717724;0.000892369002219;0.000618753470553;-3.007393884117966;...
    1.855962389046954;10;0.01213258809181310];

t = 0:400;
U_as_1 = modelfit(t,1,[par_min1(1),par_min1(2),par_min1(7)]);
U_chq_1 = modelfit(t,1,[par_min1(3),par_min1(4),par_min1(7)]);
U_chq_pmq_vhx_1 = modelfit(t,1,[par_min1(3),par_min1(4),par_min1(8)]);
U_chq_pmq_bpd_1 = modelfit(t,1,[par_min1(3),par_min1(4),par_min1(9)]);
U_dp_pmq_bpd_1 = modelfit(t,1,[par_min1(5),par_min1(6),par_min1(9)]);
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
U_as_4 = sum(modelfit(t,4,[par_min4(1),par_min4(2),par_min4(7),par_min4(11),par_min4(9),par_min4(10),par_min4(12)]),1);
U_chq_4 = sum(modelfit(t,4,[par_min4(3),par_min4(4),par_min4(7),par_min4(11),par_min4(9),par_min4(10),par_min4(12)]),1);
U_chq_pmq_vhx_4 = sum(modelfit(t,4,[par_min4(3),par_min4(4),par_min4(7),par_min4(11),0,0,0]),1);
U_chq_pmq_bpd_4 = sum(modelfit(t,4,[par_min4(3),par_min4(4),par_min4(8),par_min4(11),0,0,0]),1);
U_dp_pmq_bpd_4 = sum(modelfit(t,4,[par_min4(5),par_min4(6),par_min4(8),par_min4(11),0,0,0]),1);

subplot(4,2,1)
plot(t,U_as_1,'r','LineWidth',3)
ylim([0 1.05])
title('Model 1: 1st rec.')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_chq_1,'g','LineWidth',3)
plot(t,U_chq_pmq_vhx_1,'Color',[0.3010,0.7450,0.9330],'LineWidth',3)
plot(t,U_chq_pmq_bpd_1,'b','LineWidth',3)
plot(t,U_dp_pmq_bpd_1,'Color',[0.9290,0.6940,0.1250],'LineWidth',3)
text(250,0.3,['AIC = ',num2str(round(2*nllh_min1+2*9))])
hold off

subplot(4,2,2)
plot(t,U_as_1,'r','LineWidth',3)
ylim([0 1.05])
title('Model 1: 2nd rec.')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_chq_1,'g','LineWidth',3)
plot(t,U_chq_pmq_vhx_1,'Color',[0.3010,0.7450,0.9330],'LineWidth',3)
plot(t,U_chq_pmq_bpd_1,'b','LineWidth',3)
plot(t,U_dp_pmq_bpd_1,'Color',[0.9290,0.6940,0.1250],'LineWidth',3)
text(250,0.3,['AIC = ',num2str(round(2*nllh_min1+2*9))])
hold off

subplot(4,2,3)
plot(t,U_as_2,'r','LineWidth',3)
ylim([0 1.05])
title('Model 2: 1st rec.')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_chq_2,'g','LineWidth',3)
plot(t,U_chq_pmq_vhx_2,'Color',[0.3010,0.7450,0.9330],'LineWidth',3)
plot(t,U_chq_pmq_bpd_2,'b','LineWidth',3)
plot(t,U_dp_pmq_bpd_2,'Color',[0.9290,0.6940,0.1250],'LineWidth',3)
text(250,0.3,['AIC = ',num2str(round(2*nllh_min2+2*10))])
hold off

subplot(4,2,4)
plot(t,U_as_2,'r','LineWidth',3)
ylim([0 1.05])
title('Model 2: 2nd rec.')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_chq_2,'g','LineWidth',3)
plot(t,U_chq_pmq_vhx_2,'Color',[0.3010,0.7450,0.9330],'LineWidth',3)
plot(t,U_chq_pmq_bpd_2,'b','LineWidth',3)
plot(t,U_dp_pmq_bpd_2,'Color',[0.9290,0.6940,0.1250],'LineWidth',3)
text(250,0.3,['AIC = ',num2str(round(2*nllh_min2+2*10))])
hold off

subplot(4,2,5)
plot(t,U_as_3,'r','LineWidth',3)
ylim([0 1.05])
title('Model 3: 1st rec.')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_chq_3,'g','LineWidth',3)
plot(t,U_chq_pmq_vhx_3,'Color',[0.3010,0.7450,0.9330],'LineWidth',3)
plot(t,U_chq_pmq_bpd_3,'b','LineWidth',3)
plot(t,U_dp_pmq_bpd_3,'Color',[0.9290,0.6940,0.1250],'LineWidth',3)
text(250,0.3,['AIC = ',num2str(round(2*nllh_min3+2*10))])
hold off

subplot(4,2,6)
plot(t,U_as_3,'r','LineWidth',3)
ylim([0 1.05])
title('Model 3: 2nd rec.')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_chq_3,'g','LineWidth',3)
plot(t,U_chq_pmq_vhx_3,'Color',[0.3010,0.7450,0.9330],'LineWidth',3)
plot(t,U_chq_pmq_bpd_3,'b','LineWidth',3)
plot(t,U_dp_pmq_bpd_3,'Color',[0.9290,0.6940,0.1250],'LineWidth',3)
text(250,0.3,['AIC = ',num2str(round(2*nllh_min3+2*10))])
hold off

subplot(4,2,7)
plot(t,U_as_4,'r','LineWidth',3)
ylim([0 1.05])
title('Model 4: 1st rec.')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_chq_4,'g','LineWidth',3)
plot(t,U_chq_pmq_vhx_4,'Color',[0.3010,0.7450,0.9330],'LineWidth',3)
plot(t,U_chq_pmq_bpd_4,'b','LineWidth',3)
plot(t,U_dp_pmq_bpd_4,'Color',[0.9290,0.6940,0.1250],'LineWidth',3)
text(250,0.3,['AIC = ',num2str(round(2*nllh_min4+2*11))])
hold off

subplot(4,2,8)
plot(t,U_as_4,'r','LineWidth',3)
ylim([0 1.05])
title('Model 4: 2nd rec.')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_chq_4,'g','LineWidth',3)
plot(t,U_chq_pmq_vhx_4,'Color',[0.3010,0.7450,0.9330],'LineWidth',3)
plot(t,U_chq_pmq_bpd_4,'b','LineWidth',3)
plot(t,U_dp_pmq_bpd_4,'Color',[0.9290,0.6940,0.1250],'LineWidth',3)
text(250,0.3,['AIC = ',num2str(round(2*nllh_min4+2*11))])
hold off


% Comparison of model fits with 1 or 2 reinfection rates (Fig. S17)
% f=1; n=10;

% load survival curve data:
surv_data = cell2mat(struct2cell(load('surv_curve_data_by_drug_and_study.mat')));

% parameter values and negative loglikelihoods
% for 1 reinfection rate in the first row and for 2 reinfection rates in
% the second row:
par1 = [3.4391297,1.61612896,4.04232847,1.14716182,3.89350259,0.05831517,0.0809168,0.00073168,nan;...
    3.44260377,1.61060976,4.051098296,1.15040291,3.88587443,0.05723774,0.08141406,0.001056392,0.00064541822];
nllh1 = [4.3913e+03;4.3881e+03];
par2 = [2.7814,0.23586,3.3785,0.31353,3.8893,0.057723,0.00066443,0.06531,0,0.022986,nan;...
    2.7914,0.23686,3.3847,0.3071,3.882,0.056692,0.00097171,0.00056813,0.070277,0,0.024911];
nllh2 = [4.2405e+03;4.2362e+03];
par3 = [2.9447,0.23984,3.5277,0.329,3.8791,0.059409,0.00063404,-3.7175,2.7055,0,0,10,nan;...
    2.95051812,0.24074922,3.53368025,0.32837685,3.88237672,0.05992895,0.00082751,0.00057377,-3.72783179,2.80439599,0,0,10];
nllh3 = [4.2048e+03;4.2030e+03];

subplot(2,2,1)
stairs([0;surv_data.time(surv_data.drug_study=="AS_VHX")],[1;surv_data.surv(surv_data.drug_study=="AS_VHX")],'r','LineWidth',2)
ylim([0 1.05])
title('Model 1: constant relapse rate')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
stairs([0;surv_data.time(surv_data.drug_study=="CHQ_VHX")],[1;surv_data.surv(surv_data.drug_study=="CHQ_VHX")],'g','LineWidth',2)
stairs([0;surv_data.time(surv_data.drug_study=="CHQ/PMQ_VHX")],[1;surv_data.surv(surv_data.drug_study=="CHQ/PMQ_VHX")],'Color',[0.3010,0.7450,0.9330],'LineWidth',2)
stairs([0;surv_data.time(surv_data.drug_study=="CHQ/PMQ_BPD")],[1;surv_data.surv(surv_data.drug_study=="CHQ/PMQ_BPD")],'b','LineWidth',2)
stairs([0;surv_data.time(surv_data.drug_study=="DP/PMQ_BPD")],[1;surv_data.surv(surv_data.drug_study=="DP/PMQ_BPD")],'Color',[0.9290,0.6940,0.1250],'LineWidth',2)
t = (0:400)';
U_as1 = modelfit(t,1,par1(1,[1,2,7]));
U_chq1 = modelfit(t,1,par1(1,[3,4,7]));
U_cpmq_vhx1 = modelfit(t,1,par1(1,[3,4,8]));
U_cpmq_bpd1 = modelfit(t,1,par1(1,[3,4,8]));
U_dpmq1 = modelfit(t,1,par1(1,[5,6,8]));
plot(t,U_as1,'Color',[1,0,0],'LineWidth',2)
plot(t,U_chq1,'Color',[0,1,0],'LineWidth',2)
plot(t,U_cpmq_vhx1,'Color',[0.3010,0.7450,0.9330],'LineWidth',2)
plot(t,U_cpmq_bpd1,'Color',[0,0,1],'LineWidth',2)
plot(t,U_dpmq1,'Color',[0.9290,0.6940,0.1250],'LineWidth',2)
U_as2 = modelfit(t,1,par1(2,[1,2,7]));
U_chq2 = modelfit(t,1,par1(2,[3,4,7]));
U_cpmq_vhx2 = modelfit(t,1,par1(2,[3,4,8]));
U_cpmq_bpd2 = modelfit(t,1,par1(2,[3,4,9]));
U_dpmq2 = modelfit(t,1,par1(2,[5,6,9]));
plot(t,U_as2,'--','Color',0.8*[1,0,0],'LineWidth',2)
plot(t,U_chq2,'--','Color',0.8*[0,1,0],'LineWidth',2)
plot(t,U_cpmq_vhx2,'--','Color',0.8*[0.3010,0.7450,0.9330],'LineWidth',2)
plot(t,U_cpmq_bpd2,'--','Color',0.8*[0,0,1],'LineWidth',2)
plot(t,U_dpmq2,'--','Color',0.8*[0.9290,0.6940,0.1250],'LineWidth',2)
legend('AS (VHX)','CHQ (VHX)','CHQ/PMQ (VHX)','CHQ/PMQ (BPD)','DP/PMQ (BPD)','Location','east')
hold off

subplot(2,2,2)
stairs([0;surv_data.time(surv_data.drug_study=="AS_VHX")],[1;surv_data.surv(surv_data.drug_study=="AS_VHX")],'r','LineWidth',2)
ylim([0 1.05])
title('Model 2: temporal heterogeneity')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
stairs([0;surv_data.time(surv_data.drug_study=="CHQ_VHX")],[1;surv_data.surv(surv_data.drug_study=="CHQ_VHX")],'g','LineWidth',2)
stairs([0;surv_data.time(surv_data.drug_study=="CHQ/PMQ_VHX")],[1;surv_data.surv(surv_data.drug_study=="CHQ/PMQ_VHX")],'Color',[0.3010,0.7450,0.9330],'LineWidth',2)
stairs([0;surv_data.time(surv_data.drug_study=="CHQ/PMQ_BPD")],[1;surv_data.surv(surv_data.drug_study=="CHQ/PMQ_BPD")],'b','LineWidth',2)
stairs([0;surv_data.time(surv_data.drug_study=="DP/PMQ_BPD")],[1;surv_data.surv(surv_data.drug_study=="DP/PMQ_BPD")],'Color',[0.9290,0.6940,0.1250],'LineWidth',2)
t = (0:400)';
U_as1 = modelfit(t,2,par2(1,[1,2,7,8,10]));
U_chq1 = modelfit(t,2,par2(1,[3,4,7,8,10]));
U_cpmq_vhx1 = modelfit(t,2,par2(1,[3,4,7,9,10]));
U_cpmq_bpd1 = modelfit(t,2,par2(1,[3,4,7,9,10]));
U_dpmq1 = modelfit(t,2,par2(1,[5,6,7,9,10]));
plot(t,U_as1,'Color',[1,0,0],'LineWidth',2)
plot(t,U_chq1,'Color',[0,1,0],'LineWidth',2)
plot(t,U_cpmq_vhx1,'Color',[0.3010,0.7450,0.9330],'LineWidth',2)
plot(t,U_cpmq_bpd1,'Color',[0,0,1],'LineWidth',2)
plot(t,U_dpmq1,'Color',[0.9290,0.6940,0.1250],'LineWidth',2)
U_as2 = modelfit(t,2,par2(2,[1,2,7,9,11]));
U_chq2 = modelfit(t,2,par2(2,[3,4,7,9,11]));
U_cpmq_vhx2 = modelfit(t,2,par2(2,[3,4,7,10,11]));
U_cpmq_bpd2 = modelfit(t,2,par2(2,[3,4,8,10,11]));
U_dpmq2 = modelfit(t,2,par2(2,[5,6,8,10,11]));
plot(t,U_as2,'--','Color',0.8*[1,0,0],'LineWidth',2)
plot(t,U_chq2,'--','Color',0.8*[0,1,0],'LineWidth',2)
plot(t,U_cpmq_vhx2,'--','Color',0.8*[0.3010,0.7450,0.9330],'LineWidth',2)
plot(t,U_cpmq_bpd2,'--','Color',0.8*[0,0,1],'LineWidth',2)
plot(t,U_dpmq2,'--','Color',0.8*[0.9290,0.6940,0.1250],'LineWidth',2)
legend('AS (VHX)','CHQ (VHX)','CHQ/PMQ (VHX)','CHQ/PMQ (BPD)','DP/PMQ (BPD)','Location','east')
hold off

subplot(2,2,3)
stairs([0;surv_data.time(surv_data.drug_study=="AS_VHX")],[1;surv_data.surv(surv_data.drug_study=="AS_VHX")],'r','LineWidth',2)
ylim([0 1.05])
title('Model 3: population heterogeneity')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
stairs([0;surv_data.time(surv_data.drug_study=="CHQ_VHX")],[1;surv_data.surv(surv_data.drug_study=="CHQ_VHX")],'g','LineWidth',2)
stairs([0;surv_data.time(surv_data.drug_study=="CHQ/PMQ_VHX")],[1;surv_data.surv(surv_data.drug_study=="CHQ/PMQ_VHX")],'Color',[0.3010,0.7450,0.9330],'LineWidth',2)
stairs([0;surv_data.time(surv_data.drug_study=="CHQ/PMQ_BPD")],[1;surv_data.surv(surv_data.drug_study=="CHQ/PMQ_BPD")],'b','LineWidth',2)
stairs([0;surv_data.time(surv_data.drug_study=="DP/PMQ_BPD")],[1;surv_data.surv(surv_data.drug_study=="DP/PMQ_BPD")],'Color',[0.9290,0.6940,0.1250],'LineWidth',2)
t = (0:400)';
U_as1 = sum(modelfit(t,3,par3(1,[1,2,7,12,8,9])),1);
U_chq1 = sum(modelfit(t,3,par3(1,[3,4,7,12,8,9])),1);
U_cpmq_vhx1 = sum(modelfit(t,3,par3(1,[3,4,7,12,10,11])),1);
U_cpmq_bpd1 = sum(modelfit(t,3,par3(1,[3,4,7,12,10,11])),1);
U_dpmq1 = sum(modelfit(t,3,par3(1,[5,6,7,12,10,11])),1);
plot(t,U_as1,'Color',[1,0,0],'LineWidth',2)
plot(t,U_chq1,'Color',[0,1,0],'LineWidth',2)
plot(t,U_cpmq_vhx1,'Color',[0.3010,0.7450,0.9330],'LineWidth',2)
plot(t,U_cpmq_bpd1,'Color',[0,0,1],'LineWidth',2)
plot(t,U_dpmq1,'Color',[0.9290,0.6940,0.1250],'LineWidth',2)
U_as2 = sum(modelfit(t,3,par3(2,[1,2,7,13,9,10])),1);
U_chq2 = sum(modelfit(t,3,par3(2,[3,4,7,13,9,10])),1);
U_cpmq_vhx2 = sum(modelfit(t,3,par3(2,[3,4,7,13,11,12])),1);
U_cpmq_bpd2 = sum(modelfit(t,3,par3(2,[3,4,8,13,11,12])),1);
U_dpmq2 = sum(modelfit(t,3,par3(2,[5,6,8,13,11,12])),1);
plot(t,U_as2,'--','Color',0.8*[1,0,0],'LineWidth',2)
plot(t,U_chq2,'--','Color',0.8*[0,1,0],'LineWidth',2)
plot(t,U_cpmq_vhx2,'--','Color',0.8*[0.3010,0.7450,0.9330],'LineWidth',2)
plot(t,U_cpmq_bpd2,'--','Color',0.8*[0,0,1],'LineWidth',2)
plot(t,U_dpmq2,'--','Color',0.8*[0.9290,0.6940,0.1250],'LineWidth',2)
legend('AS (VHX)','CHQ (VHX)','CHQ/PMQ (VHX)','CHQ/PMQ (BPD)','DP/PMQ (BPD)','AS (VHX)','CHQ (VHX)','CHQ/PMQ (VHX)','CHQ/PMQ (BPD)','DP/PMQ (BPD)',...
    'AS (VHX)','CHQ (VHX)','CHQ/PMQ (VHX)','CHQ/PMQ (BPD)','DP/PMQ (BPD)','Location','east')
hold off

% AICs: (Table S23)
fun_aic = @(n_par,nllh) 2.*n_par+2.*nllh; % for computing the AIC
n_par = [9,10,10];

% model 1:
aic1_1 = fun_aic(n_par(1)-1,nllh1(1)); % 1 reinfection rate
aic1_2 = fun_aic(n_par(1),nllh1(2)); % 2 reinfection rates
% model 2:
aic2_1 = fun_aic(n_par(2)-1,nllh2(1)); % 1 reinfection rate
aic2_2 = fun_aic(n_par(2),nllh2(2)); % 2 reinfection rates
% model 3:
aic3_1 = fun_aic(n_par(3)-1,nllh3(1)); % 1 reinfection rate
aic3_2 = fun_aic(n_par(3),nllh3(2)); % 2 reinfection rates



% Comparison of model fits with different follow-up schemes (Fig. S18-20)
% reinf = 2; n=10;

% load survival curve data:
surv_data = cell2mat(struct2cell(load('surv_curve_data_by_drug_and_study.mat')));

% parameter values and negative loglikelihoods
% for each row contains the parameter values for a different follow-up scheme
par1 = [3.44260377,1.61060976,4.051098296,1.15040291,3.88587443,0.05723774,0.08141406,0.001056392,0.00064541822;...
    3.25528984,1.83093938,3.977038737,1.22871697,3.82570173,0.04870433,0.082020754,0.001056794,0.000640362;...
    2.88256452,2.30270699,3.87114951,1.34928499,3.7143667,0.03271831,0.08110879,0.001046854,0.000630766;...
    1.54875928,4.10050329,3.64224885,1.64977028,3.52952104,0.02078084,0.0942147,0.001036436,0.000619133;...
    2.7744,2.3836,3.7592,1.4554,3.5353,0.070517,0.083922,0.0010353,0.00061911;...
    2.94936003,2.16422512,3.810031807,1.392310901,3.540521788,0.02519966,0.082028755,0.001038405,0.000619971;...
    3.07594585,2.01133962,3.853021409,1.34118233,3.53641982,0.069025375,0.07991755,0.001041173,0.000620325];
nllh1 = [4.3881e+03,2.9236e+03,2.3843e+03,1.8472e+03,2.1356e+03,2.2920e+03,2.5024e+03];
par2 = [2.7914,0.23686,3.3847,0.3071,3.882,0.056692,0.00097171,0.00056813,0.070277,0,0.024911;
    2.5860,0.2290,3.2696,0.3272,3.8220,0.0482,0.000956467,0.00056354,0.0645,0,0.0249;
    2.2770,0.2223,3.1155,0.3194,3.704856984,0.03132212197,0.000946449,0.000553096,0.0584,0,0.0250;
    -28.8705,22.6540,2.7168,0.2372,3.4457,0.0577,0.000888602,0.000536003,0.0438,0,0.0213;
    2.1242,0.0959,2.9930,0.2305,3.4428,0.1889,0.000973178,0.000541365,0.0639,0,0.0285;
    2.4012,0.1367,3.0927,0.2538,3.4861,0.1482,0.001015644,0.000545369,0.0751,0,0.0304;
    2.5713,0.1688,3.1755,0.2739,3.4600,0.0482,0.001025447,0.000545646,0.0829,0,0.0310];
nllh2 = [4.2362e+03,2.7949e+03,2.2879e+03,1.8164e+03,2.0100e+03,2.1493e+03,2.3505e+03];
par3 = [2.95051812,0.24074922,3.53368025,0.32837685,3.88237672,0.05992895,0.00082751,0.00057377,-3.72783179,2.80439599,0,0,10;
    2.7681,0.2368,3.4288,0.3417,3.8176,0.0503,0.000811145,0.000568073,-3.7833,2.7445,0,0,10;
    2.5477,0.2227,3.3021,0.3394,3.7045,0.0331,0.0008028323,0.000559069,-3.8062,2.7774,0,0,10;
    2.009593696,0.16045247,2.993471699,0.2963044078,3.445498662,0.02797017,0.000785771,0.000540561,-3.835265919,2.856245202,0,0,10;
    2.2650,0.1284,3.1318,0.2226,3.4564,0.0633,0.000778126,0.000544682,-3.8855,2.6541,0,0,10;
    2.5167,0.1360,3.2231,0.2511,-20.9493,39.7260,0.000792029,0.000597037,-3.8185,2.7628,0,0,10;
    2.6949,0.1585,3.3016,0.2753,-10.2251,24.5367,0.0008042049,0.0006039908,-3.7529,2.8629,0,0,10];
nllh3 = [4.2030e+03,2.7620e+03,2.2562e+03,1.7929e+03,1.9850e+03,2.1219e+03,2.3217e+03];

colors = {'k','b','g','r','c','m','y'};

m = 3; % model to be plotted

subplot(3,2,1)
stairs([0;surv_data.time(surv_data.drug_study=="AS_VHX")],[1;surv_data.surv(surv_data.drug_study=="AS_VHX")],'Color',[200 200 200]/255,'LineWidth',3)
ylim([0 1.05])
title('AS (VHX)')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
t = (0:400)';
for f=1:7
    if m==1
        U_as = modelfit(t,m,par1(f,[1,2,7]));
    elseif m==2
        U_as = modelfit(t,m,par2(f,[1,2,7,9,11]));
    elseif m==3
        U_as = sum(modelfit(t,m,par3(f,[1,2,7,13,9,10])),1);
    end
    plot(t,U_as,'-','LineWidth',2,'Color',colors{f})
end
% legend('data','daily follow-up','weekly','fortnightly','4-weekly',...
%     'beginning of weeks 2,4,8,12,...','middle of weeks 2,4,8,12,...',...
%     'end of weeks 2,4,8,12,...')
hold off

subplot(3,2,2)
stairs([0;surv_data.time(surv_data.drug_study=="CHQ_VHX")],[1;surv_data.surv(surv_data.drug_study=="CHQ_VHX")],'Color',[200 200 200]/255,'LineWidth',2)
ylim([0 1.05])
title('CHQ (VHX)')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
t = (0:400)';
for f=1:7
    if m==1
        U_chq = modelfit(t,m,par1(f,[3,4,7]));
    elseif m==2
        U_chq = modelfit(t,m,par2(f,[3,4,7,9,11]));
    elseif m==3
        U_chq = sum(modelfit(t,m,par3(f,[3,4,7,13,9,10])),1);
    end
    plot(t,U_chq,'-','LineWidth',2,'Color',colors{f})
end
% legend('data','daily follow-up','weekly','fortnightly','4-weekly',...
%     'beginning of weeks 2,4,8,12,...','middle of weeks 2,4,8,12,...',...
%     'end of weeks 2,4,8,12,...')
hold off

subplot(3,2,3)
stairs([0;surv_data.time(surv_data.drug_study=="CHQ/PMQ_VHX")],[1;surv_data.surv(surv_data.drug_study=="CHQ/PMQ_VHX")],'Color',[200 200 200]/255,'LineWidth',2)
ylim([0 1.05])
title('CHQ/PMQ (VHX)')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
t = (0:400)';
for f=1:7
    if m==1
        U_cpmq = modelfit(t,m,par1(f,[3,4,8]));
    elseif m==2
        U_cpmq = modelfit(t,m,par2(f,[3,4,7,10,11]));
    elseif m==3
        U_cpmq = sum(modelfit(t,m,par3(f,[3,4,7,13,11,12])),1);
    end
    plot(t,U_cpmq,'-','LineWidth',2,'Color',colors{f})
end
% legend('data','daily follow-up','weekly','fortnightly','4-weekly',...
%     'beginning of weeks 2,4,8,12,...','middle of weeks 2,4,8,12,...',...
%     'end of weeks 2,4,8,12,...')
hold off

subplot(3,2,4)
stairs([0;surv_data.time(surv_data.drug_study=="CHQ/PMQ_BPD")],[1;surv_data.surv(surv_data.drug_study=="CHQ/PMQ_BPD")],'Color',[200 200 200]/255,'LineWidth',2)
ylim([0 1.05])
title('CHQ/PMQ (BPD)')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
t = (0:400)';
for f=1:7
    if m==1
        U_cpmq = modelfit(t,m,par1(f,[3,4,9]));
    elseif m==2
        U_cpmq = modelfit(t,m,par2(f,[3,4,8,10,11]));
    elseif m==3
        U_cpmq = sum(modelfit(t,m,par3(f,[3,4,8,13,11,12])),1);
    end
    plot(t,U_cpmq,'-','LineWidth',2,'Color',colors{f})
end
% legend('data','daily follow-up','weekly','fortnightly','4-weekly',...
%     'beginning of weeks 2,4,8,12,...','middle of weeks 2,4,8,12,...',...
%     'end of weeks 2,4,8,12,...')
hold off

subplot(3,2,5)
stairs([0;surv_data.time(surv_data.drug_study=="DP/PMQ_BPD")],[1;surv_data.surv(surv_data.drug_study=="DP/PMQ_BPD")],'Color',[200 200 200]/255,'LineWidth',2)
ylim([0 1.05])
title('DP/PMQ (BPD)')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
t = (0:400)';
for f=1:7
    if m==1
        U_dpmq = modelfit(t,m,par1(f,[5,6,9]));
    elseif m==2
        U_dpmq = modelfit(t,m,par2(f,[5,6,8,10,11]));
    elseif m==3
        U_dpmq = sum(modelfit(t,m,par3(f,[5,6,8,13,11,12])),1);
    end
    plot(t,U_dpmq,'-','LineWidth',2,'Color',colors{f})
end
legend('data','daily follow-up','weekly','fortnightly','4-weekly',...
    'beginning of weeks 2,4,8,12,...','middle of weeks 2,4,8,12,...',...
    'end of weeks 2,4,8,12,...')
hold off

% AICs: (Table S24)
fun_aic = @(n_par,nllh) 2.*n_par+2.*nllh; % for computing the AIC
n_par = [9,10,10];

% model 1:
aic1_1 = fun_aic(n_par(1),nllh1(1)); % follow-up scheme 1
aic1_2 = fun_aic(n_par(1),nllh1(2)); % follow-up scheme 2
aic1_3 = fun_aic(n_par(1),nllh1(3)); % follow-up scheme 3
aic1_4 = fun_aic(n_par(1),nllh1(4)); % follow-up scheme 4
aic1_5 = fun_aic(n_par(1),nllh1(5)); % follow-up scheme 5
aic1_6 = fun_aic(n_par(1),nllh1(6)); % follow-up scheme 6
aic1_7 = fun_aic(n_par(1),nllh1(7)); % follow-up scheme 7
% model 2:
aic2_1 = fun_aic(n_par(2),nllh2(1)); % follow-up scheme 1
aic2_2 = fun_aic(n_par(2),nllh2(2)); % follow-up scheme 2
aic2_3 = fun_aic(n_par(2),nllh2(3)); % follow-up scheme 3
aic2_4 = fun_aic(n_par(2),nllh2(4)); % follow-up scheme 4
aic2_5 = fun_aic(n_par(2),nllh2(5)); % follow-up scheme 5
aic2_6 = fun_aic(n_par(2),nllh2(6)); % follow-up scheme 6
aic2_7 = fun_aic(n_par(2),nllh2(7)); % follow-up scheme 7
% model 3:
aic3_1 = fun_aic(n_par(3),nllh3(1)); % follow-up scheme 1
aic3_2 = fun_aic(n_par(3),nllh3(2)); % follow-up scheme 2
aic3_3 = fun_aic(n_par(3),nllh3(3)); % follow-up scheme 3
aic3_4 = fun_aic(n_par(3),nllh3(4)); % follow-up scheme 4
aic3_5 = fun_aic(n_par(3),nllh3(5)); % follow-up scheme 5
aic3_6 = fun_aic(n_par(3),nllh3(6)); % follow-up scheme 6
aic3_7 = fun_aic(n_par(3),nllh3(7)); % follow-up scheme 7



% Comparison of model 3 fit with different number of relapse risk groups
% reinf = 2; f=1; (Fig. S21)

par = table('Size',[16,3],'VariableTypes',{'double','double','double'},...
    'VariableNames',{'model1','model2','model3'},...
    'rowNames',{'washout1_AS','washout2_AS','washout1_CHQ','washout2_CHQ','washout1_DP',...
    'washout2_DP','recurrence_AS/CHQ_VHX','recurrence_PMQ_VHX','recurrence_PMQ_BPD',...
    'reinfection_VHX','reinfection_BPD','relapse_risk1_AS/CHQ','relapse_risk2_AS/CHQ',...
    'number_relapse_groups','I_AS/CHQ','d'});
n_par = [9,10,10]; % number of estimated parameters for each model

subplot(3,2,1)
ylim([0 1.05])
title('AS (VHX)')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
t = (0:400)';
% n=5;
par{[1:6,10:14],3} = [2.9639;0.2429;3.5502;0.3361;3.8775;0.0673;0.000838419;0.000575624;-3.6376;3.0019;5];
U_as1 = sum(modelfit(t,3,par.model3([1,2,10,14,12,13])),1);
plot(t,U_as1,'k-','LineWidth',2)
% n=10;
par{[1:6,10:14],3} = [2.95051812;0.24074922;3.53368025;0.32837685;3.88237672;0.05992895;0.00082751;0.00057377;-3.72783179;2.80439599;10];
U_as2 = sum(modelfit(t,3,par.model3([1,2,10,14,12,13])),1);
plot(t,U_as2,'--','Color',[0,1,0],'LineWidth',2)
% n=15;
par{[1:6,10:14],3} = [2.9484;0.2388;3.5306;0.3267;3.8794;0.0607;0.000822403;0.000574456;-3.7413;2.7767;15];
U_as3 = sum(modelfit(t,3,par.model3([1,2,10,14,12,13])),1);
plot(t,U_as3,':','Color',[1,0,0],'LineWidth',2)
% n=20;
par{[1:6,10:14],3} = [2.9482;0.2390;3.5315;0.3271;3.8830;0.0611;0.00082275;0.00057491;-3.7423;2.7718;20];
U_as4 = sum(modelfit(t,3,par.model3([1,2,10,14,12,13])),1);
plot(t,U_as4,'-.','Color',[0,0,1],'LineWidth',2)
legend('k=5','k=10','k=15','k=20','Location','northeast')
hold off

subplot(3,2,2)
ylim([0 1.05])
title('CHQ (VHX)')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
t = (0:400)';
% n=5;
par{[1:6,10:14],3} = [2.9639;0.2429;3.5502;0.3361;3.8775;0.0673;0.000838419;0.000575624;-3.6376;3.0019;5];
U_chq1 = sum(modelfit(t,3,par.model3([3,4,10,14,12,13])),1);
plot(t,U_chq1,'k-','LineWidth',2)
% n=10;
par{[1:6,10:14],3} = [2.95051812;0.24074922;3.53368025;0.32837685;3.88237672;0.05992895;0.00082751;0.00057377;-3.72783179;2.80439599;10];
U_chq2 = sum(modelfit(t,3,par.model3([3,4,10,14,12,13])),1);
plot(t,U_chq2,'--','Color',[0,1,0],'LineWidth',2)
% n=15;
par{[1:6,10:14],3} = [2.9484;0.2388;3.5306;0.3267;3.8794;0.0607;0.000822403;0.000574456;-3.7413;2.7767;15];
U_chq3 = sum(modelfit(t,3,par.model3([3,4,10,14,12,13])),1);
plot(t,U_chq3,':','Color',[1,0,0],'LineWidth',2)
% n=20;
par{[1:6,10:14],3} = [2.9482;0.2390;3.5315;0.3271;3.8830;0.0611;0.00082275;0.00057491;-3.7423;2.7718;20];
U_chq4 = sum(modelfit(t,3,par.model3([3,4,10,14,12,13])),1);
plot(t,U_chq4,'-.','Color',[0,0,1],'LineWidth',2)
legend('k=5','k=10','k=15','k=20','Location','northeast')
hold off

subplot(3,2,3)
ylim([0 1.05])
title('CHQ/PMQ (VHX)')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
t = (0:400)';
% n=5;
par{[1:6,10:14],3} = [2.9639;0.2429;3.5502;0.3361;3.8775;0.0673;0.000838419;0.000575624;-3.6376;3.0019;5];
U_cpmq_vhx1 = sum(modelfit(t,3,[par.model3([3,4,10,14])',0,0]),1);
plot(t,U_cpmq_vhx1,'k-','LineWidth',2)
% n=10;
par{[1:6,10:14],3} = [2.95051812;0.24074922;3.53368025;0.32837685;3.88237672;0.05992895;0.00082751;0.00057377;-3.72783179;2.80439599;10];
U_cpmq_vhx2 = sum(modelfit(t,3,[par.model3([3,4,10,14])',0,0]),1);
plot(t,U_cpmq_vhx2,'--','Color',[0,1,0],'LineWidth',2)
% n=15;
par{[1:6,10:14],3} = [2.9484;0.2388;3.5306;0.3267;3.8794;0.0607;0.000822403;0.000574456;-3.7413;2.7767;15];
U_cpmq_vhx3 = sum(modelfit(t,3,[par.model3([3,4,10,14])',0,0]),1);
plot(t,U_cpmq_vhx3,':','Color',[1,0,0],'LineWidth',2)
% n=20;
par{[1:6,10:14],3} = [2.9482;0.2390;3.5315;0.3271;3.8830;0.0611;0.00082275;0.00057491;-3.7423;2.7718;20];
U_cpmq_vhx4 = sum(modelfit(t,3,[par.model3([3,4,10,14])',0,0]),1);
plot(t,U_cpmq_vhx4,'-.','Color',[0,0,1],'LineWidth',2)
legend('k=5','k=10','k=15','k=20','Location','northeast')
hold off

subplot(3,2,4)
ylim([0 1.05])
title('CHQ/PMQ (BPD)')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
t = (0:400)';
% n=5;
par{[1:6,10:14],3} = [2.9639;0.2429;3.5502;0.3361;3.8775;0.0673;0.000838419;0.000575624;-3.6376;3.0019;5];
U_cpmq_bpd1 = sum(modelfit(t,3,[par.model3([3,4,11,14])',0,0]),1);
plot(t,U_cpmq_bpd1,'k-','LineWidth',2)
% n=10;
par{[1:6,10:14],3} = [2.95051812;0.24074922;3.53368025;0.32837685;3.88237672;0.05992895;0.00082751;0.00057377;-3.72783179;2.80439599;10];
U_cpmq_bpd2 = sum(modelfit(t,3,[par.model3([3,4,11,14])',0,0]),1);
plot(t,U_cpmq_bpd2,'--','Color',[0,1,0],'LineWidth',2)
% n=15;
par{[1:6,10:14],3} = [2.9484;0.2388;3.5306;0.3267;3.8794;0.0607;0.000822403;0.000574456;-3.7413;2.7767;15];
U_cpmq_bpd3 = sum(modelfit(t,3,[par.model3([3,4,11,14])',0,0]),1);
plot(t,U_cpmq_bpd3,':','Color',[1,0,0],'LineWidth',2)
% n=20;
par{[1:6,10:14],3} = [2.9482;0.2390;3.5315;0.3271;3.8830;0.0611;0.00082275;0.00057491;-3.7423;2.7718;20];
U_cpmq_bpd4 = sum(modelfit(t,3,[par.model3([3,4,11,14])',0,0]),1);
plot(t,U_cpmq_bpd4,'-.','Color',[0,0,1],'LineWidth',2)
legend('k=5','k=10','k=15','k=20','Location','northeast')
hold off

subplot(3,2,5)
ylim([0 1.05])
title('DP/PMQ (BPD)')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
t = (0:400)';
% n=5;
par{[1:6,10:14],3} = [2.9639;0.2429;3.5502;0.3361;3.8775;0.0673;0.000838419;0.000575624;-3.6376;3.0019;5];
U_dpmq_bpd1 = sum(modelfit(t,3,[par.model3([5,6,11,14])',0,0]),1);
plot(t,U_dpmq_bpd1,'k-','LineWidth',2)
% n=10;
par{[1:6,10:14],3} = [2.95051812;0.24074922;3.53368025;0.32837685;3.88237672;0.05992895;0.00082751;0.00057377;-3.72783179;2.80439599;10];
U_dpmq_bpd2 = sum(modelfit(t,3,[par.model3([5,6,11,14])',0,0]),1);
plot(t,U_dpmq_bpd2,'--','Color',[0,1,0],'LineWidth',2)
% n=15;
par{[1:6,10:14],3} = [2.9484;0.2388;3.5306;0.3267;3.8794;0.0607;0.000822403;0.000574456;-3.7413;2.7767;15];
U_dpmq_bpd3 = sum(modelfit(t,3,[par.model3([5,6,11,14])',0,0]),1);
plot(t,U_dpmq_bpd3,':','Color',[1,0,0],'LineWidth',2)
% n=20;
par{[1:6,10:14],3} = [2.9482;0.2390;3.5315;0.3271;3.8830;0.0611;0.00082275;0.00057491;-3.7423;2.7718;20];
U_dpmq_bpd4 = sum(modelfit(t,3,[par.model3([5,6,11,14])',0,0]),1);
plot(t,U_dpmq_bpd4,'-.','Color',[0,0,1],'LineWidth',2)
legend('k=5','k=10','k=15','k=20','Location','northeast')
hold off

% AICs: (Table S25)
fun_aic = @(n_par,nllh) 2.*n_par+2.*nllh; % for computing the AIC

% n=5;
nllh5 = 4.2028e+03;
aic5 = fun_aic(n_par(3),nllh5);
% n=10;
nllh10 = 4.2030e+03;
aic10 = fun_aic(n_par(3),nllh10);
% n=15;
nllh15 = 4.2030e+03;
aic15 = fun_aic(n_par(3),nllh15);
% n=20;
nllh20 = 4.2030e+03;
aic20 = fun_aic(n_par(3),nllh20);


% Comparison of the population heterogeneity model with different relapse
% risk distributions (Fig. S22)

% Parameter values:
nllh_min3 = 4.2030e+03;
par_min3 = [2.95051812;0.24074922;3.53368025;0.32837685;3.88237672;0.05992895;0.00082751;0.00057377;-3.72783179;2.80439599;10];
nllh_min3_gam = 4.2041e+03;
par_min3_gam = [2.903686703943764,0.241085890062690,3.513597074580805,0.338705364820449,...
    3.875813718206521,0.058915640688840,0.000903125587265,0.000578465253674,0.303098845785234,0.384639708871805,10];
nllh_min3_exp = 4.2567e+03;
par_min3_exp = [2.718041005749698,0.220700541966448,3.338564771709158,0.311841142475142,...
    3.886725287979828,0.060571307252208,0.000642750709578,0.000567349056870,0.035443379089622,10];

t = 0:400;
U_as_3 = sum(modelfit(t,3,[par_min3(1),par_min3(2),par_min3(7),par_min3(11),par_min3(9),par_min3(10)]),1);
U_chq_3 = sum(modelfit(t,3,[par_min3(3),par_min3(4),par_min3(7),par_min3(11),par_min3(9),par_min3(10)]),1);
U_chq_pmq_vhx_3 = sum(modelfit(t,3,[par_min3(3),par_min3(4),par_min3(7),par_min3(11),0,0]),1);
U_chq_pmq_bpd_3 = sum(modelfit(t,3,[par_min3(3),par_min3(4),par_min3(8),par_min3(11),0,0]),1);
U_dp_pmq_bpd_3 = sum(modelfit(t,3,[par_min3(5),par_min3(6),par_min3(8),par_min3(11),0,0]),1);
U_as_3_gam = sum(modelfit(t,3,[par_min3_gam(1),par_min3_gam(2),par_min3_gam(7),par_min3_gam(11),par_min3_gam(9),par_min3_gam(10)],"gamma"),1);
U_chq_3_gam = sum(modelfit(t,3,[par_min3_gam(3),par_min3_gam(4),par_min3_gam(7),par_min3_gam(11),par_min3_gam(9),par_min3_gam(10)],"gamma"),1);
U_chq_pmq_vhx_3_gam = sum(modelfit(t,3,[par_min3_gam(3),par_min3_gam(4),par_min3_gam(7),par_min3_gam(11),0,0],"gamma"),1);
U_chq_pmq_bpd_3_gam = sum(modelfit(t,3,[par_min3_gam(3),par_min3_gam(4),par_min3_gam(8),par_min3_gam(11),0,0],"gamma"),1);
U_dp_pmq_bpd_3_gam = sum(modelfit(t,3,[par_min3_gam(5),par_min3_gam(6),par_min3_gam(8),par_min3_gam(11),0,0],"gamma"),1);
U_as_3_exp = sum(modelfit(t,3,[par_min3_exp(1),par_min3_exp(2),par_min3_exp(7),par_min3_exp(10),par_min3_exp(9),0],"exp"),1);
U_chq_3_exp = sum(modelfit(t,3,[par_min3_exp(3),par_min3_exp(4),par_min3_exp(7),par_min3_exp(10),par_min3_exp(9),0],"exp"),1);
U_chq_pmq_vhx_3_exp = sum(modelfit(t,3,[par_min3_exp(3),par_min3_exp(4),par_min3_exp(7),par_min3_exp(10),0,0],"exp"),1);
U_chq_pmq_bpd_3_exp = sum(modelfit(t,3,[par_min3_exp(3),par_min3_exp(4),par_min3_exp(8),par_min3_exp(10),0,0],"exp"),1);
U_dp_pmq_bpd_3_exp = sum(modelfit(t,3,[par_min3_exp(5),par_min3_exp(6),par_min3_exp(8),par_min3_exp(10),0,0],"exp"),1);

subplot(3,2,1)
plot(t,U_as_3,'r','LineWidth',3)
ylim([0 1.05])
title('Model 3: 1st rec. (lognormal)')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_chq_3,'g','LineWidth',3)
plot(t,U_chq_pmq_vhx_3,'Color',[0.3010,0.7450,0.9330],'LineWidth',3)
plot(t,U_chq_pmq_bpd_3,'b','LineWidth',3)
plot(t,U_dp_pmq_bpd_3,'Color',[0.9290,0.6940,0.1250],'LineWidth',3)
text(250,0.3,['AIC = ',num2str(round(2*nllh_min3+2*10))])
hold off

subplot(3,2,2)
plot(t,U_as_3,'r','LineWidth',3)
ylim([0 1.05])
title('Model 3: 2nd rec. (lognormal)')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_chq_3,'g','LineWidth',3)
plot(t,U_chq_pmq_vhx_3,'Color',[0.3010,0.7450,0.9330],'LineWidth',3)
plot(t,U_chq_pmq_bpd_3,'b','LineWidth',3)
plot(t,U_dp_pmq_bpd_3,'Color',[0.9290,0.6940,0.1250],'LineWidth',3)
text(250,0.3,['AIC = ',num2str(round(2*nllh_min3+2*10))])
hold off

subplot(3,2,3)
plot(t,U_as_3_gam,'r','LineWidth',3)
ylim([0 1.05])
title('Model 3: 1st rec. (gamma)')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_chq_3_gam,'g','LineWidth',3)
plot(t,U_chq_pmq_vhx_3_gam,'Color',[0.3010,0.7450,0.9330],'LineWidth',3)
plot(t,U_chq_pmq_bpd_3_gam,'b','LineWidth',3)
plot(t,U_dp_pmq_bpd_3_gam,'Color',[0.9290,0.6940,0.1250],'LineWidth',3)
text(250,0.3,['AIC = ',num2str(round(2*nllh_min3_gam+2*10))])
hold off

subplot(3,2,4)
plot(t,U_as_3_gam,'r','LineWidth',3)
ylim([0 1.05])
title('Model 3: 2nd rec. (gamma)')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_chq_3_gam,'g','LineWidth',3)
plot(t,U_chq_pmq_vhx_3_gam,'Color',[0.3010,0.7450,0.9330],'LineWidth',3)
plot(t,U_chq_pmq_bpd_3_gam,'b','LineWidth',3)
plot(t,U_dp_pmq_bpd_3_gam,'Color',[0.9290,0.6940,0.1250],'LineWidth',3)
text(250,0.3,['AIC = ',num2str(round(2*nllh_min3_gam+2*10))])
hold off

subplot(3,2,5)
plot(t,U_as_3_exp,'r','LineWidth',3)
ylim([0 1.05])
title('Model 3: 1st rec. (exponential)')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_chq_3_exp,'g','LineWidth',3)
plot(t,U_chq_pmq_vhx_3_exp,'Color',[0.3010,0.7450,0.9330],'LineWidth',3)
plot(t,U_chq_pmq_bpd_3_exp,'b','LineWidth',3)
plot(t,U_dp_pmq_bpd_3_exp,'Color',[0.9290,0.6940,0.1250],'LineWidth',3)
text(250,0.3,['AIC = ',num2str(round(2*nllh_min3_exp+2*9))])
hold off

subplot(3,2,6)
plot(t,U_as_3_exp,'r','LineWidth',3)
ylim([0 1.05])
title('Model 3: 2nd rec. (exponential)')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_chq_3_exp,'g','LineWidth',3)
plot(t,U_chq_pmq_vhx_3_exp,'Color',[0.3010,0.7450,0.9330],'LineWidth',3)
plot(t,U_chq_pmq_bpd_3_exp,'b','LineWidth',3)
plot(t,U_dp_pmq_bpd_3_exp,'Color',[0.9290,0.6940,0.1250],'LineWidth',3)
text(250,0.3,['AIC = ',num2str(round(2*nllh_min3_exp+2*9))])
hold off


% Plot the relapse risk distributions for lognormal, gamma and exponential
% distribution (Fig. S23)

par_ln = [-3.72783179,2.80439599]; % parameters for the relapse rate distribution model fit with lognormal distribution
par_gam = [0.303098845785234,0.384639708871805]; % parameters for the relapse rate distribution model fit with lognormal distribution
par_exp = 0.035443379089622; % parameter for the relapse rate distribution model fit with exponential distribution


u=0.01;
du=0.00001;
plot(0:du:u,lognpdf(0:du:u,par_ln(1),par_ln(2)),'k','LineWidth',2)
hold on
plot(0:du:u,gampdf(0:du:u,par_gam(1),par_gam(2)),'r','LineWidth',2)
plot(0:du:u,exppdf(0:du:u,par_exp),'b','LineWidth',2)
ylim([0,500])
xlim([0,u])
xlabel('Relapse rate')
ylabel('Density')
title('Distribution of the relapse rates in the population heterogeneity model')
legend('Lognormal distribution','Gamma distribution','Exponential distribution','Location','east')

% plot on a log-scale
u=0.01;
du=0.00001;
semilogy(0:du:u,lognpdf(0:du:u,par_ln(1),par_ln(2)),'k','LineWidth',2)
hold on
semilogy(0:du:u,gampdf(0:du:u,par_gam(1),par_gam(2)),'r','LineWidth',2)
semilogy(0:du:u,exppdf(0:du:u,par_exp(1)),'b','LineWidth',2)
ylim([0,500])
xlim([0,u])
xlabel('Relapse rate')
ylabel('Density')
title('Distribution of the relapse rates in the population heterogeneity model')
legend('Lognormal distribution','Gamma distribution','Exponential distribution','Location','east')
