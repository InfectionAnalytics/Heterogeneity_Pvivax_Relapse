%%% Fitting all models to the first recurrence time in the PNG data by
%%% Robinson et al.


%%% Overview:
% 1. Fitting to all data
% 2. CIs for fitting to all data
% 3. Fitting to data by village
% 4. CIs for fitting to data by village


% This script uses the following functions:
% modelfit.m
% modelfit_1rec_1w.m
% modelfit_1rec_vill.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% 1. Fitting to all data %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data for time to first P. vivax infection (PCR positive):
% table1 = readtable('Albinama_Table1.csv');
% table1 = table1(:,2:15);
table2 = readtable('Albinama_Table2_Pv_PCR.csv');
table2 = table2(:,2:9);

% extract data of time to first P. vivax recurrence:
data = table2(:,[1,4,5,6,7]);
data.Properties.VariableNames = {'id','pq','village','event','time'};


% Model fit:
% fit models 1 to 4 with 1 reinfection rate, daily follow-up, and 10 risk
% groups for models 3 and 4
% use 1 drug washout distribution for all individuals regardless of
% treatment

n = 10; % number of relapse risk groups for model 3
f = 1; % use daily follow-up scheme

n_par = [4,5,5,6]; % number of estimated parameters for each model

k = 100; % number of optimisations
opts = optimoptions('fmincon','Display','iter','UseParallel',true);

for m=1:4
    if m==1
        nllh_1 = zeros(k,1);
        par_1 = zeros(k,n_par(m));
        
        fun =@(par) modelfit_1rec_1w(data,m,f,par);
        
        for i=1:k
            [par_tmp,nllh_tmp] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10],...
                [],[],[],[],[-Inf,0,0,0],[Inf,Inf,Inf,Inf],[],opts);
            nllh_1(i) = nllh_tmp;
            par_1(i,:) = par_tmp;
            disp(i)
        end
        [nllh_min1,ind_min1] = min(nllh_1);
        par_min1 = par_1(ind_min1,:);
        
    elseif m==2
        nllh_2 = zeros(k,1);
        par_2 = zeros(k,n_par(m));
        
        fun =@(par) modelfit_1rec_1w(data,m,f,par);
        
        warning off
        for i=1:k
            [par_tmp,nllh_tmp] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10],...
                [],[],[],[],[-Inf,0,0,0,0],[Inf,Inf,Inf,Inf,Inf],[],opts);
            nllh_2(i) = nllh_tmp;
            par_2(i,:) = par_tmp;
            disp(i)
        end
        [nllh_min2,ind_min2] = min(nllh_2);
        par_min2 = par_2(ind_min2,:);
        warning on
        
    elseif m==3
        nllh_3 = zeros(k,1);
        par_3 = zeros(k,n_par(m));
        
        fun =@(par) modelfit_1rec_1w(data,m,f,[par,n]);
        
        warning off
        for i=1:k
            [par_tmp,nllh_tmp] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10],...
                [],[],[],[],[-Inf,0,0,-Inf,0],[Inf,Inf,Inf,Inf,Inf],[],opts);
            nllh_3(i) = nllh_tmp;
            par_3(i,:) = par_tmp;
            disp(i)
            disp(nllh_tmp)
        end
        [nllh_min3,ind_min3] = min(nllh_3);
        par_min3 = par_3(ind_min3,:);
        warning on
        
    elseif m==4
        nllh_4 = zeros(k,1);
        par_4 = zeros(k,n_par(m));
        
        fun =@(par) modelfit_1rec_1w(data,m,f,[par,n]);
        
        warning off
        for i=1:k
            [par_tmp,nllh_tmp] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10,rand(1)*10],...
                [],[],[],[],[-Inf,0,0,-Inf,0,0],[Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
            nllh_4(i) = nllh_tmp;
            par_4(i,:) = par_tmp;
            disp(i)
        end
        [nllh_min4,ind_min4] = min(nllh_4);
        par_min4 = par_4(ind_min4,:);
        warning on
    end
    
end

%%% Visualize model fits with no relapses for PMQ: (Fig. 2, Fig. S2)

% Model fit results:
nllh_min1 = 1.4258e+03;
par_min1 = [3.968444968280213;1.814423923403087;0.044233859539881;0.003543223391707];
nllh_min2 = 1.402786092821062e+03;
par_min2 = [3.009398435588818;0.459727046381152;0.001932779709706;0.033289366409188;0.021559635498117];
nllh_min3 = 1.4113e+03;
par_min3 = [3.114275562021493;0.480188030482425;0.001925892692931;-4.728788522912486;2.261067653951518;10];
nllh_min4 = 1.402787391422036e+03;
par_min4 = [3.010736580936364;0.460333407270969;0.001933266689448;-3.399990688185307;0.001780871558982;...
    0.021586716594162;10];

n_par = [4,5,5,6]; % number of estimated parameters for each model

% Plot all models: 
t = 0:250;
U_pl_1 = modelfit(t,1,[par_min1(1),par_min1(2),par_min1(3)]);
U_pq_1 = modelfit(t,1,[par_min1(1),par_min1(2),par_min1(4)]);
U_pl_2 = modelfit(t,2,[par_min2(1),par_min2(2),par_min2(3),par_min2(4),par_min2(5)]);
U_pq_2 = modelfit(t,2,[par_min2(1),par_min2(2),par_min2(3),0,0]);
U_pl_3 = sum(modelfit(t,3,[par_min3(1),par_min3(2),par_min3(3),par_min3(6),par_min3(4),par_min3(5)]),1);
U_pq_3 = sum(modelfit(t,3,[par_min3(1),par_min3(2),par_min3(3),par_min3(6),0,0]),1);
U_pl_4 = sum(modelfit(t,4,[par_min4(1),par_min4(2),par_min4(3),par_min4(7),par_min4(4),par_min4(5),par_min4(6)]),1);
U_pq_4 = sum(modelfit(t,4,[par_min4(1),par_min4(2),par_min4(3),par_min4(7),0,0,0]),1);

subplot(2,2,1)
plot(t,U_pq_1,'r','LineWidth',3)
ylim([0 1.05])
title('Model 1: constant relapse rate')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_1,'b','LineWidth',3)
text(25,0.15,['AIC = ',num2str(round(2*nllh_min1+2*n_par(1)))])
hold off

subplot(2,2,2)
plot(t,U_pq_2,'r','LineWidth',3)
ylim([0 1.05])
title('Model 2: temporal heterogeneity')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_2,'b','LineWidth',3)
text(25,0.15,['AIC = ',num2str(round(2*nllh_min2+2*n_par(2)))])
hold off

subplot(2,2,3)
plot(t,U_pq_3,'r','LineWidth',3);
ylim([0 1.05])
title('Model 3: population heterogeneity')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_3,'b','LineWidth',3);
text(25,0.15,['AIC = ',num2str(round(2*nllh_min3+2*n_par(3)))])
hold off

subplot(2,2,4)
h1 = plot(t,U_pq_4,'r','LineWidth',3);
ylim([0 1.05])
title('Model 4: temporal & population heterogeneity')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
h2 = plot(t,U_pl_4,'b','LineWidth',3);
text(25,0.15,['AIC = ',num2str(round(2*nllh_min4+2*n_par(4)))])
hold off

% add legend
% hl = subplot(2,2,4);
% poshL = get(hl,'position');     % Getting its position
% lgd = legend(hl,[h1;h2],'Primaquine','Placebo');
% set(lgd,'position',poshL);      % Adjusting legend's position
% axis(hl,'off');                 % Turning its axis off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% 2. CIs for fitting to all data %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data for time to first P. vivax infection (PCR positive):
% table1 = readtable('Albinama_Table1.csv');
% table1 = table1(:,2:15);
table2 = readtable('Albinama_Table2_Pv_PCR.csv');
table2 = table2(:,2:9);

% extract data of time to first P. vivax recurrence:
data = table2(:,[1,4,5,6,7]);
data.Properties.VariableNames = {'id','pq','village','event','time'};

m = 4; % model (1, 2, 3, or 4)
n_par = [4,5,5,6]; % number of estimated parameters for each model
opts = optimoptions('fmincon','Display','iter','UseParallel',true);

n = 10; % number of relapse risk groups for models 3 and 4
f = 1; % use daily follow-up scheme

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
        fun =@(par) modelfit_1rec_1w(data,m,f,par);
        for j=1:(n_rand+1)
            if j==1 % use best fit as initial parameter values
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[3.968444968280213;1.814423923403087;...
                    0.044233859539881;0.003543223391707]',...
                    [],[],[],[],[-Inf,0,0,0],[Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            else
                pctRunOnAll warning off
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10],...
                    [],[],[],[],[-Inf,0,0,0],[Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            end
        end
    elseif m==2
        fun =@(par) modelfit_1rec_1w(data,m,f,par);
        for j=1:(n_rand+1)
            if j==1 % use best fit as initial parameter values
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[3.009398435588818;0.459727046381152;...
                    0.001932779709706;0.033289366409188;0.021559635498117]',...
                    [],[],[],[],[-Inf,0,0,0,0],[Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            else
                pctRunOnAll warning off
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10],...
                    [],[],[],[],[-Inf,0,0,0,0],[Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            end
        end
    elseif m==3
        fun =@(par) modelfit_1rec_1w(data,m,f,[par,n]);
        for j=1:(n_rand+1)
            if j==1 % use best fit as initial parameter values
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[3.114275562021493;0.480188030482425;0.001925892692931;...
                    -4.728788522912486;2.261067653951518]',...
                    [],[],[],[],[-Inf,0,0,-Inf,0],[Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            else
                pctRunOnAll warning off
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10],...
                    [],[],[],[],[-Inf,0,0,-Inf,0],[Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            end
        end
    elseif m==4
        fun =@(par) modelfit_1rec_1w(data,m,f,[par,n]);
        for j=1:(n_rand+1)
            if j==1 % use best fit as initial parameter values
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[3.010736580936364;0.460333407270969;0.001933266689448;...
                    -3.399990688185307;0.001780871558982;0.021586716594162]',...
                    [],[],[],[],[-Inf,0,0,-Inf,0,0],[Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            else
                pctRunOnAll warning off
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10,rand(1)*10],...
                    [],[],[],[],[-Inf,0,0,-Inf,0,0],[Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
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
    save('b-CI(11opt)m1-PNG-all.mat','nllhs_all','param_all')
elseif m==2
    save('b-CI(11opt)m2-PNG-all.mat','nllhs_all','param_all')
elseif m==3
    save('b-CI(11opt)m3-PNG-all.mat','nllhs_all','param_all')
elseif m==4
    save('b-CI(11opt)m4-PNG-all.mat','nllhs_all','param_all')
end


% load bootstrapped data:
load('b-CI(11opt)m1-PNG-all.mat')
param1 = param_all;
load('b-CI(11opt)m2-PNG-all.mat')
param2 = param_all;
load('b-CI(11opt)m3-PNG-all.mat')
param3 = param_all;
load('b-CI(11opt)m4-PNG-all.mat')
param4 = param_all;

% CI using the percentile method:
alpha = 95;
for m=1:4 % for each model
    if m==1
        ci_1 = zeros(size(param1,2),2);
        for i=1:size(param1,2)
            ci_1(i,1) = prctile(param1(:,i),(100-alpha)/2);
            ci_1(i,2) = prctile(param1(:,i),alpha+(100-alpha)/2);
        end
    elseif m==2
        ci_2 = zeros(size(param2,2),2);
        for i=1:size(param2,2)
            ci_2(i,1) = prctile(param2(:,i),(100-alpha)/2);
            ci_2(i,2) = prctile(param2(:,i),alpha+(100-alpha)/2);
        end
    elseif m==3
        ci_3 = zeros(size(param3,2),2);
        for i=1:size(param3,2)
            ci_3(i,1) = prctile(param3(:,i),(100-alpha)/2);
            ci_3(i,2) = prctile(param3(:,i),alpha+(100-alpha)/2);
        end
    elseif m==4
        ci_4 = zeros(size(param4,2),2);
        for i=1:size(param4,2)
            ci_4(i,1) = prctile(param4(:,i),(100-alpha)/2);
            ci_4(i,2) = prctile(param4(:,i),alpha+(100-alpha)/2);
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% 3. Fitting to data by village %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fit models 1 to 4 to the data with the same drug washout distribution for
% all villages but different reinfection rates and relapse time distribution

% load data for time to first P. vivax infection (PCR positive):
table2 = readtable('Albinama_Table2_Pv_PCR.csv');
table2 = table2(:,2:9);

% extract data of time to first P. vivax recurrence:
data = table2(:,[1,4,5,6,7]);
data.Properties.VariableNames = {'id','pq','village','event','time'};

% Model fits:
% fit models 1 to 4 with 1 reinfection rate, daily follow-up, and 10 risk
% groups for models 3 and 4
n = 10; % number of relapse risk groups for model 3
f = 1; % use daily follow-up scheme

% Parameters for the model fits:
% Model 1: 2 parameters for drug washout distribution, recurrence rates for
% PQ and PL treated patients for each of the 5 villages, i.e., 10
% recurrence rates => overall 12 parameters
% Model 2: 2 parameters for drug washout distribution, reinfection rate,
% intitial relapse rate, and decay of the relapse rate for each village =>
% overall 17 parameters
% Model 3: 2 parameters for drug washout distribution, reinfection rates
% and 2 parameters for the relapse rate distribution for each village =>
% overall 17 parameters
% Model 4: 2 parameters for drug washout distribution, reinfection rates, 2
% parameters for the relapse rate distribution, and a decay rate of the
% relapse rate for each village => overall 22 parameters
n_par = [12,17,17,22]; % number of estimated parameters for each model

k = 100; % number of optimisations
opts = optimoptions('fmincon','Display','iter','UseParallel',true);

for m=1:4
    if m==1
        nllh_1 = zeros(k,1);
        par_1 = zeros(k,n_par(m));
        
        fun =@(par) modelfit_1rec_vill(data,m,f,par);
        
        for i=1:k
            [par_tmp,nllh_tmp] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,...
                rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,...
                rand(1)*10],[],[],[],[],[-Inf,0,0,0,0,0,0,0,0,0,0,0],...
                [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
            nllh_1(i) = nllh_tmp;
            par_1(i,:) = par_tmp;
            disp(i)
        end
        [nllh_min1,ind_min1] = min(nllh_1);
        par_min1 = par_1(ind_min1,:);
        
    elseif m==2
        nllh_2 = zeros(k,1);
        par_2 = zeros(k,n_par(m));
        
        fun =@(par) modelfit_1rec_vill(data,m,f,par);
        
        warning off
        for i=1:k
            [par_tmp,nllh_tmp] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,...
                rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,...
                rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10],[],[],[],[],...
                [-Inf,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,...
                Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
            nllh_2(i) = nllh_tmp;
            par_2(i,:) = par_tmp;
            disp(i)
        end
        [nllh_min2,ind_min2] = min(nllh_2);
        par_min2 = par_2(ind_min2,:);
        warning on
        
    elseif m==3
        nllh_3 = zeros(k,1);
        par_3 = zeros(k,n_par(m));
        
        fun =@(par) modelfit_1rec_vill(data,m,f,[par,n]);
        
        warning off
        for i=1:k
            [par_tmp,nllh_tmp] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10,...
                rand(1)*10,rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10,rand(1)*10,...
                rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10],[],[],[],[],...
                [-Inf,0,0,-Inf,0,0,-Inf,0,0,-Inf,0,0,-Inf,0,0,-Inf,0],[Inf,Inf,Inf,Inf,...
                Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
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
        
        fun =@(par) modelfit_1rec_vill(data,m,f,[par,n]);
        
        warning off
        for i=1:k
            [par_tmp,nllh_tmp] = fmincon(fun,[rand(1)*10-5,rand(1)*10,...
                rand(1)*10,rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10,rand(1)*10,...
                rand(1)*10,rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10,rand(1)*10,...
                rand(1)*10,rand(1)*10-5,rand(1)*10,rand(1)*10],[],[],[],[],...
                [-Inf,0,0,-Inf,0,0,0,-Inf,0,0,0,-Inf,0,0,0,-Inf,0,0,0,-Inf,0,0],...
                [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
            nllh_4(i) = nllh_tmp;
            par_4(i,:) = par_tmp;
            disp(i)
        end
        [nllh_min4,ind_min4] = min(nllh_4);
        par_min4 = par_4(ind_min4,:);
        warning on
    end
    
end


%%% Visualize model prediction for each village (Fig. S3, Fig. S4)

% Model fit results:
n_par = [12,17,17,22];
n = 10;

% Parameters:
nllh_min1 = 1.359768351566670e+03;
par_min1 = [3.700155427262943,1.443013215269731,0.028470860431871,0.003141683070781,...
    0.007048157709338,0.000144231673975,0.054764786765070,0.009216807968129,...
    0.043855388547067,0.001273162370849,0.124504939113478,0.030839662925882];
nllh_min2 = 1.320622721416308e+03;
par_min2 = [3.153300222608564,0.489190433336640,0.002020264562845,0.027094026113808,...
    0.018616586301897,0.000110104363797,0.008964252270973,0.008006867974303,...
    0.005068228319300,0.043136249637375,0.010873198337848,0.000911910956375,...
    0.044135113124644,0.022153594732631,0.009981864178311,0.241979319347036,...
    0.014764582609550];
nllh_min3 = 1.325515162253067e+03;
par_min3 = [3.195952269453012,0.488718857944258,0.002073439849889,-5.029239751189614,...
    1.845373879451467,0.000107921583196,-5.694800411642519,1.632290670765690,...
    0.005095707640667,-3.574578993884039,0.935013179874915,0.000881235324295,...
    -4.393291614830225,2.258679027251658,0.010107520276450,-1.235438233197847,...
    1.127999032231824];
nllh_min4 = 1.320488747060091e+03;
par_min4 = [3.179119144463923;0.509942679230476;0.002031181165949;-3.558805481375394;...
    0.014631662468136;0.019007870501204;0.000108849467981;-5.611510001555497;1.555147740784638;...
    0.000743807140125;0.005107254227906;-3.081461402702949;0.037684818068242;0.011308476686592;...
    0.000918381996604;-3.067167304398428;0.041320843011641;0.022612271608578;0.010099986083958;...
    -1.364601710534829;0.085097209059369;0.011397905583331];

t = 0:230;
% Model 1 fit for different villages:
U_pl_1_1 = modelfit(t,1,[par_min1(1),par_min1(2),par_min1(3)]);
U_pq_1_1 = modelfit(t,1,[par_min1(1),par_min1(2),par_min1(4)]);
U_pl_1_2 = modelfit(t,1,[par_min1(1),par_min1(2),par_min1(5)]);
U_pq_1_2 = modelfit(t,1,[par_min1(1),par_min1(2),par_min1(6)]);
U_pl_1_3 = modelfit(t,1,[par_min1(1),par_min1(2),par_min1(7)]);
U_pq_1_3 = modelfit(t,1,[par_min1(1),par_min1(2),par_min1(8)]);
U_pl_1_4 = modelfit(t,1,[par_min1(1),par_min1(2),par_min1(9)]);
U_pq_1_4 = modelfit(t,1,[par_min1(1),par_min1(2),par_min1(10)]);
U_pl_1_5 = modelfit(t,1,[par_min1(1),par_min1(2),par_min1(11)]);
U_pq_1_5 = modelfit(t,1,[par_min1(1),par_min1(2),par_min1(12)]);
% Model 2 fit for different villages:
U_pl_2_1 = modelfit(t,2,[par_min2(1),par_min2(2),par_min2(3),par_min2(4),par_min2(5)]);
U_pq_2_1 = modelfit(t,2,[par_min2(1),par_min2(2),par_min2(3),0,0]);
U_pl_2_2 = modelfit(t,2,[par_min2(1),par_min2(2),par_min2(6),par_min2(7),par_min2(8)]);
U_pq_2_2 = modelfit(t,2,[par_min2(1),par_min2(2),par_min2(6),0,0]);
U_pl_2_3 = modelfit(t,2,[par_min2(1),par_min2(2),par_min2(9),par_min2(10),par_min2(11)]);
U_pq_2_3 = modelfit(t,2,[par_min2(1),par_min2(2),par_min2(9),0,0]);
U_pl_2_4 = modelfit(t,2,[par_min2(1),par_min2(2),par_min2(12),par_min2(13),par_min2(14)]);
U_pq_2_4 = modelfit(t,2,[par_min2(1),par_min2(2),par_min2(12),0,0]);
U_pl_2_5 = modelfit(t,2,[par_min2(1),par_min2(2),par_min2(15),par_min2(16),par_min2(17)]);
U_pq_2_5 = modelfit(t,2,[par_min2(1),par_min2(2),par_min2(15),0,0]);
% Model 3 fit for different villages:
U_pl_3_1 = sum(modelfit(t,3,[par_min3(1),par_min3(2),par_min3(3),n,par_min3(4),par_min3(5)]),1);
U_pq_3_1 = sum(modelfit(t,3,[par_min3(1),par_min3(2),par_min3(3),n,0,0]),1);
U_pl_3_2 = sum(modelfit(t,3,[par_min3(1),par_min3(2),par_min3(6),n,par_min3(7),par_min3(8)]),1);
U_pq_3_2 = sum(modelfit(t,3,[par_min3(1),par_min3(2),par_min3(6),n,0,0]),1);
U_pl_3_3 = sum(modelfit(t,3,[par_min3(1),par_min3(2),par_min3(9),n,par_min3(10),par_min3(11)]),1);
U_pq_3_3 = sum(modelfit(t,3,[par_min3(1),par_min3(2),par_min3(9),n,0,0]),1);
U_pl_3_4 = sum(modelfit(t,3,[par_min3(1),par_min3(2),par_min3(12),n,par_min3(13),par_min3(14)]),1);
U_pq_3_4 = sum(modelfit(t,3,[par_min3(1),par_min3(2),par_min3(12),n,0,0]),1);
U_pl_3_5 = sum(modelfit(t,3,[par_min3(1),par_min3(2),par_min3(15),n,par_min3(16),par_min3(17)]),1);
U_pq_3_5 = sum(modelfit(t,3,[par_min3(1),par_min3(2),par_min3(15),n,0,0]),1);
% Model 4 fit for different villages:
U_pl_4_1 = sum(modelfit(t,4,[par_min4(1),par_min4(2),par_min4(3),n,par_min4(4),par_min4(5),par_min4(6)]),1);
U_pq_4_1 = sum(modelfit(t,4,[par_min4(1),par_min4(2),par_min4(3),n,0,0,0]),1);
U_pl_4_2 = sum(modelfit(t,4,[par_min4(1),par_min4(2),par_min4(7),n,par_min4(8),par_min4(9),par_min4(10)]),1);
U_pq_4_2 = sum(modelfit(t,4,[par_min4(1),par_min4(2),par_min4(7),n,0,0,0]),1);
U_pl_4_3 = sum(modelfit(t,4,[par_min4(1),par_min4(2),par_min4(11),n,par_min4(12),par_min4(13),par_min4(14)]),1);
U_pq_4_3 = sum(modelfit(t,4,[par_min4(1),par_min4(2),par_min4(11),n,0,0,0]),1);
U_pl_4_4 = sum(modelfit(t,4,[par_min4(1),par_min4(2),par_min4(15),n,par_min4(16),par_min4(17),par_min4(18)]),1);
U_pq_4_4 = sum(modelfit(t,4,[par_min4(1),par_min4(2),par_min4(15),n,0,0,0]),1);
U_pl_4_5 = sum(modelfit(t,4,[par_min4(1),par_min4(2),par_min4(19),n,par_min4(20),par_min4(21),par_min4(22)]),1);
U_pq_4_5 = sum(modelfit(t,4,[par_min4(1),par_min4(2),par_min4(19),n,0,0,0]),1);


% Plot of models 1 to 3 (Fig. S3): (model 4 is plotted separately below)
sgtitle('Model fit to data grouped by village')
subplot(5,3,1)
plot(t,U_pq_1_1,'b','LineWidth',3)
ylim([0 1.05])
title('Village cluster 1, Model 1')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_1_1,'r','LineWidth',3)
text(10,0.1,['AIC = ',num2str(round(2*nllh_min1+2*n_par(1)))])
hold off

subplot(5,3,2)
plot(t,U_pq_2_1,'b','LineWidth',3)
ylim([0 1.05])
title('Village cluster 1, Model 2')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_2_1,'r','LineWidth',3)
text(10,0.1,['AIC = ',num2str(round(2*nllh_min2+2*n_par(2)))])
hold off

subplot(5,3,3)
plot(t,U_pq_3_1,'b','LineWidth',3)
ylim([0 1.05])
title('Village cluster 1, Model 3')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_3_1,'r','LineWidth',3)
text(10,0.1,['AIC = ',num2str(round(2*nllh_min3+2*n_par(3)))])
hold off

subplot(5,3,4)
plot(t,U_pq_1_2,'b','LineWidth',3)
ylim([0 1.05])
title('Village cluster 2, Model 1')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_1_2,'r','LineWidth',3)
% text(150,0.75,['AIC = ',num2str(round(2*nllh_min3(1)+2*n_par(m)))])
hold off

subplot(5,3,5)
plot(t,U_pq_2_2,'b','LineWidth',3)
ylim([0 1.05])
title('Village cluster 2, Model 2')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_2_2,'r','LineWidth',3)
% text(150,0.75,['AIC = ',num2str(round(2*nllh_min3(1)+2*n_par(m)))])
hold off

subplot(5,3,6)
plot(t,U_pq_3_2,'b','LineWidth',3)
ylim([0 1.05])
title('Village cluster 2, Model 3')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_3_2,'r','LineWidth',3)
% text(150,0.75,['AIC = ',num2str(round(2*nllh_min3(1)+2*n_par(m)))])
hold off

subplot(5,3,7)
plot(t,U_pq_1_3,'b','LineWidth',3)
ylim([0 1.05])
title('Village cluster 3, Model 1')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_1_3,'r','LineWidth',3)
% text(150,0.75,['AIC = ',num2str(round(2*nllh_min3(1)+2*n_par(m)))])
hold off

subplot(5,3,8)
plot(t,U_pq_2_3,'b','LineWidth',3)
ylim([0 1.05])
title('Village cluster 3, Model 2')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_2_3,'r','LineWidth',3)
% text(150,0.75,['AIC = ',num2str(round(2*nllh_min3(1)+2*n_par(m)))])
hold off

subplot(5,3,9)
plot(t,U_pq_3_3,'b','LineWidth',3)
ylim([0 1.05])
title('Village cluster 3, Model 3')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_3_3,'r','LineWidth',3)
% text(150,0.75,['AIC = ',num2str(round(2*nllh_min3(1)+2*n_par(m)))])
hold off

subplot(5,3,10)
plot(t,U_pq_1_4,'b','LineWidth',3)
ylim([0 1.05])
title('Village cluster 4, Model 1')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_1_4,'r','LineWidth',3)
% text(150,0.75,['AIC = ',num2str(round(2*nllh_min3(1)+2*n_par(m)))])
hold off

subplot(5,3,11)
plot(t,U_pq_2_4,'b','LineWidth',3)
ylim([0 1.05])
title('Village cluster 4, Model 2')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_2_4,'r','LineWidth',3)
% text(150,0.75,['AIC = ',num2str(round(2*nllh_min3(1)+2*n_par(m)))])
hold off

subplot(5,3,12)
plot(t,U_pq_3_4,'b','LineWidth',3)
ylim([0 1.05])
title('Village cluster 4, Model 3')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_3_4,'r','LineWidth',3)
% text(150,0.75,['AIC = ',num2str(round(2*nllh_min3(1)+2*n_par(m)))])
hold off

subplot(5,3,13)
plot(t,U_pq_1_5,'b','LineWidth',3)
ylim([0 1.05])
title('Village cluster 5, Model 1')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_1_5,'r','LineWidth',3)
% text(150,0.75,['AIC = ',num2str(round(2*nllh_min3(1)+2*n_par(m)))])
hold off

subplot(5,3,14)
plot(t,U_pq_2_5,'b','LineWidth',3)
ylim([0 1.05])
title('Village cluster 5, Model 2')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_2_5,'r','LineWidth',3)
% text(150,0.75,['AIC = ',num2str(round(2*nllh_min3(1)+2*n_par(m)))])
hold off

subplot(5,3,15)
plot(t,U_pq_3_5,'b','LineWidth',3)
ylim([0 1.05])
title('Village cluster 5, Model 3')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_3_5,'r','LineWidth',3)
% text(150,0.75,['AIC = ',num2str(round(2*nllh_min3(1)+2*n_par(m)))])
hold off

% Plot model 4 separately: (Fig. S4)
sgtitle('Model 4 (temporal and population heterogeneity) fit to data grouped by village')
subplot(2,3,1)
plot(t,U_pq_4_1,'b','LineWidth',3)
ylim([0 1.05])
title('Village cluster 1')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_4_1,'r','LineWidth',3)
text(10,0.1,['AIC = ',num2str(round(2*nllh_min4+2*n_par(4)))])
hold off

subplot(2,3,2)
plot(t,U_pq_4_2,'b','LineWidth',3)
ylim([0 1.05])
title('Village cluster 2')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_4_2,'r','LineWidth',3)
% text(10,0.1,['AIC = ',num2str(round(2*nllh_min4+2*n_par(4)))])
hold off

subplot(2,3,3)
plot(t,U_pq_4_3,'b','LineWidth',3)
ylim([0 1.05])
title('Village cluster 3')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_4_3,'r','LineWidth',3)
% text(10,0.1,['AIC = ',num2str(round(2*nllh_min4+2*n_par(4)))])
hold off

subplot(2,3,4)
plot(t,U_pq_4_4,'b','LineWidth',3)
ylim([0 1.05])
title('Village cluster 4')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_4_4,'r','LineWidth',3)
% text(10,0.1,['AIC = ',num2str(round(2*nllh_min4+2*n_par(4)))])
hold off

subplot(2,3,5)
plot(t,U_pq_4_5,'b','LineWidth',3)
ylim([0 1.05])
title('Village cluster 5')
xlabel('Time [days]')
ylabel('Fraction without recurrence')
hold on
plot(t,U_pl_4_5,'r','LineWidth',3)
% text(10,0.1,['AIC = ',num2str(round(2*nllh_min4+2*n_par(4)))])
hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% 4. CIs for fitting to data by village %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data for time to first P. vivax infection (PCR positive):
% table1 = readtable('Albinama_Table1.csv');
% table1 = table1(:,2:15);
table2 = readtable('Albinama_Table2_Pv_PCR.csv');
table2 = table2(:,2:9);

% extract data of time to first P. vivax recurrence:
data = table2(:,[1,4,5,6,7]);
data.Properties.VariableNames = {'id','pq','village','event','time'};

m = 4; % model (1, 2, 3, or 4)
n_par = [12,17,17,22]; % number of estimated parameters for each model
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
        fun =@(par) modelfit_1rec_vill(data,m,f,par);
        for j=1:(n_rand+1)
            if j==1 % use best fit as initial parameter values
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[3.700155427262943,1.443013215269731,0.028470860431871,...
                    0.003141683070781,0.007048157709338,0.000144231673975,0.054764786765070,0.009216807968129,...
                    0.043855388547067,0.001273162370849,0.124504939113478,0.030839662925882],...
                    [],[],[],[],[-Inf,0,0,0,0,0,0,0,0,0,0,0],...
                    [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            else
                pctRunOnAll warning off
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,...
                    rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,...
                    rand(1)*10],[],[],[],[],[-Inf,0,0,0,0,0,0,0,0,0,0,0],...
                    [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            end
        end
    elseif m==2
        fun =@(par) modelfit_1rec_vill(data,m,f,par);
        for j=1:(n_rand+1)
            if j==1 % use best fit as initial parameter values
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[3.153300222608564,0.489190433336640,0.002020264562845,...
                    0.027094026113808,0.018616586301897,0.000110104363797,0.008964252270973,0.008006867974303,...
                    0.005068228319300,0.043136249637375,0.010873198337848,0.000911910956375,0.044135113124644,...
                    0.022153594732631,0.009981864178311,0.241979319347036,0.014764582609550],[],[],[],[],...
                    [-Inf,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,...
                    Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            else
                pctRunOnAll warning off
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,...
                    rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,...
                    rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10],[],[],[],[],...
                    [-Inf,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,...
                    Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            end
        end
    elseif m==3
        fun =@(par) modelfit_1rec_vill(data,m,f,[par,n]);
        for j=1:(n_rand+1)
            if j==1 % use best fit as initial parameter values
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[3.195952269453012,0.488718857944258,0.002073439849889,...
                    -5.029239751189614,1.845373879451467,0.000107921583196,-5.694800411642519,1.632290670765690,...
                    0.005095707640667,-3.574578993884039,0.935013179874915,0.000881235324295,-4.393291614830225,...
                    2.258679027251658,0.010107520276450,-1.235438233197847,1.127999032231824],[],[],[],[],...
                    [-Inf,0,0,-Inf,0,0,-Inf,0,0,-Inf,0,0,-Inf,0,0,-Inf,0],[Inf,Inf,Inf,Inf,...
                    Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            else
                pctRunOnAll warning off
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10,...
                    rand(1)*10,rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10,rand(1)*10,...
                    rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10],[],[],[],[],...
                    [-Inf,0,0,-Inf,0,0,-Inf,0,0,-Inf,0,0,-Inf,0,0,-Inf,0],[Inf,Inf,Inf,Inf,...
                    Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            end
        end
    elseif m==4
        fun =@(par) modelfit_1rec_vill(data,m,f,[par,n]);
        for j=1:(n_rand+1)
            if j==1 % use best fit as initial parameter values
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[3.179119144463923;0.509942679230476;0.002031181165949;...
                    -3.558805481375394;0.014631662468136;0.019007870501204;0.000108849467981;-5.611510001555497;...
                    1.555147740784638;0.000743807140125;0.005107254227906;-3.081461402702949;0.037684818068242;...
                    0.011308476686592;0.000918381996604;-3.067167304398428;0.041320843011641;0.022612271608578;...
                    0.010099986083958;-1.364601710534829;0.085097209059369;0.011397905583331]',[],[],[],[],...
                    [-Inf,0,0,-Inf,0,0,0,-Inf,0,0,0,-Inf,0,0,0,-Inf,0,0,0,-Inf,0,0],...
                    [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
                warning on                
                nllh_tmp(j) = nllh_tmp2;
                par_tmp(j,:) = par_tmp2;
            else
                pctRunOnAll warning off
                warning('off','all')
                [par_tmp2,nllh_tmp2] = fmincon(fun,[rand(1)*10-5,rand(1)*10,...
                    rand(1)*10,rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10,rand(1)*10,...
                    rand(1)*10,rand(1)*10-5,rand(1)*10,rand(1)*10,rand(1)*10,rand(1)*10-5,rand(1)*10,rand(1)*10,...
                    rand(1)*10,rand(1)*10-5,rand(1)*10,rand(1)*10],[],[],[],[],...
                    [-Inf,0,0,-Inf,0,0,0,-Inf,0,0,0,-Inf,0,0,0,-Inf,0,0,0,-Inf,0,0],...
                    [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf],[],opts);
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
    save('b-CI(11opt)m1-PNG-vill.mat','nllhs_all','param_all')
elseif m==2
    save('b-CI(11opt)m2-PNG-vill.mat','nllhs_all','param_all')
elseif m==3
    save('b-CI(11opt)m3-PNG-vill.mat','nllhs_all','param_all')
elseif m==4
    save('b-CI(11opt)m4-PNG-vill.mat','nllhs_all','param_all')
end


% load bootstrapped data:
load('b-CI(11opt)m1-PNG-vill.mat')
param1 = param_all;
load('b-CI(11opt)m2-PNG-vill.mat')
param2 = param_all;
load('b-CI(11opt)m3-PNG-vill.mat')
param3 = param_all;
load('b-CI(11opt)m4-PNG-vill.mat')
param4 = param_all;

% CI using the percentile method:
alpha = 95;
for m=1:4 % for each model
    if m==1
        ci_1 = zeros(size(param1,2),2);
        for i=1:size(param1,2)
            ci_1(i,1) = prctile(param1(:,i),(100-alpha)/2);
            ci_1(i,2) = prctile(param1(:,i),alpha+(100-alpha)/2);
        end
    elseif m==2
        ci_2 = zeros(size(param2,2),2);
        for i=1:size(param2,2)
            ci_2(i,1) = prctile(param2(:,i),(100-alpha)/2);
            ci_2(i,2) = prctile(param2(:,i),alpha+(100-alpha)/2);
        end
    elseif m==3
        ci_3 = zeros(size(param3,2),2);
        for i=1:size(param3,2)
            ci_3(i,1) = prctile(param3(:,i),(100-alpha)/2);
            ci_3(i,2) = prctile(param3(:,i),alpha+(100-alpha)/2);
        end
    elseif m==4
        ci_4 = zeros(size(param4,2),2);
        for i=1:size(param4,2)
            ci_4(i,1) = prctile(param4(:,i),(100-alpha)/2);
            ci_4(i,2) = prctile(param4(:,i),alpha+(100-alpha)/2);
        end
    end
end


