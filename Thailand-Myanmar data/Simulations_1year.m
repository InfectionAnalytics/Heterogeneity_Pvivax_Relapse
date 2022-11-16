%%% Simulating the models for a 1 year


% Simulate a population of individuals with reinfection and relapse
% rates/distributions determined by the best model fit parameters to the
% data for first and second recurrence


% Parameter values for reinf = 2; f = 1; n = 10;
par = table('Size',[16,3],'VariableTypes',{'double','double','double'},...
    'VariableNames',{'model1','model2','model3'},...
    'rowNames',{'washout1_AS','washout2_AS','washout1_CHQ','washout2_CHQ','washout1_DP',...
    'washout2_DP','recurrence_AS/CHQ_VHX','recurrence_PMQ_VHX','recurrence_PMQ_BPD',...
    'reinfection_VHX','reinfection_BPD','relapse_risk1_AS/CHQ','relapse_risk2_AS/CHQ',...
    'number_relapse_groups','I_AS/CHQ','d'});
par{1:9,1} = [3.44260377;1.61060976;4.051098296;1.15040291;3.88587443;0.05723774;0.08141406;0.001056392;0.00064541822];
par{[1:6,10,11,15:16],2} = [2.7914;0.23686;3.3847;0.3071;3.882;0.056692;0.00097171;0.00056813;0.070277;0.024911];
par{[1:6,10:14],3} = [2.95051812;0.24074922;3.53368025;0.32837685;3.88237672;0.05992895;0.00082751;0.00057377;-3.72783179;2.80439599;10];


% other parameter values:
n_sim = 1000; % number of simulations of n_ind individuals for each model
n_ind = 1000; % number of individuals to be simulated
ind_drug_study = "AS_VHX"; % drug and study for the simulated individuals
% drug and study combinations:
% AS_VHX, CHQ_VHX, CHQ/PMQ_VHX, CHQ/PMQ_BPD, DP/PMQ_BPD
tmax = 365; % simulate individuals to tmax, i.e. for 1 year

% simulation:
sim1 = cell(n_sim,1);
sim2 = cell(n_sim,1);
sim3 = cell(n_sim,1);
for l=1:n_sim % for each simulation
    for i=1:n_ind % for each individual
        % for each individual save number of individual (as and ID), time
        % of event, relapse (0 or 1), and event number
        ind_tmp_1 = [i,0,0,0]; % if the individual has no recurrences
        ind_tmp_2 = [i,0,0,0];
        ind_tmp_3 = [i,0,0,0];
        for m=1:3
            % For model 3: draw the individual's relapse rate from a
            % lognormal distribution. The individual keeps the same relapse
            % rate for the entire simulation (1 year).
            % If the individual is treated with primaquine, then the
            % relapse rate is 0 as all liver-stage parasites are assumed to
            % be killed by primaquine.
            if m==3
                if ind_drug_study=="AS_VHX"
                    rel_rate = lognrnd(par{12,m},par{13,m},1);
                elseif ind_drug_study=="CHQ_VHX"
                    rel_rate = lognrnd(par{12,m},par{13,m},1);
                elseif ind_drug_study=="CHQ/PMQ_VHX"
                    rel_rate = 0;
                elseif ind_drug_study=="CHQ/PMQ_BPD"
                    rel_rate = 0;
                elseif ind_drug_study=="DP/PMQ_BPD"
                    rel_rate = 0;
                end
            end
            cur_time = 0; % keep track of the time to simulate to tmax
            while cur_time<tmax
                if m==1
                    % Drug washout time for each individual: draw a random
                    % drug washout time from a lognormal distribution for
                    % each recurrence and treatment.
                    if ind_drug_study=="AS_VHX"
                        washout_time = lognrnd(par{1,m},par{2,m},1);
                    elseif ind_drug_study=="DP/PMQ_BPD"
                        washout_time = lognrnd(par{5,m},par{6,m},1);
                    else % "CHQ_VHX", "CHQ/PMQ_VHX", or "CHQ/PMQ_BPD"
                        washout_time = lognrnd(par{3,m},par{4,m},1);
                    end
                    % For the recurrence, we assume that the recurrence
                    % rate of primaquine treated individuals is the
                    % reinfection rate (as those individuals are assumed to
                    % not have any relapses). Thus, the relapse rate is the
                    % recurrence rate minus the recurrence rate of
                    % primaquine treated patients.
                    if ismember(ind_drug_study,["AS_VHX","CHQ_VHX"])
                        reinf_time = washout_time + exprnd(1/par{8,m},1);
                        rel_time = washout_time + exprnd(1/(par{7,m}-par{8,m}),1);
                        recurrence_time = min([reinf_time,rel_time]);
                        rel = rel_time < reinf_time;
                    elseif ind_drug_study=="CHQ/PMQ_VHX"
                        recurrence_time = washout_time + exprnd(1/par{8,m},1);
                        rel = 0;
                    else
                        recurrence_time = washout_time + exprnd(1/par{9,m},1);
                        rel = 0;
                    end
                    if cur_time + recurrence_time <= tmax
                        if cur_time==0
                            ind_tmp_1 = [i,recurrence_time,rel,1];
                        else
                            ind_tmp_1 = [ind_tmp_1;i,recurrence_time,rel,max(ind_tmp_1(:,4))+1];
                        end
                    end
                    cur_time = cur_time + recurrence_time;
                elseif m==2
                    % For the temporal heterogeneity model, we use inverse
                    % transform sampling to simulate relapses (for details
                    % see the supplementary methods).
                    t = 0:0.01:(tmax-cur_time);
                    x = rand(1); % random number between 0 and 1
                    if ind_drug_study=="AS_VHX"
                        washout_time = lognrnd(par{1,m},par{2,m},1);
                        reinf_time = exprnd(1./par{10,m});
                        % cdf of the time to next relapse distribution:
                        y = 1-exp(-par{15,m}.*exp(-par{16,m}*washout_time).*(1-exp(-par{16,m}.*t))./par{16,m});
                    elseif ind_drug_study=="CHQ_VHX"
                        washout_time = lognrnd(par{3,m},par{4,m},1);
                        reinf_time = exprnd(1./par{10,m});
                        y = 1-exp(-par{15,m}.*exp(-par{16,m}*washout_time).*(1-exp(-par{16,m}.*t))./par{16,m});
                    elseif ind_drug_study=="CHQ/PMQ_VHX"
                        washout_time = lognrnd(par{3,m},par{4,m},1);
                        reinf_time = exprnd(1./par{10,m});
                    elseif ind_drug_study=="CHQ/PMQ_BPD"
                        washout_time = lognrnd(par{3,m},par{4,m},1);
                        reinf_time = exprnd(1./par{11,m});
                    elseif ind_drug_study=="DP/PMQ_BPD"
                        washout_time = lognrnd(par{5,m},par{6,m},1);
                        reinf_time = exprnd(1./par{11,m});
                    end
                    if ismember(ind_drug_study,["AS_VHX","CHQ_VHX"])
                        [~,tmp] = min(abs(y-x));
                        if all(y==0) || x>max(y)
                            rel_time = inf;
                        else
                            rel_time = t(tmp);
                        end
                    else % no relapses for primaquine treated patients
                        rel_time = inf;
                    end
                    recurrence_time = washout_time + min([reinf_time,rel_time]);
                    if cur_time + recurrence_time <= tmax
                        if cur_time==0
                            rel = (rel_time<reinf_time);
                            ind_tmp_2 = [i,recurrence_time,rel,1];
                        else
                            rel = (rel_time<reinf_time);
                            ind_tmp_2 = [ind_tmp_2;i,recurrence_time,rel,max(ind_tmp_2(:,4))+1];
                        end
                    end
                    cur_time = cur_time + recurrence_time;
                elseif m==3
                    % Drug washout time and reinfection time:
                    if ind_drug_study=="AS_VHX"
                        washout_time = lognrnd(par{1,m},par{2,m},1);
                        reinf_time = exprnd(1./par{10,m});
                    elseif ind_drug_study=="CHQ_VHX"
                        washout_time = lognrnd(par{3,m},par{4,m},1);
                        reinf_time = exprnd(1./par{10,m});
                    elseif ind_drug_study=="CHQ/PMQ_VHX"
                        washout_time = lognrnd(par{3,m},par{4,m},1);
                        reinf_time = exprnd(1./par{10,m});
                    elseif ind_drug_study=="CHQ/PMQ_BPD"
                        washout_time = lognrnd(par{3,m},par{4,m},1);
                        reinf_time = exprnd(1./par{11,m});
                    elseif ind_drug_study=="DP/PMQ_BPD"
                        washout_time = lognrnd(par{5,m},par{6,m},1);
                        reinf_time = exprnd(1./par{11,m});
                    end
                    % Relapse time:
                    rel_time = exprnd(1./rel_rate);
                    recurrence_time = washout_time + min([reinf_time,rel_time]);
                    if cur_time + recurrence_time <= tmax
                        if cur_time==0
                            rel = (rel_time<reinf_time);
                            ind_tmp_3 = [i,recurrence_time,rel,1];
                        else
                            rel = (rel_time<reinf_time);
                            ind_tmp_3 = [ind_tmp_3;i,recurrence_time,rel,max(ind_tmp_3(:,4))+1];
                        end
                    end
                    cur_time = cur_time + recurrence_time;
                end
            end
        end
        
        if i==1
            sim_tmp_1 = ind_tmp_1;
            sim_tmp_2 = ind_tmp_2;
            sim_tmp_3 = ind_tmp_3;
        else
            sim_tmp_1 = [sim_tmp_1;ind_tmp_1];
            sim_tmp_2 = [sim_tmp_2;ind_tmp_2];
            sim_tmp_3 = [sim_tmp_3;ind_tmp_3];
        end
        
    end
    
    sim1{l} = sim_tmp_1;
    sim2{l} = sim_tmp_2;
    sim3{l} = sim_tmp_3;
    
    if mod(l,100)==0
        disp(strcat("l = ",num2str(l)))
    end
    
end

% Save simulated datadata:
filename = strcat('Simulation_1year_',ind_drug_study,' (1000 sim, mod 1 to 3).mat');
% filename = strcat('Simulation_1year_CHQ_PMQ_BPD (1000 sim, mod 1 to 3).mat');
% filename = strcat('Simulation_1year_CHQ_PMQ_VHX (1000 sim, mod 1 to 3).mat');
% filename = strcat('Simulation_1year_DP_PMQ_BPD (1000 sim, mod 1 to 3).mat');
save(filename,'sim1','sim2','sim3')
% Load the simulated data:
% load(filename)
% Note that the data files for the simulated data were too large to be uploaded to a public repository.
% If you want to reproduce the results, the simulations can be repeated with the provided code or you can send an email to estadler@kirby.unsw.edu.au.

% The simulated data were analyzed in R.
