function [nllh] = modelfit_1rec(data,model,followup,par,reinf)
%%% fit a model using an ode solver (ode15s)
%%% fit to data for one recurrence, data grouped by drug and study, no
%%% relapses for PMQ treated patients


%%% Input:
% data      structure containing the data with:
%               data.id
%               data.time1
%               data.event1
%               data.time2
%               data.event2
%               data.drug
%               data.study
%               data.drug_study
%               data.followup
% model     model used for fitting:
%               1: constant relapse rate
%               2: temporal heterogeneity
%               3: population heterogeneity
%               4: temporal & population heterogeneity
% followup  follow-up scheme:
%               1: daily follow-up
%               2: weekly follow-up
%               3: fortnightly follow-up
%               4: 4-weekly follow-up
%               5: follow-up at the beginning of weeks 2,4,8,12,...
%               6: follow-up in the middle of weeks 2,4,8,12,...
%               7: follow-up at the end of weeks 2,4,8,12,...
% par       vector of parameters, depends on the model:
%           model==1:
%               par(1): washout distribution parameter 1, AS group
%               par(2): washout distribution parameter 2, AS group
%               par(3): washout distribution parameter 1, CHQ group
%               par(4): washout distribution parameter 2, CHQ group
%               par(5): washout distribution parameter 1, DP group
%               par(6): washout distribution parameter 2, DP group
%               par(7): recurrence rate, AS/CHQ & VHX group
%               par(8): recurrence rate, PMQ & VHX group
%               par(9): recurrence rate, PMQ & BPD group
%                   if reinf==1: par(9)=par(8)
%           model==2:
%               par(1): washout distribution parameter 1, AS group
%               par(2): washout distribution parameter 2, AS group
%               par(3): washout distribution parameter 1, CHQ group
%               par(4): washout distribution parameter 2, CHQ group
%               par(5): washout distribution parameter 1, DP group
%               par(6): washout distribution parameter 2, DP group
%               par(7): new infection rate, VHX study
%               par(8): new infection rate, BPD study
%                   if reinf==1: par(7)=par(8) & decrease all following
%                   parameter numbers by one
%               par(9): I, initial relapse rate, AS/CHQ group
%               par(10): d, relapse rate decay rate
%           model==3:
%               par(1): washout distribution parameter 1, AS group
%               par(2): washout distribution parameter 2, AS group
%               par(3): washout distribution parameter 1, CHQ group
%               par(4): washout distribution parameter 2, CHQ group
%               par(5): washout distribution parameter 1, DP group
%               par(6): washout distribution parameter 2, DP group
%               par(7): new infection rate, VHX study
%               par(8): new infection rate, BPD study
%                   if reinf==1: par(7)=par(8) & decrease all following
%                   parameter numbers by one
%               par(9): relapse risk parameter 1, AS/CHQ group
%               par(10): relapse risk parameter 2, AS/CHQ group
%               par(11): n, number of relapse risk groups -> integer
%           model==4:
%               par(1): washout distribution parameter 1, AS group
%               par(2): washout distribution parameter 2, AS group
%               par(3): washout distribution parameter 1, CHQ group
%               par(4): washout distribution parameter 2, CHQ group
%               par(5): washout distribution parameter 1, DP group
%               par(6): washout distribution parameter 2, DP group
%               par(7): new infection rate, VHX study
%               par(8): new infection rate, BPD study
%                   if reinf==1: par(7)=par(8) & decrease all following
%                   parameter numbers by one
%               par(9): initial relapse risk parameter 1, AS/CHQ group
%               par(10): initial relapse risk parameter 2, AS/CHQ group
%               par(11): d, relapse rate decay rate
%               par(12): n, number of relapse risk groups -> integer
% reinf     number of new infection rates to use in the model fit
%               1: same new infection rate for both studies
%               2: different new infection rate for the two different studies
%%% Output:
% nllh      negative loglikelihood for the model


%%% Uses the following functions:
% modelfit.m    to compute the fraction of uninfected individuals


%%% Function:
drug_study_groups = ["AS_VHX","CHQ_VHX","CHQ/PMQ_VHX","CHQ/PMQ_BPD","DP/PMQ_BPD"];
nllh_tmp = zeros(5,1); % negative loglikelihood for each drug & study group

for a = 1:5 % for the antimalarial treatments and study groups
    data_tmp = [data.time1(data.drug_study==drug_study_groups(a)),data.event1(data.drug_study==drug_study_groups(a))];
    
    t_max = round(max(data.time1))+1;
    t_all = 0:t_max; % vector for all times from 0 to the maximal time
    times0 = data_tmp(data_tmp(:,2)==0,1); % time of censoring for individuals with 0 recurrences
    times1 = data_tmp(data_tmp(:,2)==1,1); % recurrence time for individuals with at least 1 recurrence
    
    % different follow-up schemes (last sample):
    if followup == 1 % daily sampling
        ls = 0:(max(t_max)-1);
    elseif followup == 2 % weekly sampling
        ls = max((1:max(t_max))-7,0);
    elseif followup == 3 % fornightly sampling
        ls = max((1:max(t_max))-14,0);
    elseif followup == 4 % four-weekly sampling
        ls = max((1:max(t_max))-28,0);
    elseif followup == 5 % sampling at the beginning of the week
        ls = [zeros(1,14),repmat(8,1,14),(ceil((29:max(t_max))/28)-1)*28-6];
    elseif followup == 6 % sampling in the middle of the week
        ls = [zeros(1,14),repmat(11,1,14),(ceil((29:max(t_max))/28)-1)*28-3];
    elseif followup == 7 % sampling at the end of the week
        ls = [zeros(1,14),repmat(14,1,14),(ceil((29:max(t_max))/28)-1)*28];
    else
        ls = [];
    end
    
    if model==1 % constant relapse rate
        if a==1
            U_all = modelfit(t_all,model,[par(1),par(2),par(7)]);
        elseif a==2
            U_all = modelfit(t_all,model,[par(3),par(4),par(7)]);
        elseif a==3
            U_all = modelfit(t_all,model,[par(3),par(4),par(8)]);
        elseif a==4
            if reinf==1
                U_all = modelfit(t_all,model,[par(3),par(4),par(8)]);
            elseif reinf==2
                U_all = modelfit(t_all,model,[par(3),par(4),par(9)]);
            end
        elseif a==5
            if reinf==1
                U_all = modelfit(t_all,model,[par(5),par(6),par(8)]);
            elseif reinf==2
                U_all = modelfit(t_all,model,[par(5),par(6),par(9)]);
            end
        end
        
    elseif model==2 % temporal heterogeneity
        if a==1
            if reinf==1
                U_all = modelfit(t_all,model,[par(1),par(2),par(7),par(8),par(9)]);
            elseif reinf==2
                U_all = modelfit(t_all,model,[par(1),par(2),par(7),par(9),par(10)]);
            end
        elseif a==2
            if reinf==1
                U_all = modelfit(t_all,model,[par(3),par(4),par(7),par(8),par(9)]);
            elseif reinf==2
                U_all = modelfit(t_all,model,[par(3),par(4),par(7),par(9),par(10)]);
            end
        elseif a==3
            if reinf==1
                U_all = modelfit(t_all,model,[par(3),par(4),par(7),0,par(9)]);
            elseif reinf==2
                U_all = modelfit(t_all,model,[par(3),par(4),par(7),0,par(10)]);
            end
        elseif a==4
            if reinf==1
                U_all = modelfit(t_all,model,[par(3),par(4),par(7),0,par(9)]);
            elseif reinf==2
                U_all = modelfit(t_all,model,[par(3),par(4),par(8),0,par(10)]);
            end
        elseif a==5
            if reinf==1
                U_all = modelfit(t_all,model,[par(5),par(6),par(7),0,par(9)]);
            elseif reinf==2
                U_all = modelfit(t_all,model,[par(5),par(6),par(8),0,par(10)]);
            end
        end
        
    elseif model==3 % population heterogeneity
        if a==1
            if reinf==1
                U_all = modelfit(t_all,model,[par(1),par(2),par(7),par(10),par(8),par(9)]);
            elseif reinf==2
                U_all = modelfit(t_all,model,[par(1),par(2),par(7),par(11),par(9),par(10)]);
            end
        elseif a==2
            if reinf==1
                U_all = modelfit(t_all,model,[par(3),par(4),par(7),par(10),par(8),par(9)]);
            elseif reinf==2
                U_all = modelfit(t_all,model,[par(3),par(4),par(7),par(11),par(9),par(10)]);
            end
        elseif a==3
            if reinf==1
                U_all = modelfit(t_all,model,[par(3),par(4),par(7),par(10),0,0]);
            elseif reinf==2
                U_all = modelfit(t_all,model,[par(3),par(4),par(7),par(11),0,0]);
            end
        elseif a==4
            if reinf==1
                U_all = modelfit(t_all,model,[par(3),par(4),par(7),par(10),0,0]);
            elseif reinf==2
                U_all = modelfit(t_all,model,[par(3),par(4),par(8),par(11),0,0]);
            end
        elseif a==5
            if reinf==1
                U_all = modelfit(t_all,model,[par(5),par(6),par(7),par(10),0,0]);
            elseif reinf==2
                U_all = modelfit(t_all,model,[par(5),par(6),par(8),par(11),0,0]);
            end
        end
        U_all = sum(U_all,1);
        
        elseif model==4 % temporal & population heterogeneity
        if a==1
            if reinf==1
                U_all = modelfit(t_all,model,[par(1),par(2),par(7),par(11),par(8),par(9),par(10)]);
            elseif reinf==2
                U_all = modelfit(t_all,model,[par(1),par(2),par(7),par(12),par(9),par(10),par(11)]);
            end
        elseif a==2
            if reinf==1
                U_all = modelfit(t_all,model,[par(3),par(4),par(7),par(11),par(8),par(9),par(10)]);
            elseif reinf==2
                U_all = modelfit(t_all,model,[par(3),par(4),par(7),par(12),par(9),par(10),par(11)]);
            end
        elseif a==3
            if reinf==1
                U_all = modelfit(t_all,model,[par(3),par(4),par(7),par(11),0,0,0]);
            elseif reinf==2
                U_all = modelfit(t_all,model,[par(3),par(4),par(7),par(12),0,0,0]);
            end
        elseif a==4
            if reinf==1
                U_all = modelfit(t_all,model,[par(3),par(4),par(7),par(11),0,0,0]);
            elseif reinf==2
                U_all = modelfit(t_all,model,[par(3),par(4),par(8),par(12),0,0,0]);
            end
        elseif a==5
            if reinf==1
                U_all = modelfit(t_all,model,[par(5),par(6),par(7),par(11),0,0,0]);
            elseif reinf==2
                U_all = modelfit(t_all,model,[par(5),par(6),par(8),par(12),0,0,0]);
            end
        end
        U_all = sum(U_all,1);
        
    end
    
    ind_0 = arrayfun(@(i) find(t_all==times0(i)),1:length(times0)); % individuals with 0 recurrences, time of censoring
    U_0 = U_all(:,ind_0);
    ind_1 = arrayfun(@(i) find(t_all==times1(i)),1:length(times1)); % individuals with at least 1 recurrence, time of 1st recurrence
    U_1 = U_all(:,ind_1);
    ind_1_d = arrayfun(@(i) find(t_all==ls(times1(i))),1:length(times1)); % individuals with 1 recurrence, time of last visit before 1st recurrence
    U_1_d = U_all(:,ind_1_d);
    
    G_1 = U_1_d-U_1;
    if any(any(G_1<=0))
        G_1(G_1<=0) = 1e-10;
    end
    if any(any(U_0<=0))
        U_0(U_0<=0) = 1e-10;
    end
    
    nllh_tmp(a,1) =(-1)*(sum(log(U_0))+sum(log(G_1)));
    
end

% The negative loglikelihood is the sum of the negative loglikelihoods for
% all antimalarial treatments and studies:
nllh = sum(nllh_tmp);

end

