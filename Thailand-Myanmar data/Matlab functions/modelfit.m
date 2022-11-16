function [U] = modelfit(time,model,par,rel_type)
%%% compute the fraction of uninfected individuals at different time points

%%% Input:
% time      time points for computing the fraction of uninfected
% model     model used for fitting:
%               1: constant relapse rate
%               2: temporal heterogeneity
%               3: population heterogeneity
%               4: temporal & population heterogeneity
% par       vector of parameters, depends on the model:
%           model==1:
%               par(1): washout distribution parameter 1
%               par(2): washout distribution parameter 2
%               par(3): recurrence rate
%           model==2:
%               par(1): washout distribution parameter 1
%               par(2): washout distribution parameter 2
%               par(3): new infection rate
%               par(4): I, initial relapse rate
%               par(5): d, relapse rate decay rate
%           model==3:
%               par(1): washout distribution parameter 1
%               par(2): washout distribution parameter 2
%               par(3): new infection rate
%               par(4): number of recurrence risk groups
%               par(5): relapse risk parameter 1
%               par(6): relapse risk parameter 2
%           model==4:
%               par(1): washout distribution parameter 1
%               par(2): washout distribution parameter 2
%               par(3): new infection rate
%               par(4): number of recurrence risk groups
%               par(5): initial relapse risk parameter 1
%               par(6): initial relapse risk parameter 2
%               par(7): d, relapse rate decay rate
% rel_type  type of the relapse rate distribution in model 3:
%               lognormal: lognormal distribution (default)
%               gamma:     Gamma distribution
%               exp:       Exponential distribution
%%% Output:
% U         vector with the fraction of uninfected individuals at the time
%           points specified by "time"

%%% Function:
if nargin<4 % if the relapse rate distribution is not specified use "lognormal" as the default
    rel_type = "lognormal";
end

opts = odeset('RelTol',1e-6,'AbsTol',1e-9);

t_all = 0:(round(max(time))+1);

if model==1 % constant relapse rate
    if par(1)==0 && par(2)==0 % if there is no drug washout
        [~,S_all] = ode15s(@(t,y) (-par(3))*y,t_all,1,opts);
    else
        [~,S_all] = ode15s(@(t,y) (-par(3))*y+lognpdf(t,par(1),par(2)),t_all,zeros(1,1),opts);
    end
elseif model==2 % temporal heterogeneity
    [~,S_all] = ode15s(@(t,y) (-1).*(par(3)+par(4).*exp(-par(5).*t))*y+lognpdf(t,par(1),par(2)),...
        t_all, zeros(1,1),opts);
elseif model==3 % population heterogeneity
    % relapse rates for the different relapse risk groups:
    if rel_type=="lognormal" % Lognormal distribution of risk of recurrence
        if par(5)==0 && par(6)==0 % if there are no relapses
            r = zeros(par(4),1);
        else
            r = logninv(1/(2*par(4))+(0:(par(4)-1))/par(4),par(5),par(6));
        end
    elseif rel_type=="gamma" % Gamma distribution of risk of recurrence
        if par(5)==0 && par(6)==0 % if there are no relapses
            r = zeros(par(4),1);
        else
            r = gaminv(1/(2*par(4))+(0:(par(4)-1))/par(4),par(5),par(6));
        end
    elseif rel_type=="exp" % Exponential distribution of risk of recurrences
        if par(5)==0
            r = zeros(par(4),1);
        else
            r = expinv(1/(2*par(4))+(0:(par(4)-1))/par(4),par(5));
        end
    end
    rec = r + par(3); % recurrences = relapses + new infections
    m = diag(-rec);
    [~,S_all] = ode15s(@(t,y) m*y+lognpdf(t,par(1),par(2)).*ones(par(4),1)./par(4),t_all,zeros(par(4),1),opts);
elseif model==4 % temporal & population heterogeneity
    % initial relapse rates for the different relapse risk groups:
    if par(5)==0 && par(6)==0 % if there are no relapses
        i = zeros(par(4),1);
    else
        i = logninv(1/(2*par(4))+(0:(par(4)-1))/par(4),par(5),par(6));
    end
    S_all = zeros(length(t_all),par(4));
    for k=1:par(4)
        [~,tmp] = ode15s(@(t,y) (-1).*(par(3)+i(k).*exp(-par(7).*t))*y+...
            lognpdf(t,par(1),par(2))/par(4),t_all, zeros(1,1),opts);
        S_all(:,k) = tmp;
    end
end

S_all(S_all<0) = 0;
if model==1 || model==2
    if par(1)==0 && par(2)==0
        U_all = S_all';
    else
        U_all = (1-logncdf(t_all,par(1),par(2)))+S_all';
    end
elseif model==3 || model==4
    % for the sum over all risk groups:
    %     U_all = (1-logncdf(t_all,par(1),par(2)))+(sum(S_all,2))';
    % for indvidual risk groups:
    U_all = (1-logncdf(t_all,par(1),par(2)))/par(4)+S_all';
end

ind = arrayfun(@(i) find(t_all==time(i)),1:length(time));
U = U_all(:,ind);

end

