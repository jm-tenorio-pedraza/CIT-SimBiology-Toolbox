function [params,logP, accept]=adaptive_mcmc(p0,logL_fun,prior_fun,n_samples,varargin)
if nargin<3
    error('GWMCMC:toofewinputs','AMCMC requires atleast 3 inputs.')
end
p=inputParser;
p.addParameter('StepSize',.2,@isnumeric);
p.addParameter('Thinning',1,@isnumeric);
p.addParameter('AdaptSteps',100,@isnumeric);
p.addParameter('BurnIn',1,@isnumeric);
p.addParameter('PreviousSamples',p0,@isnumeric);
p.addParameter('SigmaThinning',10,@isnumeric);

p.parse(varargin{:});
p=p.Results;

% intermed = floor((n_samples*p.Thinning+p.BurnIn)/20); %estimate number for 
% displaying progression and intermediate acceptance rate every 5%
partition=floor(n_samples/p.AdaptSteps);
%% Proposal step
n_col=length(p0);

initial_sigma=eye(n_col)*((p.StepSize/n_col));
curr_p=p0;
curr_L=(logL_fun(curr_p)+prior_fun(curr_p));
curr_Prior=prior_fun(curr_p);
prop_sigma=(initial_sigma);
U=log(rand(p.BurnIn+n_samples*p.Thinning,1));

accept=0;
if size(p.PreviousSamples,1)==1
proposal_samples=randn(p.BurnIn,n_col)*prop_sigma';

params=NaN(n_samples*p.Thinning+p.BurnIn, n_col);
logP=NaN(n_samples*p.Thinning+p.BurnIn,1);

%draw steps for proposals outside of the loop
datestr(clock)
disp('Starting burn-in phase')
    for i=1:p.BurnIn
    % Get new proposal
        prop_p=curr_p+proposal_samples(i,:);
        % Calculate prior prob
        prop_Prior=prior_fun(prop_p);
        % First criteria for rejection
        if isinf(prop_Prior) || ~isreal(prop_Prior) || isnan(prop_Prior)
            params(i,:)=curr_p;
            logP(i)=curr_L;    
            continue        
        else
            % Second criteria for rejection
%             r_prior=prop_Prior-curr_Prior+sum(prop_p)-sum(curr_p);
%             if U(i)> min(0,r_prior)
%                 params(i,:)=curr_p;
%                 logP(i)=curr_L;
%                 continue
%             else
                % Third criteria for rejection
                prop_L=logL_fun(prop_p)+prop_Prior;
                r=prop_L-curr_L+sum(prop_p)-sum(curr_p);
                if U(i)<=min(r,0)
                    curr_L=prop_L;
                    curr_p=prop_p;
                    curr_Prior=prop_Prior;
                    accept=accept+1;
                else
                end
                params(i,:)=curr_p;
                logP(i)=curr_L;
%             end
                if mod(i,p.BurnIn/20)==0
                    perc=i/(p.BurnIn+n_samples*p.Thinning);
                    fprintf('%.1f percent completed\n',perc*100)
                    acceptInter = accept/(i)*100;
                    fprintf('%.1f percent accepted\n',acceptInter)
                    display(datestr(clock))
                end
  
        end
    end
    p_chol=(chol(cov(params(1:p.SigmaThinning:p.BurnIn,:),'omitrows')));
else
    p_chol=(chol(cov(p.PreviousSamples)));
end
%% Adaptive step
posterior_samples=randn(n_samples*p.Thinning,n_col);
U=log(rand(n_samples*p.Thinning,1));

disp('Starting adaptive phase')
for i=1:n_samples
    for j=1:p.Thinning
        % Calculate index for storage
        indx=(i-1)*p.Thinning+j;
        % Calculate weight of sigma matrix according to iteration number
        % along the sample step
        scale_weight=p.BurnIn+(size(p.PreviousSamples,1))*p.Thinning;
        % Calculate weighted covariance matrix decomposition
        prop_sigma=((p.BurnIn/(scale_weight+indx)*initial_sigma+...
            (scale_weight-p.BurnIn+indx)/(scale_weight+indx)*p_chol*p.StepSize));
        % Calculate proposal vector
        prop_p=curr_p+posterior_samples(indx,:)*prop_sigma';
        % Calculate prior prob of proposal
        prop_Prior=prior_fun(prop_p);
        % First rejection criteria
        if isinf(prop_Prior) || ~isreal(prop_Prior) || isnan(prop_Prior)
            params(indx+p.BurnIn,:)=curr_p;
            logP(indx+p.BurnIn)=curr_L;
            continue
        else
            % Second rejection criteria
%             r_prior=prop_Prior-curr_Prior+sum(prop_p)-sum(curr_p);
%             if U(i)> min(0,r_prior)
%                 params(indx+p.BurnIn,:)=curr_p;
%                 logP(indx+p.BurnIn)=curr_L;
%                 continue
%             else
                % Third rejection criteria
                prop_L=logL_fun(prop_p)+prop_Prior;
                r=prop_L-curr_L+ sum(prop_p)-sum(curr_p);
                if U(indx)<=min(r,0)
                    curr_L=prop_L;
                    curr_p=prop_p;
%                     curr_Prior=prop_Prior;
                    accept=accept+1;
                end
%             end
        params(indx+p.BurnIn,:)=curr_p;
        logP(indx+p.BurnIn)=curr_L;
        end
    if mod(i,partition)==0
        try
        p_chol=(chol(cov([p.PreviousSamples;params(1:p.SigmaThinning:indx+p.BurnIn,:)], 'omitrows')));
        catch
        lastwarn('')
        end
        perc=(indx+p.BurnIn)/(n_samples*p.Thinning+p.BurnIn);
        fprintf('%.1f percent completed\n',perc*100)
        acceptInter = accept/(indx+p.BurnIn)*100;
        fprintf('%.1f percent accepted\n',acceptInter);
        display(datestr(clock))
    end
    end
end
accept = accept/(n_samples*p.Thinning+p.BurnIn)*100; 

return