function [params,logP, accept]=mcmc_mh(p0,logL_fun,prior_fun,n_samples,varargin)
if nargin<3
    error('GWMCMC:toofewinputs','AMCMC requires atleast 3 inputs.')
end
p=inputParser;
p.addParameter('StepSize',2.38,@isnumeric);
p.addParameter('Thinning',1,@isnumeric);
p.addParameter('BurnIn',0,@isnumeric);
p.addParameter('ProposalSigma',diag(ones(size(p0))),@isnumeric);

p.parse(varargin{:});
p=p.Results;

partition=floor(n_samples*p.Thinning/20);
n_params=length(p0);

% Initial values
curr_p=p0;
curr_L=(logL_fun(curr_p)+prior_fun(curr_p));

prop_sigma=p.ProposalSigma*p.StepSize^2/(length(p0));


% Allocating space
params=NaN(p.BurnIn+n_samples*p.Thinning,n_params);
U=log(rand(p.BurnIn+n_samples*p.Thinning,1));

% Generating multivariate normal proposals
proposals=mvnrnd(zeros(1,n_params),prop_sigma,p.BurnIn+n_samples*p.Thinning);
logP=NaN(n_params,1);
burnIn = p.BurnIn*n_samples;

accept=0;
%% Burn In step
datestr(clock)
disp('Starting burn-in phase')
    for i=1:burnIn
        % Propose new sample
        prop_p=curr_p+proposals(i,:);
        % Evaluate with prior
        prior_L=prior_fun(prop_p);
        if isinf(prior_L) || ~isreal(prior_L) || isnan(prior_L)
            params(i,:)=curr_p;
            continue
        else
            % Evaluate numerator of MH ratio
            prop_L=logL_fun(prop_p)+prior_L;
            r=prop_L-(curr_L);
            if U(i)<=min(r,0)
                curr_L=prop_L;
                curr_p=prop_p;
                accept=accept+1;
            else
            end
            params(i,:)=curr_p;
             if mod(i,p.burnIn/20)==0
                perc=i/(n_samples*p.Thinning);
                fprintf('%.1f percent completed',perc*100)
                acceptInter = accept/(i)*100;
                fprintf('%.1f percent accepted',acceptInter);
                display(datestr(clock))
             end
        end

    end
%% Adaptive step
disp('Starting adaptive phase')
for i=1:n_samples
    for j=1:p.Thinning
        indx=(i-1)*p.Thinning+j;
        prop_p=curr_p+proposals(indx,:);
        prior_L=prior_fun(prop_p);
        if isinf(prior_L) || ~isreal(prior_L) || isnan(prior_L)
            continue
        else

            prop_L=logL_fun(prop_p)+prior_L;
            if (prop_L)<=-1e7
                continue;
            else
                r=prop_L-(curr_L);
                if U(indx)<=min(r,0)
                    curr_L=prop_L;
                    curr_p=prop_p;
                    accept=accept+1;
                end
            end
        end
    end
    params(i,:)=curr_p;
    logP(i)=curr_L;
    if mod(i,partition)==0      
        perc=(p.BurnIn+i*p.Thinning)/(n_samples*p.Thinning);
        fprintf('%.1f percent completed\n',perc*100)
        acceptInter = accept/(i*p.Thinning+p.BurnIn)*100;
        if acceptInter<20
            proposals=proposals*0.1;
        elseif acceptInter>50
            proposals = proposals*10;
        end
        fprintf('%.1f percent accepted\n',acceptInter);
        display(datestr(clock))

    end
end
accept = accept/(n_samples*p.Thinning)*100; 
end