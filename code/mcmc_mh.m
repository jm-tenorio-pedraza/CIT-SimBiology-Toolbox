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
burnIn = p.BurnIn*n_samples;


% Allocating space
params=NaN(burnIn+(n_samples-burnIn)*p.Thinning,n_params);
U=log(rand(burnIn+(n_samples-burnIn)*p.Thinning,1));

% Generating multivariate normal proposals
proposals=mvnrnd(zeros(1,n_params),prop_sigma,size(params,1));
logP=NaN(size(params,1),1);

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
        logP(i) = curr_L;
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
        logP(i) = curr_L;
    end
    if mod(i,burnIn/20)==0
        perc=i/(burnIn+(n_samples-burnIn)*p.Thinning);
        fprintf('%.1f percent completed\n',perc*100)
        acceptInter = accept/(i)*100;
        if acceptInter<20
            proposals=proposals*0.1;
        elseif acceptInter>40
            proposals = proposals*10;
        end
        fprintf('%.1f percent accepted\n',acceptInter);
        display(datestr(clock))
    end
    
end
%% Adaptive step
disp('Starting adaptive phase')
for i=1:(n_samples-burnIn)
    for j=1:p.Thinning
        indx=(i-1)*p.Thinning+j+burnIn;
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
    params(indx,:)=curr_p;
    logP(indx)=curr_L;
    if mod(i,partition)==0      
        perc=(burnIn+i*p.Thinning)/(burnIn+(n_samples-burnIn)*p.Thinning);
        fprintf('%.1f percent completed\n',perc*100)
        acceptInter = accept/(indx)*100;
        fprintf('%.1f percent accepted\n',acceptInter);
        display(datestr(clock))

    end
end
accept = accept/((n_samples-burnIn)*p.Thinning+burnIn)*100; 
end