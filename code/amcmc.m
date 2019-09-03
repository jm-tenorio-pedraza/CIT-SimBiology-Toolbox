function [params,logP, accept,proposal_sigma]=amcmc(p0,logL_fun,prior_fun,n_samples,varargin)
if nargin<3
    error('GWMCMC:toofewinputs','AMCMC requires atleast 3 inputs.')
end
p=inputParser;
p.addParameter('StepSize',.2,@isnumeric);
p.addParameter('Thinning',10,@isnumeric);
p.addParameter('AdaptSteps',1000,@isnumeric);
p.addParameter('BurnIn',1e5,@isnumeric);
p.addParameter('PreviousSamples',p0,@isnumeric);
p.parse(varargin{:});
p=p.Results;

% intermed = floor((n_samples*p.Thinning+p.BurnIn)/20); %estimate number for 
% displaying progression and intermediate acceptance rate every 5%
partition=floor(n_samples/p.AdaptSteps);
%% Proposal step
n_col=length(p0);

proposal_sigma=NaN(p.AdaptSteps+2,n_col);


curr_sigma=repelem((p.StepSize/n_col).^2, n_col);
curr_p=p0;
prop_sigma=curr_sigma;
curr_L=(logL_fun(curr_p)+prior_fun(curr_p));
U=log(rand(p.BurnIn,1));
proposal_sigma(1,:)=prop_sigma;

accept=0;
if size(p.PreviousSamples,1)==1
proposal_samples=randn(p.BurnIn,n_col);

params=NaN(size(proposal_samples));
%draw steps for proposals outside of the loop
datestr(clock)
display('Starting burn-in phase')
    for i=1:p.BurnIn

        prop_p=curr_p+proposal_samples(i,:).*sqrt(prop_sigma);
        prior_L=prior_fun(prop_p);
        if isinf(prior_L) || ~isreal(prior_L) || isnan(prior_L)
            params(i,:)=curr_p;
            continue
        else
            prop_L=logL_fun(prop_p)+prior_L;
            r=prop_L-(curr_L);
            if U(i)<=min(r,0)
                curr_L=prop_L;
                curr_p=prop_p;
                accept=accept+1;
            else
            end
            params(i,:)=curr_p;
             if mod(i,p.BurnIn/20)==0
                perc=i/(p.BurnIn+n_samples*p.Thinning);
                display(sprintf('%.1f percent completed',perc*100))
                acceptInter = accept/(i)*100;
                display(sprintf('%.1f percent accepted',acceptInter));
                display(datestr(clock))
             end
        end

    end

    p_chol=diag(chol(cov(params,'omitrows')));
else
    p_chol=diag(chol(cov(p.PreviousSamples)));
end
proposal_sigma(2,:)=p_chol;
%% Adaptive step
posterior_samples=randn(n_samples*p.Thinning,n_col);
params=NaN(n_samples, n_col);
U=log(rand(n_samples*p.Thinning,1));
logP=NaN(n_samples,1);

display('Starting adaptive phase')
for i=1:n_samples
    for j=1:p.Thinning
        indx=(i-1)*p.Thinning+j;
        scale_weight=p.BurnIn+(size(p.PreviousSamples,1))*p.Thinning;
        prop_sigma=((p.BurnIn/(scale_weight+indx)*p.StepSize/n_col+(scale_weight-p.BurnIn+indx)/(scale_weight+indx)*p_chol'*p.StepSize)).^2;
        prop_p=curr_p+posterior_samples(indx,:).*sqrt(prop_sigma);
        prior_L=prior_fun(prop_p);
        if isinf(prior_L) || ~isreal(prior_L) || isnan(prior_L)
            continue
        else

            prop_L=logL_fun(prop_p)+prior_L;
            if (prop_L)<=-1e3
                continue;
            else
                r=prop_L+sum(log(sqrt(prop_sigma)))-(curr_L+sum(log(sqrt(curr_sigma))));
                if U(indx)<=min(r,0)
                    curr_L=prop_L;
                    curr_p=prop_p;
                    curr_sigma=prop_sigma;
                    accept=accept+1;
                end
            end
        end
    end
    params(i,:)=curr_p;
    logP(i)=curr_L;
    if mod(i,partition)==0
        try
        p_chol=diag(chol(cov([p.PreviousSamples;params(1:i,:)], 'omitrows')));
        catch
        lastwarn('')
        end
        perc=(p.BurnIn+i*p.Thinning)/(p.BurnIn+n_samples*p.Thinning);
        display(sprintf('%.1f percent completed',perc*100))
        acceptInter = accept/(i*p.Thinning)*100;
        display(sprintf('%.1f percent accepted',acceptInter));
        display(datestr(clock))
        proposal_sigma(i/partition+2,:)=p_chol;

    end
end
accept = accept/(n_samples*p.Thinning+p.BurnIn)*100; 
end