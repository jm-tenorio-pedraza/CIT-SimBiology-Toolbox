function logPrior = createPriorDistribution2(p,PI,H,varargin)
% Function to create a prior handle using information available in PI.par
% Inputs:
%           - PI: structure with a 'par' field with 'startValue',
%           'minValue', 'maxValue' fields with entries for each parameter
%           - H: structure with 'PopulationParams', "IndividualParams', and
%           'SigmaParams' fields. Each field contains a double vector with
%           indexes that correspond to the position in the parameter vector
% Outputs:
%           - prior_handle: handle to evaluate a multivariate function of
%           log-transformed parameters

if nargin<2
    error('createPriorDistribution requires at least 2 inputs')
end
if size(p,1)>size(p,2)
    p=p';
end
param=inputParser;
param.addParameter('type', 'uniform');
param.addParameter('IndividualParams',1,@isnumeric);
param.parse(varargin{:});
param=param.Results;

% Parameter indexes to parameters that vary at the population level
pop_indx=H.PopulationParams;

% Param indexes of individual params
ind_indx=[H.IndividualParams(:).Index];

% Param indexes of variance params
sigma_indx=H.SigmaParams; % first indx is error variance of tumor volume , second is error variance of immune cell fractions, the rest are individual params variance

% Extracting prior info from PI
lower=log([PI.par(:).minValue]);  
upper=log([PI.par(:).maxValue]);
mu=log([PI.par(:).mu_prior]);
sigma=([PI.par(:).sigma_prior]);

% Defining prior functions of the parameters and indexes
norm_prior=@(x,m,s)sum(log(exp(-(x-m).^2./(2*s.^2))./sqrt(2*s.^2*pi)));
unif_prior=@(x,indx)log((prod(and(x>=lower(indx),x<=upper(indx)))));
jeff_prior = @(x)sum(log(1./(exp(x).^2)));
lognorm_prior=@(x,m,s)sum(log(exp(-(log(x)-m).^2./(2*s.^2))./sqrt(x.^2*2.*s.^2*pi)));

% Evaluating handle at indexes:
if strcmp(param.type, 'uniform')
     % Uniform priors for fixed params and normal priors for individual
     % params
        logPrior=(unif_prior(p(pop_indx), pop_indx)+jeff_prior(p(sigma_indx))...
        +norm_prior(p(ind_indx), mu(ind_indx), sigma(ind_indx)));
    
elseif strcmp(param.type,'normal')
        % Normal priors for fixed params and normal priors for individual
        % params
        logPrior=(norm_prior((p(pop_indx)), mu(pop_indx),sigma(pop_indx))...
            + norm_prior((p(sigma_indx)),mu(sigma_indx),sigma(sigma_indx))...
            + norm_prior((p(ind_indx)), mu(ind_indx),sigma(ind_indx)));
else
    logPrior=(lognorm_prior((p(pop_indx)), mu(pop_indx),sigma(pop_indx))...
            + lognorm_prior((p(sigma_indx)),mu(sigma_indx),sigma(sigma_indx))...
            + lognorm_prior((p(ind_indx)), mu(ind_indx),sigma(ind_indx)));
end
    
return


