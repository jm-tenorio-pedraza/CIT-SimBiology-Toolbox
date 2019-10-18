function logPrior = createPriorDistribution3(p,PI,H,varargin)
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

% Param indexes of variance params
sigma_indx=H.SigmaParams; % first indx is error variance of tumor volume , second is error variance of immune cell fractions, the rest are individual params variance

% Extracting prior info from PI
lower=([PI.par(:).minValue]);  
upper=([PI.par(:).maxValue]);
mu=([PI.par(:).mu_prior]);
sigma=([PI.par(:).sigma_prior]);

% Defining prior functions of the parameters and indexes
norm_prior=@(x,m,s)sum(log((exp(-(x - m).^2. / (2*s.^2))) ./ sqrt(2 * s.^2 * pi)));
unif_prior=@(x,indx)log((prod(and(x>=lower(indx),x<=upper(indx)))));
%%jeff_prior = @(x)sum(log(1./(exp(x).^2)));
lognorm_prior=@(x,m,s)sum(log(exp(-0.5 * ((log(x) - m)./s).^2) ./ (x .* sqrt(2*pi) .* s)));
invgamma_prior= @(x,a,b)sum(log(b.^a./gamma(a).*x.^(-a-1).*exp(-b./x)));
% Evaluating handle at indexes:
if strcmp(param.type, 'uniform')
     % Uniform priors for fixed params and normal priors for individual
     % params
        logPrior = unif_prior(p(pop_indx), pop_indx)+...
            invgamma_prior((p(sigma_indx)).^2,(mu(sigma_indx)),sigma(sigma_indx))+...
            + sum(arrayfun(@(x) norm_prior(log(p(x.Index))+log(p(x.EtaIndex)), log(p(x.EtaIndex)),...
            p(x.OmegaIndex)), H.IndividualParams))...
            + sum(arrayfun(@(x) norm_prior(log(p(x.Index))+log(p(x.EtaIndex)), log(p(x.EtaIndex)),...
            p(x.OmegaIndex)), H.CellParams));
        
elseif strcmp(param.type,'normal')
        % Normal priors for fixed params and normal priors for individual
        % params
        logPrior = norm_prior((p(pop_indx)), mu(pop_indx),sigma(pop_indx))...
            + norm_prior((p(sigma_indx)),log(mu(sigma_indx)),sigma(sigma_indx))...
            + sum(arrayfun(@(x) lognorm_prior(p(x.Index), log(p(x.EtaIndex)),p(x.OmegaIndex)), H.IndividualParams));
else
    logPrior = lognorm_prior((p(pop_indx)), log(mu(pop_indx)),sigma(pop_indx))...
            + lognorm_prior((p(sigma_indx)),log(mu(sigma_indx)),sigma(sigma_indx))...
            + sum(arrayfun(@(x) lognorm_prior(p(x.Index), log(p(x.EtaIndex)),p(x.OmegaIndex)), H.IndividualParams));
end
    
return


