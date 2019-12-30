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
param.addParameter('type', {'uniform' 'normal' 'inverse gamma'});
param.addParameter('IndividualParams',1,@isnumeric);
param.parse(varargin{:});
param=param.Results;

% Parameter indexes to parameters that vary at the population level
pop_indx=H.PopulationParams;

% Param indexes of variance params
sigma_indx=H.SigmaParams; % first indx is error variance of tumor volume , second is error variance of immune cell fractions, the rest are individual params variance
omega_indx = [H.IndividualParams.OmegaIndex];
psi_indx = [H.CellParams.OmegaIndex];
% Extracting prior info from PI
lower=([PI.par(:).minValue]);  
upper=([PI.par(:).maxValue]);
mu=([PI.par(:).mu_prior]);
sigma=([PI.par(:).sigma_prior]);

% Defining prior functions of the parameters and indexes
norm_prior=@(x,m,s)sum(log((exp(-(x - m).^2. / (2*s.^2))) ./ sqrt(2 * s.^2 * pi)));
unif_prior=@(x,indx)log((prod(and(x>=lower(indx),x<=upper(indx)))));
invgamma_prior= @(x,a,b)sum(log(b.^a./gamma(a).*x.^(-a-1).*exp(-b./x)));
jeffreys_prior= @(x)sum(log(1./x));

% Evaluating handle at indexes:
distributions = [{'uniform/normal/inverse gamma/inverse gamma'};
    {'normal/normal/inverse gamma/inverse gamma'};
    {'uniform/normal/uniform/uniform'};
    {'uniform/normal/jeffreys/jeffreys'};
    {'uniform/normal/inverse gamma/jeffreys'};
  ];

probType = find(ismember(distributions, param.type));
switch probType
    case 1
     % Uniform priors for fixed params and normal priors for individual
     % params
        logPrior = unif_prior(p(pop_indx), pop_indx)+...
            invgamma_prior((p(sigma_indx)).^2,(mu(sigma_indx)),sigma(sigma_indx))+...
            + sum(arrayfun(@(x) norm_prior(log(p(x.Index)), zeros(size((p(x.Index)))),...
            p(x.OmegaIndex)), H.IndividualParams))...
            + sum(arrayfun(@(x) norm_prior(log(p(x.Index)), zeros(size(p(x.Index))),...
            p(x.OmegaIndex)), H.CellParams));
        
    case 2
        % Normal priors for fixed params and normal priors for individual
        % params
        logPrior = norm_prior(log(p(pop_indx)),mu(pop_indx), sigma_prior(pop_indx))+...
            invgamma_prior((p(sigma_indx)).^2,(mu(sigma_indx)),sigma(sigma_indx))+...
            + sum(arrayfun(@(x) norm_prior(log(p(x.Index)), zeros(size((p(x.Index)))),...
            p(x.OmegaIndex)), H.IndividualParams))...
            + sum(arrayfun(@(x) norm_prior(log(p(x.Index)), zeros(size(p(x.Index))),...
            p(x.OmegaIndex)), H.CellParams));
    case 3
    logPrior = unif_prior(p(pop_indx), pop_indx)+...
            unif_prior((p(sigma_indx)),sigma_indx)+...
            + sum(arrayfun(@(x) norm_prior(log(p(x.Index)), zeros(size((p(x.Index)))),...
            p(x.OmegaIndex)), H.IndividualParams))...
            + sum(arrayfun(@(x) norm_prior(log(p(x.Index)), zeros(size(p(x.Index))),...
            p(x.OmegaIndex)), H.CellParams));
    case 4
         logPrior = unif_prior(p(pop_indx), pop_indx)+...
            jeffreys_prior((p(sigma_indx)).^2)+...
            + sum(arrayfun(@(x) norm_prior(log(p(x.Index)), zeros(size((p(x.Index)))),...
            p(x.OmegaIndex)), H.IndividualParams))...
            + sum(arrayfun(@(x) norm_prior(log(p(x.Index)), zeros(size(p(x.Index))),...
            p(x.OmegaIndex)), H.CellParams));
    case 5
        logPrior = unif_prior(p(pop_indx), pop_indx)+...
            jeffreys_prior((p(setdiff(sigma_indx, [omega_indx, psi_indx]))).^2)+...
            invgamma_prior((p([omega_indx, psi_indx])).^2,(mu([omega_indx, psi_indx])),sigma([omega_indx, psi_indx]))+...
            + sum(arrayfun(@(x) norm_prior(log(p(x.Index)), zeros(size((p(x.Index)))),...
            p(x.OmegaIndex)), H.IndividualParams))...
            + sum(arrayfun(@(x) norm_prior(log(p(x.Index)), zeros(size(p(x.Index))),...
            p(x.OmegaIndex)), H.CellParams));
  
        
end
    
return


