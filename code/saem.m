function [parameters, logL] = saem(p0, likelihood,prior, H,varargin)

% Inputs:
%           % p0: Start value of all parameters
%           likelihood: likelihood function to optimize
%           H: Hierarchical structure of the parameters

% Outputs:
%           parameters: vector of optimized parameters
%           logL: logLikelihood evaluated at the optimum

p = inputParser;
p.addParameter('Prior', 'Normal')
p.addParameter('delta', 1e-6)
p.addParameter('m', 2e3)
p.addParameter('gamma', 1)
p.addParameter('StepSize', 2.38^2/length(p0))
p.addParameter('MinFunc', 1)
p.addParameter('OutputFn', [])
p.addParameter('BurnIn',500)

p.parse(varargin{:});
p=p.Results;


errTol = 1;
curr_p = p0; % Original log-transformed vector of parameters
m = p.m;
gamma = p.gamma;
q_prev = prior(curr_p)+likelihood(curr_p);
i = 1;
min_fx = p.MinFunc;
random_indx = arrayfun(@(x) x.Index, H.IndividualParams, 'UniformOutput',false);
fixed_indx = [H.PopulationParams];
fminunc_options = optimset('Display', 'iter','MaxIter',2e3,'MaxFunEvals',2e3);
while errTol>p.delta
    %% Expectation
    % Evaluate integral by MCMC MH sampling
    
    E_likelihood = @(x) likelihood([curr_p([H.PopulationParams]) x curr_p(H.SigmaParams)]);
    E_prior = @(x) prior([curr_p([H.PopulationParams]) x curr_p(H.SigmaParams)]);
    [randeffects, logP, ~] = mcmc_mh(curr_p([random_indx{:,:}]),E_likelihood, E_prior, m,'StepSize',p.StepSize);
    q_k = mean(logP);
    q_hat = q_prev + gamma*(q_k - q_prev);
    E_Z = mean(randeffects(p.BurnIn:end,:));
    curr_p([random_indx{:,:}]) = E_Z;
    if ~isempty(H.SigmaParams)
        E_prior = @(x) prior([curr_p([H.PopulationParams]) E_Z x]);
        PI = p.OutputFn(curr_p);
        sigma_indx = size(H.IndividualParams.OmegaIndex,1)+1:length(H.SigmaParams);
        E_likelihood = (@(x)sum(getErrors(PI,exp(x(sigma_indx))))*(-1));
        [sigma, logP_sigma, ~] = mcmc_mh(curr_p(H.SigmaParams),...
            E_likelihood, E_prior, m);
        E_sigma = mean(sigma(p.BurnIn:end,:));
        q_sigma = mean(logP_sigma(p.BurnIn:end));
         
         curr_p(H.SigmaParams) = E_sigma;
         q_hat = (q_hat +q_sigma)/2;
    end
    
   
    %% Maximization
     if ~isempty(H.SigmaParams)
        M_likelihood = @(x) likelihood([x(1:length(H.PopulationParams)) E_Z E_sigma]);
        M_prior = @(x) prior([x(1:length(H.PopulationParams)) E_Z E_sigma]);
     else
        M_likelihood = @(x) likelihood([x(1:length(H.PopulationParams)) E_Z]);
        M_prior = @(x) prior([x(1:length(H.PopulationParams)) E_Z]);
     end
     
    if strcmp(min_fx, 'fminunc')
        [fixedeffects, logL] = fminunc(@(x)( (M_likelihood(x)+M_prior(x))*(-1)),...
            curr_p([H.PopulationParams]),fminunc_options);
    elseif strcmp(min_fx, 'fminsearch')
        [fixedeffects, logL] = fminsearch(@(x)( (M_likelihood(x)+M_prior(x))*(-1)),...
            curr_p([H.PopulationParams]),fminunc_options);
    end
    
    % Updating new values
    curr_p([random_indx{:,:}]) = E_Z;
    curr_p(fixed_indx) = fixedeffects;
    errTol = (-logL - q_hat);
    q_prev = -logL;
    i = i+1;
    gamma = gamma/i;
end
parameters = curr_p;

return

