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

p.parse(varargin{:});
p=p.Results;


errTol = 1;
curr_p = p0; % Original log-transformed vector of parameters
m = p.m;
gamma = p.gamma;
q_prev = likelihood(curr_p);
i = 1;
random_indx = arrayfun(@(x) x.Index, H.IndividualParams, 'UniformOutput',false);
fixed_indx = [H.PopulationParams H.SigmaParams];
fminunc_options = optimset('Display', 'iter');
while errTol>p.delta
%% Expectation 
% Evaluate integral by MCMC MH sampling

E_likelihood = @(x) likelihood([curr_p([H.PopulationParams]) x curr_p(H.SigmaParams)]);
E_prior = @(x) prior([curr_p([H.PopulationParams]) x curr_p(H.SigmaParams)]);
[randeffects, logP, ~] = mcmc_mh(curr_p([random_indx{:,:}]),E_likelihood, E_prior, m);
q_k = mean(logP); 
q_hat = q_prev + gamma*(q_k - q_prev);
E_Z = mean(randeffects);
%% Maximization 

M_likelihood = @(x) likelihood([x(1:length(H.PopulationParams)) E_Z x(length(H.PopulationParams)+1:end)]);
M_prior = @(x) prior([x(1:length(H.PopulationParams)) E_Z x(length(H.PopulationParams)+1:end)]);

[fixedeffects, logL] = fminunc(@(x)( (M_likelihood(x)+M_prior(x))*(-1)),...
    curr_p([H.PopulationParams H.SigmaParams]),fminunc_options);

% Updating new values
curr_p([random_indx{:,:}]) = E_Z;
curr_p(fixed_indx) = fixedeffects;
errTol = (logL - q_hat);
q_prev = logL;
i = i+1;
gamma = gamma/i;
end
parameters = curr_p;

return

