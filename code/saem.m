function [parameters, logL] = saem(p0, residuals,prior, H,PI,varargin)

% Inputs:
%           % p0: Start value of all parameters
%           likelihood: likelihood function to optimize
%           H: Hierarchical structure of the parameters

% Outputs:
%           parameters: vector of optimized parameters
%           logL: logLikelihood evaluated at the optimum

p = inputParser;
p.addParameter('Prior', 'Normal')
p.addParameter('delta', 1e-4)
p.addParameter('m', 2e3)
p.addParameter('gamma', 1)
p.addParameter('StepSize', 2.38^2/length(p0))
p.addParameter('MinFunc', 1)
p.addParameter('OutputFn', [])
p.addParameter('SimFn', [])
p.addParameter('BurnIn',500)

p.parse(varargin{:});
p=p.Results;


errTol = 1;
curr_p = p0; % Original log-transformed vector of parameters
m = p.m;
gamma = p.gamma;
i = 1;
min_fx = p.MinFunc;
random_indx = arrayfun(@(x) x.Index, H.IndividualParams, 'UniformOutput',false);
fixed_indx = [H.PopulationParams];
fminunc_options = optimset('Display', 'iter','MaxIter',1e4,'MaxFunEvals',1e4);
new_p = curr_p;
likelihood = @(x,sigma)(sum(residuals(x,sigma,true))*(-1));

q_prev = prior(curr_p)+likelihood(curr_p([H.PopulationParams H.IndividualParams.Index]),...
    exp(curr_p(setdiff(H.SigmaParams,H.IndividualParams.OmegaIndex))));

ub = log([PI.par([H.PopulationParams H.IndividualParams.Index]).maxValue]);
lb = log([PI.par([H.PopulationParams H.IndividualParams.Index]).minValue]);
options_fminsearch=optimset('Display','iter','MaxFunEvals', 1e4, 'MaxIter',1e4, 'TolFun', 1e-4);


while errTol>p.delta
    %% Expectation
    % Evaluate integral by MCMC MH sampling
    
    E_likelihood = @(x) likelihood([curr_p([H.PopulationParams]) x], ...
        exp(curr_p(setdiff(H.SigmaParams,H.IndividualParams.OmegaIndex))));
    E_prior = @(x) prior([curr_p([H.PopulationParams]) x curr_p(H.SigmaParams)]);
    
    [randeffects, logP, ~] = mcmc_mh(curr_p([random_indx{:,:}]),E_likelihood, E_prior, m,'StepSize',p.StepSize);
    q_k = mean(logP(ceil(size(randeffects,1)/2):end));
    E_Z = mean(randeffects(ceil(size(randeffects,1)/2):end,:));
    
    if ~isempty(H.SigmaParams)
        E_prior = @(x) prior([new_p([H.PopulationParams]) E_Z x]);
        PI = p.OutputFn(new_p);
        
        sigma_indx = size(H.IndividualParams.OmegaIndex,1)+1:length(H.SigmaParams);
        E_likelihood = (@(x)sum(getErrors(PI,exp(x(sigma_indx))))*(-1));
        [sigma, logP_sigma, ~] = mcmc_mh(new_p(H.SigmaParams),...
            E_likelihood, E_prior, m,'StepSize', p.StepSize);
        
        E_sigma = mean(sigma(ceil(size(sigma,1)/2):end,:));
        q_sigma = mean(logP_sigma(ceil(size(sigma,1)/2):end));
         
         q_k = q_sigma;
    end
    
       q_hat = q_prev + gamma*(q_k - q_prev);

       if q_hat < q_prev
           q_hat = q_prev;
           E_Z = curr_p([H.IndividualParams(:).Index]);
           E_sigma = curr_p([H.SigmaParams]);

       end
          new_p([random_indx{:,:}]) = E_Z;
          new_p(H.SigmaParams) = E_sigma;

    %% Maximization
     if ~isempty(H.SigmaParams)
        M_likelihood = @(x) likelihood([x(1:length(H.PopulationParams)) E_Z], ...
            exp(new_p(setdiff(H.SigmaParams,H.IndividualParams.OmegaIndex))));
        M_prior = @(x) prior([x(1:length(H.PopulationParams)) E_Z E_sigma]);
     else
        M_likelihood = @(x) likelihood([x(1:length(H.PopulationParams)) E_Z]);
        M_prior = @(x) prior([x(1:length(H.PopulationParams)) E_Z]);
     end
     
    if strcmp(min_fx, 'fminunc')
        [fixedeffects, logL] = fminunc(@(x)( (M_likelihood(x)+M_prior(x))*(-1)),...
            curr_p([H.PopulationParams]),fminunc_options);
    elseif strcmp(min_fx, 'lsqnonlin')
        
        residuals_fn = @(x) residuals(x,exp(new_p(setdiff(H.SigmaParams,H.IndividualParams.OmegaIndex))), false);
        [p_hat, ~] = lsqnonlin(residuals_fn,new_p([H.PopulationParams [H.IndividualParams(:).Index]]), lb,ub, options_fminsearch);
        E_Z = p_hat( [H.IndividualParams(:).Index]);
        fixedeffects = p_hat(H.PopulationParams);
        M_likelihood = @(x) likelihood([x(1:length(H.PopulationParams)) E_Z], ...
            exp(new_p(setdiff(H.SigmaParams,H.IndividualParams.OmegaIndex))));
        M_prior = @(x) prior([x(1:length(H.PopulationParams)) E_Z E_sigma]);
        [fixedeffects, logL] = fminsearch(@(x)( (M_likelihood(x)+M_prior(x))*(-1)),...
            fixedeffects,fminunc_options);
    else
        strcmp(min_fx, 'fminsearch')
        [fixedeffects, logL] = fminsearch(@(x)( (M_likelihood(x)+M_prior(x))*(-1)),...
            curr_p([H.PopulationParams]),fminunc_options);

    end
    
    % Updating new values
    new_p(fixed_indx) = fixedeffects;
    errTol = abs(-logL - q_hat);
    q_prev = -logL;
    i = i+1;
    gamma = gamma/i;
end
parameters = curr_p;

return

