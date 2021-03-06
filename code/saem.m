function [parameters, logL] = saem(p0, residuals, likelihood,prior, H, PI,varargin)

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
p.addParameter('StepSize', 2.38^2)
p.addParameter('MinFunc', 1)
p.addParameter('OutputFn', [])
p.addParameter('SimFn', [])
p.addParameter('BurnIn',0.5)
p.addParameter('Thinning',20)
p.addParameter('estimateSigma',false)


p.parse(varargin{:});
p=p.Results;


errTol = 1;
curr_p = p0;                                                                % Original log-transformed vector of parameters
m = p.m;
gamma = p.gamma;
i = 1;
min_fx = p.MinFunc;
burnIn = p.BurnIn;
thin=p.Thinning;
random_indx = arrayfun(@(x) x.Index, H.CellParams, 'UniformOutput',false);  % Identify vector parameter indexes of cellular parameters
random_indx(end+1:end+length(H.IndividualParams)) = arrayfun(@(x) x.Index,...% Identify vector parameter indexes of individual parameters
    H.IndividualParams, 'UniformOutput', false);
fixed_indx = [H.PopulationParams];                                          % Identify vector parameter indexes of population parameters

sigma_indx = setdiff(H.SigmaParams, [H.CellParams.OmegaIndex...             % Identify vector parameter indexes of sigma parameters
    H.IndividualParams.OmegaIndex]);

fminunc_options = optimset('Display', 'iter','MaxIter',1e4,...
    'MaxFunEvals',1e4);
anneal_options.Verbosity=2;
anneal_options.InitTemp=100;

q_prev = prior(curr_p)+likelihood(curr_p);                                  % Evaluate log-posterior likelihood at initial point

ub = log([PI.par([H.PopulationParams]).maxValue]);                          % Upper limit on population parameter values
lb = log([PI.par([H.PopulationParams]).minValue]);
options_fminsearch=optimset('Display','iter','MaxFunEvals', 1e4,...
    'MaxIter', 1e4, 'TolFun', 1e-4);

while errTol>p.delta
    %% Expectation
    % Evaluate integral by MCMC MH sampling
    E_likelihood = @(x) likelihood([curr_p([H.PopulationParams]) x, ...    % Likelihood function of the individual parameters
        curr_p(H.SigmaParams)]);
    
    E_prior = @(x) prior([curr_p([H.PopulationParams]) x ...                % Prior function of the individual params
        curr_p(H.SigmaParams)]);
    
    [randeffects, logP, ~] = mcmc_mh(curr_p([random_indx{:,:}]),...         % Metropolis-Hastings MCMC
        E_likelihood, E_prior, m,'StepSize',p.StepSize/...
        (length([random_indx{:,:}])),'BurnIn', burnIn);
    
    q_k = mean(logP((burnIn*m)+1:thin:end));                                     % Calculate the mean log-posterior likelihood and the expected value
    E_Z = mean(randeffects(burnIn*m+1:thin:end,:));
    
    if p.estimateSigma                                              % Integrate over the sigmas
        E_prior = @(x) prior([curr_p([H.PopulationParams]) E_Z x]);
        PI = p.OutputFn(curr_p);
        
        sigma_index = size([H.CellParams.OmegaIndex H.IndividualParams.OmegaIndex ],2)...
            +1:length(H.SigmaParams);
        
        E_likelihood = (@(x)sum(getErrors(PI,exp(x(sigma_index))))*(-1));
        [sigma, logP_sigma, ~] = mcmc_mh(curr_p(H.SigmaParams),...
            E_likelihood, E_prior, m,'StepSize', p.StepSize/length(H.SigmaParams),...
            'BurnIn', burnIn);
        
        E_sigma = mean(sigma(burnIn*m+1:thin:end,:));
        q_sigma = mean(logP_sigma(burnIn*m+1:thin:end));
         
         q_k = q_sigma;
         curr_p(H.SigmaParams) = E_sigma;

    end
    
       q_hat = q_prev + gamma*(q_k - q_prev);

       if q_hat < q_prev
           q_hat = q_prev;
           E_Z = curr_p([random_indx{:,:}]);

       end
          curr_p([random_indx{:,:}]) = E_Z;

    %% Maximization
     if p.estimateSigma
        M_likelihood = @(x) likelihood([x(H.PopulationParams) E_Z, ...
            E_sigma]);
        M_prior = @(x) prior([x(H.PopulationParams) E_Z E_sigma]);
        p0 = curr_p([H.PopulationParams]);
     else
        M_likelihood = @(x) likelihood([x(1:length(H.PopulationParams)) E_Z x(length(H.PopulationParams)+1:end)]);
        M_prior = @(x) prior([x(1:length(H.PopulationParams)) E_Z x(length(H.PopulationParams)+1:end)]);
        p0 = curr_p([H.PopulationParams H.SigmaParams]);
        fixed_indx = [H.PopulationParams H.SigmaParams];
     end
     
    if strcmp(min_fx, 'fminunc')
        [fixedeffects, logL] = fminunc(@(x)( (M_likelihood(x)+M_prior(x))*(-1)),...
            p0,fminunc_options);
    elseif strcmp(min_fx, 'lsqnonlin')

        residuals_fn = @(x) residuals([x E_Z],exp(curr_p(sigma_indx)));
        
        [p_hat, ~] = lsqnonlin(residuals_fn,curr_p([H.PopulationParams...
            ]), lb,ub, options_fminsearch);
        fixedeffects = p_hat(H.PopulationParams);
       curr_p(fixed_indx) = fixedeffects;
       logL = (likelihood(curr_p)+prior(curr_p))*(-1);
    else
        if strcmp(min_fx, 'fminsearch')
        [fixedeffects, logL] = fminsearch(@(x)( (M_likelihood(x)+M_prior(x))*(-1)),...
            p0,fminunc_options);
        elseif strcmp(min_fx, 'anneal')
              [fixedeffects, logL] = anneal(@(x)( (M_likelihood(x)+M_prior(x))*(-1)),...
            p0,anneal_options);
        end
        curr_p(fixed_indx) = fixedeffects;
                
    end
    % Updating new values
    
    q_hat = q_prev + gamma*(-logL - q_prev);

    errTol = abs(q_hat -q_prev);
    q_prev = q_hat;
    i = i+1;
    gamma = gamma/i;
end
parameters = curr_p;

return

