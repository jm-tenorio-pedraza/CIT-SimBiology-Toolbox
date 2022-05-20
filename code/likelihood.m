function loglikelihood =likelihood(p,sim_fn,PI,varargin) 
par=inputParser;
par.addParameter('censoring', false)
par.addParameter('logTransform', true)
par.addParameter('errorModel', 'additive')
par.addParameter('constantVar',0.1);
par.addParameter('indivData',false);

par.parse(varargin{:})
par=par.Results;
[residuals,PI]=(getNormResiduals(p,@(x)sim_fn(x,PI.tspan(end),PI.u,PI.tspan),PI,...
    @(x)getPhi2(x,PI.H,size(PI.data,1),'initialValue',PI.x_0),...
   PI.normIndx,'log', par.logTransform,'errorModel',par.errorModel, 'constantVar',...
   par.constantVar,'indivData',par.indivData));
loglikelihood = sum(residuals*(-1));

% Censoring correction
if par.censoring
    sigma = p(setdiff(PI.H.SigmaParams,[PI.H.IndividualParams(:).OmegaIndex]));
    loglikelihood_cens = arrayfun(@(x)reshape(likelihood_censored(log(x.TV(end,1)),...
        log(x.simValue(x.simTime>=x.dataTime(end),1)),(sigma(1)),...
        x.censoring),1,[]),PI.data,'UniformOutput',false);
    loglikelihood_cens = [loglikelihood_cens{:}];
    loglikelihood = loglikelihood+sum(loglikelihood_cens);
else
end

return
