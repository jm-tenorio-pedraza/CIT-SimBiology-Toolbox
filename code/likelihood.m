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
    if par.indivData
        loglikelihood_rightCens = arrayfun(@(x)reshape(likelihood_censored(log(x.dataValue(x.rightCensIndx,1)),...
            log(x.simValue(ismember(x.simTime,x.dataTime(x.rightCensIndx)),1)),sigma(1),...
            'right'),1,[]),PI.IndivData,'UniformOutput',false);
        loglikelihood_leftCens = arrayfun(@(x)reshape(likelihood_censored(log(x.dataValue(x.leftCensIndx,1)),...
            log(x.simValue(ismember(x.simTime,x.dataTime(x.leftCensIndx)),1)),sigma(1),...
            'left'),1,[]),PI.IndivData,'UniformOutput',false);
    else
        loglikelihood_rightCens = arrayfun(@(x)reshape(likelihood_censored(log(x.dataValue(x.rightCensIndx,1)),...
            log(x.simValue(ismember(x.simTime,x.dataTime(x.rightCensIndx)),1)),sigma(1),...
            'right'),1,[]),PI.data,'UniformOutput',false);
        loglikelihood_leftCens = arrayfun(@(x)reshape(likelihood_censored(log(x.dataValue(x.leftCensIndx,1)),...
            log(x.simValue(ismember(x.simTime,x.dataTime(x.leftCensIndx)),1)),sigma(1),...
            'left'),1,[]),PI.data,'UniformOutput',false);
    end
    loglikelihood_rightCens = [loglikelihood_rightCens{:}];
        loglikelihood_leftCens = [loglikelihood_leftCens{:}];

    loglikelihood = loglikelihood+sum(loglikelihood_rightCens)+ sum(loglikelihood_leftCens);
else
end

return
