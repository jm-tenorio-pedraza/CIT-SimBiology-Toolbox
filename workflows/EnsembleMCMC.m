%% Ensemble MCMC
% Initial walkers
% priorL=rowfun(prior_fun,table(w0));
% w0=w0(priorL{:,:}>-Inf,:);
w0=[finalValues;randn(300,size(finalValues,2))*0.1+repmat(finalValues,300,1)];
logL=rowfun(obj_fun,table(w0));
[~,I]=sort(logL{:,:});

w0=w0(I(1:size(w0,2)*2),:);

% MCMC sampler
tic
[models,logP]=gwmcmc(w0',{prior_fun likelihood_fun},1e5, 'StepSize', 2.5,...
'ThinChain',1,'BurnIn',0);
toc

tic
[models2,logP2]=gwmcmc(models(:,:,end),{prior_fun likelihood_fun},9e5, 'StepSize', 2.5,...
'ThinChain',1,'BurnIn',0);
toc

tic
[models3,logP3]=gwmcmc(models2(:,:,end),{prior_fun likelihood_fun},9e5, 'StepSize', 2.0,...
'ThinChain',1,'BurnIn',0);
toc

tic
[models4,logP4]=gwmcmc(models3(:,:,end),{prior_fun likelihood_fun},9e5, 'StepSize', 1.5,...
'ThinChain',1,'BurnIn',0);
toc

tic
[models5,logP5]=gwmcmc(models4(:,:,end),{prior_fun likelihood_fun},9e5, 'StepSize', 1.25,...
'ThinChain',1,'BurnIn',0);
toc

x= cat(3,models,models2,models3,models4,models5);
logL=cat(3,logP,logP2,logP3,logP4,logP5);

% Plotting entire chain
plotMCMCDiagnostics(x,logL,PI,'model','Kosinsky (General)')

% Plotting thinned out samples
postSample=models_array(:,:,1.8e4:1000:end);
postSample=postSample(:,:)';

logPSample=logL(:,:,1.8e4:1:end);
logPSample=sum(logPSample(:,:))';

plotMCMCDiagnostics(postSample,logPSample,PI,'model', 'Kosisnky (kpro)')

% Bivariate marginal plots
plotBivariateMarginals_2(exp(postSample(:,H.PopulationParams)),PI)
plotBivariateMarginals_2(exp(postSample(:,H.IndividualParams(1).Index)),PI,'names',...
    {PI.par(H.IndividualParams(1).Index).name})
plotBivariateMarginals_2(exp(postSample(:,H.SigmaParams)),PI,'names', {PI.par(H.SigmaParams).name})

save('models_Kosinksy_MOC1_general.mat','x')
save('postlogP_Kosinksy_MOC1_general.mat','logL')

%% Posterior predictions
simFun=@(x)getOutput(PI,@(p)sim(p,100,u,1:1:100),x,...
    @(p)getPhi2(p,H,length(u)), 4:6);
tic
models_array=models(:,:)';
PI=getPosteriorPredictions(exp(models_array(6e5:400:end,:)),PI,simFun,observables);
toc
PI=getCredibleIntervals(PI,observables, exp(models_array(6e5:400:end,:)),H);
plotPosteriorPredictions(PI,observables)
plotGroupPosteriorPredictions(PI,observables(1),exp(postSample),H)


tic 
randsample(1:3,1,true)
toc
