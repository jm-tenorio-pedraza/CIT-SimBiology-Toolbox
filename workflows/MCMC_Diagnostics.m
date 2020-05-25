%% MCMC Diagnostics
%% Diagnostics
plotMCMCDiagnostics(x,p_x,'name', paramNames,'model',...
    PI.model,'interpreter', 'tex')
plotMCMCDiagnostics(x11([PI.H.PopulationParams PI.H.SigmaParams],:,:),...
    p_x11,'name', paramNames([PI.H.PopulationParams PI.H.SigmaParams]),...
    'model', PI.model, 'interpreter', 'tex');
plotMCMCDiagnostics(x9([PI.H.CellParams(:).Index PI.H.IndividualParams(:).Index],:,:),...
    p_x9,'name', paramNames([PI.H.CellParams(:).Index PI.H.IndividualParams(:).Index]),...
    'model', PI.model, 'interpreter', 'tex');

%% Plotting results
delta = 1e3;
burnIn=10e5;
indx = ceil(burnIn/size(x,1)+1):delta:size(x,3);

[mean_w, w_indx] = sort(mean(p_x(indx,:)));

postSamples =x(:,w_indx(1:end),indx);
logP_thinned = p_x(indx,w_indx(1:end));

plotMCMCDiagnostics(postSamples([PI.H.PopulationParams PI.H.SigmaParams],:,:),...
    logP_thinned,'name', paramNames([PI.H.PopulationParams PI.H.SigmaParams]),...
    'model', PI.model, 'interpreter', 'tex');

postSamples=postSamples(:,:)';
logP_thinned=reshape(logP_thinned',1,[]);
%% 
% Population Parameters
plotBivariateMarginals_2((postSamples(:,[PI.H.PopulationParams PI.H.SigmaParams])),...
    'names',paramNames([PI.H.PopulationParams PI.H.SigmaParams]),'interpreter', 'tex')
% Individual and population parameters
plotBivariateMarginals_2((postSamples(:, [PI.H.CellParams.Index PI.H.IndividualParams.Index])),...
    'names',paramNames([PI.H.CellParams.Index PI.H.IndividualParams.Index]))
plotIIVParams(postSamples, PI,'name', paramNames)
%% Posterior predictions
simFun=@(x)getOutput(PI,@(p)sim(p,PI.tspan(end),PI.u,PI.tspan),x,...
    @(p)getPhi2(p,PI.H,length(PI.u),'initialValue',PI.x_0),PI.normIndx, PI.H);
tic
PI=getPosteriorPredictions(exp(postSamples),PI,simFun,PI.observablesPlot);
toc
PI=getCredibleIntervals(PI,PI.observablesPlot, exp(postSamples),PI.H, 'logit_indx', []);
plotPosteriorPredictions(PI,PI.observablesPlot,'output','group')

%% Posterior credible intervals
PI=mcmcCI(PI, (postSamples), logP_thinned', 0.95,'method', 'symmetric');
plotCI(PI, 'TwoComp', 'name', paramNames, 'interpreter', 'tex')
PI.postSamples = postSamples;
plotHistogram(PI.postSamples(:,[PI.H.PopulationParams]), paramNames([PI.H.PopulationParams]))
