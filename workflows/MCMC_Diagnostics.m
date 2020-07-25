%% MCMC Diagnostics
%% Diagnostics
plotMCMCDiagnostics(x,p_x,'name', paramNames,'model',...
    PI.model,'interpreter', 'tex')
plotMCMCDiagnostics(x1([PI.H.PopulationParams PI.H.SigmaParams],:,:),...
    p_x1,'name', paramNames([PI.H.PopulationParams PI.H.SigmaParams]),...
    'model', PI.model, 'interpreter', 'tex');
plotMCMCDiagnostics(x([PI.H.CellParams(:).Index PI.H.IndividualParams(:).Index],:,:),...
    p_x,'name', paramNames([PI.H.CellParams(:).Index PI.H.IndividualParams(:).Index]),...
    'model', PI.model, 'interpreter', 'tex');

%% Plotting results
delta = 3e2;
burnIn=1e5;
indx = ceil(burnIn/size(x1,1)+1):delta:size(x1,3);

[mean_w, w_indx] = sort(mean(p_x1(indx,:)));

postSamples =x1(:,w_indx(1:end),indx);
logP_thinned = p_x1(indx,w_indx(1:end));

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
plotIIVParams(postSamples, PI,'name', paramNames,'newFig', false,...
    'n_row', 1,'n_col',1,'figIndx', 1,'panel', true,'dim',true)
%% Posterior predictions
simTime = unique([PI.tspan', 1:PI.tspan(end)]);
simFun=@(x)getOutput(PI,@(p)sim(p,PI.tspan(end),PI.u,simTime),x,...
    @(p)getPhi2(p,PI.H,length(PI.u),'initialValue',PI.x_0),PI.normIndx, PI.H,...
    'output', 'simoutput', 'simTime', simTime);
tic
PI=getPosteriorPredictions2(exp(postSamples(1:end,:)),PI,simFun,PI.observablesFields,...
    'simTime', simTime);
toc
PI=getCredibleIntervals(PI,PI.observablesFields, exp(postSamples(1:end,:)),PI.H,...
    'logit_indx', [],'simTime', simTime);
%% Plot all
figure('Position', [10 10 1.5e3 1e3])
ncol = ceil(sqrt(length(PI.observablesFields)));
nrow = ceil(length(PI.observablesFields)/ncol);

for i=1:length(PI.observablesFields)
     subplot(nrow,ncol,i)
    plotPosteriorPredictions(PI,i,'outputs','group',...
       'newFig', false, 'TimeUnit', 'days','color', 'dataset',...
        'simTime', simTime, 'YScale', 'linear', 'interpreter', 'tex','plot','data')
end

for i=1:9
    subplot(3,3,i)
    set(gca, 'YScale', 'linear')
end
%% Plot individual variables
for i =1:length(PI.observablesPlot)
 plotPosteriorPredictions(PI,i,'outputs','indiv', 'all', false,...
        'newFig', true, 'TimeUnit', 'days','color', 'cell','simTime', simTime)
end
%% Plot prediction errors
for i =1:length(PI.observablesPlot)
plotPosteriorErrors(PI, i, 'outputs','indiv', 'newFig', true, 'TimeUnit' ,...
    'hours', 'color', 'dataset')
end
%% Posterior credible intervals
PI=mcmcCI(PI, (postSamples), logP_thinned', 0.95,'method', 'symmetric');
plotCI(PI, 'TwoComp', 'name', paramNames, 'interpreter', 'tex')
PI.postSamples = postSamples;
PI.logP = logP_thinned;
plotHistogram(PI.postSamples(:,[PI.H.PopulationParams]), paramNames([PI.H.PopulationParams]))
