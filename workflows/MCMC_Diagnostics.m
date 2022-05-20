%% MCMC Diagnostics

x = cat(3, x3,x4(:,:,2:end), x5(:,:,2:end),x6(:,:,2:end), x7(:,:,2:end));
p_x = [p_x3;p_x4(2:end,:);p_x5(2:end,:); p_x6(2:end,:);p_x7(2:end,:)];
%% Diagnostics

plotMCMCDiagnostics(x([PI.H.PopulationParams PI.H.SigmaParams],:,:),...
    p_x,'name', PI.paramNames([PI.H.PopulationParams PI.H.SigmaParams]),...
    'model', PI.model, 'interpreter', 'tex');
plotMCMCDiagnostics(x([PI.H.PopulationParams],:,:),...
    p_x,'name', PI.paramNames([PI.H.PopulationParams ]),...
    'model', PI.model, 'interpreter', 'tex');

plotMCMCDiagnostics(x([PI.H.CellParams(:).Index PI.H.IndividualParams(:).Index],:,:),...
    p_x,'name', paramNames([PI.H.CellParams(:).Index PI.H.IndividualParams(:).Index]),...
    'model', PI.model, 'interpreter', 'tex');

%% Plotting results
delta =1e2;
burnIn=0;
indx = ceil(burnIn/size(x,2)+1):delta:size(x,3);

[mean_w, w_indx] = sort(mean(p_x(indx,:)));

postSamples =x(:,w_indx(1:end),indx);
logP_thinned = p_x(indx,w_indx(1:end));
plotMCMCDiagnostics(postSamples,logP_thinned,'name', paramNames,'model',...
    'MOC1 tumor-bearing mice','interpreter', 'tex','Thinning',1,'plots',{'trace','corr','autocorr'})

postSamples=postSamples(:,:)';
logP_thinned=reshape(logP_thinned',1,[]);
%% 
etaIndx=[PI.H.CellParams.EtaIndex PI.H.IndividualParams.EtaIndex PI.H.RespParams.EtaIndex];
psiIndx=[PI.H.CellParams.OmegaIndex PI.H.IndividualParams.OmegaIndex PI.H.RespParams.OmegaIndex];
xIndx=[PI.H.CellParams.Index PI.H.IndividualParams.Index PI.H.RespParams.Index];
% Population Parameters
plotBivariateMarginals_2(exp(postSamples(:,[PI.H.PopulationParams])),...
    'names',paramNames([PI.H.PopulationParams ]),'interpreter', 'tex')
% Individual and population parameters
plotBivariateMarginals_2((postSamples(:, [etaIndx xIndx psiIndx])),...
    'names',paramNames([etaIndx xIndx psiIndx]),'interpreter', 'tex')
% Sigma Parameters
plotBivariateMarginals_2((postSamples(:,[PI.H.SigmaParams])),...
    'names',paramNames([PI.H.SigmaParams ]),'interpreter', 'tex')

plotIIVParams(postSamples, PI,'name', paramNames,'newFig', false,...
    'n_row', 2,'n_col',2,'figIndx', 1,'panel', true,'dim',true)
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
    'logit_indx', [],'simTime', simTime,'errorModel','additive');

%% Plot all
figure('Position', [10 10 1.5e3 1e3])
ncol = ceil(sqrt(length(PI.observablesFields)));
nrow = ceil(length(PI.observablesFields)/ncol);

for i=1:length(PI.observablesFields)
      subplot(nrow,ncol,i)

    plotPosteriorPredictions(PI,i,'outputs','group',...
       'newFig', false, 'TimeUnit', 'hours','color', 'dataset',...
        'simTime', simTime, 'YScale', 'linear', 'interpreter', 'tex','plot','data','indivData',false)
end

for i=1:9
    subplot(3,3,i)
    set(gca, 'YScale', 'linear')
end
%% Plot individual variables
for i =1:length(PI.observablesPlot)
 plotPosteriorPredictions(PI,i,'outputs','indiv', ...
        'newFig', true, 'TimeUnit', 'hours','color', 'cell','simTime', ...
        simTime, 'YScale', 'linear','indivData',false)


end
%% Plot prediction errors
for i =1:length(PI.observablesPlot)
plotPosteriorErrors(PI, i, 'outputs','indiv', 'newFig', true, 'TimeUnit' ,...
    'hours', 'color', 'dataset')
end
%% Posterior credible intervals
PI=mcmcCI(PI, (postSamples), logP_thinned', 0.95,'method', 'symmetric');
plotCI(PI, 'TwoComp', 'name', PI.paramNames, 'interpreter', 'tex')
PI.postSamples = postSamples;
PI.logP = logP_thinned;
plotHistogram(PI.postSamples(:,[PI.H.PopulationParams]), PI.paramNames([PI.H.PopulationParams]))
