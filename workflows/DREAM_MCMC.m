%% DREAM MCMC
finalValues = log([PI.par(:).finalValue]);
N = length(finalValues);

X0 =[ (finalValues); randn(100,length(finalValues))*0.2 + finalValues];
logL=rowfun(obj_fun,table(X0));
[L,I]=sort(logL{:,:});

w0=X0(I(1:N),:);

h.IndividualParams=[];
tic
[x, p_x,accept,pCR] = dreamHParallel(w0,likelihood_fun,prior_fun_MCMC,...
    size(w0,1),ceil(1e6/size(w0,1)), length(finalValues), 'BurnIn', ...
    2e5,'StepSize',2.38,'H', h);
toc

tic
[x2, p_x2,accept2,pCR2] = dreamHParallel(x(:,:,end)',likelihood_fun,prior_fun_MCMC,...
    size(w0,1),ceil(1e6/size(w0,1)), length(finalValues), 'BurnIn', ...
    2e5,'StepSize',2.38,'H', h);
toc
x=cat(3,x,x2);
p_x = [p_x; p_x2];
%% Diagnostics
plotMCMCDiagnostics(x,p_x,'name', paramNames,'model',...
    PI.model,'interpreter', 'tex')

plotMCMCDiagnostics(x([PI.H.PopulationParams PI.H.SigmaParams],:,:),...
    p_x,'name', paramNames([PI.H.PopulationParams PI.H.SigmaParams]),...
    'model', PI.model, 'interpreter', 'tex');

%% Plotting results
delta = 1e3;
burnIn=2e5;
indx = ceil(burnIn/size(x,1)+1):delta:size(x,3);

[mean_w, w_indx] = sort(mean(p_x(indx,:)));

postSamples =x(:,w_indx(1:end),indx);
logP_thinned = p_x(indx,w_indx(1:end));

plotMCMCDiagnostics(postSamples([PI.H.PopulationParams PI.H.SigmaParams],:,:),...
    logP_thinned,'name', paramNames([PI.H.PopulationParams PI.H.SigmaParams]),...
    'model', PI.model, 'interpreter', 'tex');
plotMCMCDiagnostics(postSamples([PI.H.CellParams.Index PI.H.IndividualParams.Index],:,:),...
    logP_thinned,'name', paramNames([PI.H.CellParams.Index PI.H.IndividualParams.Index]),...
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
% Sigma parameters
plotBivariateMarginals_2(exp(postSamples(:, PI.H.SigmaParams)),'names',...
    {PI.par(H.SigmaParams).name})
plotIIVParams(postSamples, PI,'name', paramNames)
%% Posterior predictions
simFun=@(x)getOutput(PI,@(p)sim(p,PI.tspan(end),PI.u,PI.tspan),x,...
    @(p)getPhi3(p,PI.H,length(PI.u),'initialValue',PI.x_0),PI.normIndx, PI.H);
tic
PI=getPosteriorPredictions(exp(postSamples),PI,simFun,PI.observablesPlot);
toc
PI=getCredibleIntervals(PI,PI.observablesPlot, exp(postSamples),PI.H, 'logit_indx', []);
plotPosteriorPredictions(PI,PI.observablesPlot,'output','indiv','all', true)

%% Posterior credible intervals
PI=mcmcCI(PI, (postSamples), logP_thinned', 0.95,'method', 'symmetric');
plotCI(PI, 'TwoComp', 'name', paramNames, 'interpreter', 'tex')
PI.postSamples = postSamples;
plotHistogram(PI.postSamples(:,[PI.H.PopulationParams]),...
    paramNames([PI.H.PopulationParams]))
plotHistogram(PI.postSamples(:,[PI.H.CellParams.Index]),...
    paramNames([PI.H.CellParams.Index ]))
plotHistogram(PI.postSamples(:,[PI.H.IndividualParams.Index]),...
    paramNames([PI.H.IndividualParams.Index ]))
%% Save results
save(strjoin({cd '/PI_CIM7_ICB_1_DREAM_MCMC_x.mat'},''), 'x')
save(strjoin({cd '/PI_CIM7_ICB_1_DREAM_MCMC_p_x.mat'},''), 'p_x')

load(strjoin({cd '/CIM_red2_DREAM_MCMC_x.mat'},''))
load(strjoin({cd '/CIM_red2_DREAM_MCMC_p_x.mat'},''))
