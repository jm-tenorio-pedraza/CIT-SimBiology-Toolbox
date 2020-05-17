%% DREAM MCMC
finalValues = log([PI.par(:).finalValue]);
N = length(finalValues);

X0 =[ (finalValues); randn(100,length(finalValues))*0.001 + finalValues];
logL=rowfun(obj_fun,table(X0));
[L,I]=sort(logL{:,:});

w0=X0(I(1:N),:);

h.IndividualParams=[];
p=parpool('local')
tic
[x, p_x,accept,pCR,stepSize, J, n_id] = dreamHParallel(w0,likelihood_fun,prior_fun_MCMC,...
    size(w0,1),ceil(3e5/size(w0,1)), length(finalValues), 'BurnIn', ...
    5e5,'StepSize',2.38,'H', h);
toc
w0 = x(:,:,end);
tic
[x2, p_x2,accept2,pCR2,stepSize2,J2,n_id2] = dreamHParallel(w0',likelihood_fun,prior_fun_MCMC,...
    size(w0,1),ceil(2.5e5/size(w0,1)), length(finalValues), 'BurnIn', ...
    3e5,'StepSize',stepSize,'H', h,'pCR', pCR, 'J', J, 'n_id', n_id);
toc
w0 = x2(:,:,end);
poolobj = gcp('nocreate');
delete(poolobj)
clearvars x2 p_x2
tic
[x3, p_x3,accept3,pCR3,stepSize3,J3,n_id3] = dreamHParallel(w0',likelihood_fun,prior_fun_MCMC,...
    size(w0,1),ceil(3e5/size(w0,1)), length(finalValues), 'BurnIn', ...
    0,'StepSize',stepSize2,'H', h,'pCR', pCR2,'J', J2, 'n_id', n_id2);
toc

w0 = x3(:,:,end);
tic
[x4, p_x4,accept4,pCR4,stepSize4,J4, n_id4] = dreamHParallel(w0',likelihood_fun,...
    prior_fun_MCMC,size(w0,1),ceil(3e5/size(w0,1)), length(finalValues),...
    'BurnIn', 5e5,'StepSize',stepSize3,'H', h,'pCR', pCR3, 'J', J3, 'n_id', n_id3);
toc

w0 = x4(:,:,end);
tic
[x5, p_x5,accept5,pCR5,stepSize5,J5, n_id5] = dreamHParallel(w0',likelihood_fun,prior_fun_MCMC,...
    size(w0,1),ceil(2.5e5/size(w0,1)), length(finalValues), 'BurnIn', ...
    3e5,'StepSize',stepSize4,'H', h, 'pCR', pCR4, 'J', J4, 'n_id', n_id4);
toc
x=cat(3,x,x2);
p_x = [p_x; p_x2];
%% Diagnostics
plotMCMCDiagnostics(x,p_x,'name', paramNames,'model',...
    PI.model,'interpreter', 'tex')
plotMCMCDiagnostics(x3([PI.H.PopulationParams PI.H.SigmaParams],:,:),...
    p_x3,'name', paramNames([PI.H.PopulationParams PI.H.SigmaParams]),...
    'model', PI.model, 'interpreter', 'tex');
plotMCMCDiagnostics(x([PI.H.CellParams(:).Index PI.H.IndividualParams(:).Index],:,:),...
    p_x,'name', paramNames([PI.H.CellParams(:).Index PI.H.IndividualParams(:).Index]),...
    'model', PI.model, 'interpreter', 'tex');

%% Plotting results
delta = 1e3;
burnIn=0;
indx = ceil(burnIn/size(x3,1)+1):delta:size(x3,3);

[mean_w, w_indx] = sort(mean(p_x3(indx,:)));

postSamples =x3(:,w_indx(1:end),indx);
logP_thinned = p_x3(indx,w_indx(1:end));

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
%% Save results
save(strjoin({cd '/PI_CIM21_Control_11_2_DREAM_MCMC_x.mat'},''), 'x')
save(strjoin({cd '/PI_CIM21_Control_11_2_DREAM_MCMC_p_x.mat'},''), 'p_x')

load(strjoin({cd '/CIM_red2_DREAM_MCMC_x.mat'},''))
load(strjoin({cd '/CIM_red2_DREAM_MCMC_p_x.mat'},''))
