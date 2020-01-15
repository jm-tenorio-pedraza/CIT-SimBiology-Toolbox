%% DREAM MCMC
finalValues = log([PI.par(:).finalValue]);
N = length(finalValues);

X0 =[ (finalValues); randn(100,length(finalValues))*0.2 + finalValues];
logL=rowfun(obj_fun,table(X0));
[L,I]=sort(logL{:,:});

w0=X0(I(1:N),:);

h.IndividualParams=[];
tic
[x, p_x,accept,pCR] = dreamHParallel(w0,likelihood_fun,prior_fun,...
    size(w0,1),ceil(2e6/size(w0,1)), length(finalValues), 'BurnIn', ...
    4e5,'StepSize',2.38,'H', h);
toc


%% Diagnostics
plotMCMCDiagnostics(x,p_x,'name', paramNames,'model',...
    PI.model,'interpreter', 'tex')
%% Plotting results
delta = 1500;
indx1 = ceil(4e5/size(x,1)+1):delta:size(x,3);
indx2 = (size(x,3)+ceil(2e5/size(x2,1))+1):delta:(size(x2,3)+size(x,3));
indx3 = (indx2(end)+ceil(2e5/size(x3,1))+1):delta:(size(x3,3)+indx2(end));

indx = [indx1];
[mean_w, w_indx] = sort(mean(p_x(indx,:)));

postSamples =x(:,w_indx(1:end),indx);
logP_thinned = p_x(indx,w_indx(1:end));
plotMCMCDiagnostics(postSamples,logP_thinned,'name',...
    paramNames,'model', PI.model,'interpreter','tex');

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
% Sigma parameters
plotBivariateMarginals_2(exp(postSamples(:, PI.H.SigmaParams)),'names',...
    {PI.par(H.SigmaParams).name})
plotIIVParams(postSamples, PI)
%% Posterior predictions
simFun=@(x)getOutput(PI,@(p)sim(p,PI.tspan(end),PI.u,PI.tspan),x,...
    @(p)getPhi2(p,PI.H,length(PI.u),'initialValue',PI.x_0),PI.normIndx, PI.H);
tic
PI=getPosteriorPredictions(exp(postSamples),PI,simFun,observables);
toc
PI=getCredibleIntervals(PI,observables, exp(postSamples),PI.H, 'logit_indx', 2:6);
plotPosteriorPredictions(PI,observables,'output','indiv')

%% Posterior credible intervals
 PI=mcmcCI(PI, (postSamples), logP_thinned', 0.95,'method', 'symmetric');
 plotCI(PI, 'TwoComp')
%% Save results
save(strjoin({cd '/DREAM_MCMC_x.mat'},''), 'xa')
save(strjoin({cd '/DREAM_MCMC_p_x.mat'},''), 'p_xa')

load(strjoin({cd '/DREAM_MCMC_x.mat'},''))
load(strjoin({cd '/DREAM_MCMC_p_x.mat'},''))
