%% DREAM MCMC
N = 32;
finalValues = log([PI.par(:).finalValue]);
X0 =[ (finalValues); randn(100,length(finalValues))*0.1 + finalValues];
logL=rowfun(obj_fun,table(X0));
[L,I]=sort(logL{:,:});

w0=X0(I(1:N),:);

h.IndividualParams=[];
tic
[x, p_x,accept,pCR] = dreamHParallel(w0,likelihood_fun,prior_fun,...
    size(w0,1),ceil(1e6/size(w0,1)), length(finalValues), 'BurnIn', 2e5,'StepSize',...
    1.38,'H', h);
toc

tic
[x2, p_x2,accept2,pCR2] = dreamH(x(:,:,end)',likelihood_fun,prior_fun,...
    size(w0,1),ceil(1e6/size(w0,1)), length(finalValues), 'BurnIn', 2e5,'StepSize',...
    1.38,'H', h);
toc


%% Concatenating chains
x_a=cat(3,x,x2);
logP=[p_x; p_x2];
x_mat=x(:,:)';

%% Diagnostics
plotMCMCDiagnostics(x_a(:,:,1:ceil(1e6/13)), logP(1:ceil(1e6/13),:),'name',...
    {PI.par(:).name},'model', ' Kosinsky','plots', 'trace')
plotMCMCDiagnostics(x, p_x,'name', {PI.par(:).name})
getGelmanRubinStatistic(x_a,{PI.par(:).name})
%% Plotting results

postSamples =x(:,:,ceil(size(x,3)/5):200:2e4);
logP_thinned = p_x(ceil(size(x,3)/5):200:2e4,:);
plotMCMCDiagnostics(postSamples,logP_thinned,'name', {PI.par(:).name},'model', 'PK model (ThreeComp)');

plotMCMCDiagnostics(postSamples([H.PopulationParams],:,:),logP_thinned,...
    'name', {PI.par([H.PopulationParams]).name},'model', 'PK model (ThreeComp)');

postSamples=postSamples(:,:)';

% Population Parameters
plotBivariateMarginals_2(exp(postSamples(:,H.PopulationParams)),'names',{PI.par(H.PopulationParams).name})
% Individual and population parameters
plotBivariateMarginals_2(exp(postSamples(:,[H.IndividualParams.EtaIndex...
    H.IndividualParams.OmegaIndex H.IndividualParams.Index])),'names', {PI.par([H.IndividualParams.EtaIndex...
    H.IndividualParams.OmegaIndex H.IndividualParams.Index]).name})
% Sigma parameters
plotBivariateMarginals_2(exp(postSamples(:, H.SigmaParams)),'names', {PI.par(H.SigmaParams).name})
%% Posterior predictions
simFun=@(x)getOutput(PI,@(p)sim(p,PI.tspan(end),u,PI.tspan),x,...
    @(p)getPhi2(p,H,length(u),'initialValue',x_0),normIndx, H);
tic
PI=getPosteriorPredictions(exp(postSamples),PI,simFun,observablesPlot);
toc
PI=getCredibleIntervals(PI,observables, exp(postSamples),H);
plotPosteriorPredictions(PI,observables)

%% Save results
save(strjoin({cd '/DREAM_MCMC_x.mat'},''), 'x')
save(strjoin({cd '/DREAM_MCMC_p_x.mat'},''), 'p_x')

save('/Users/migueltenorio/Documents/MATLAB/SimBiology/CIM/output/PI/DREAM_MCMC_logP2.mat', 'p_x')