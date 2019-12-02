%% DREAM MCMC
finalValues = log([PI.par(:).finalValue]);
N = length(finalValues);

X0 =[ (finalValues); randn(100,length(finalValues))*0.1 + finalValues];
logL=rowfun(obj_fun,table(X0));
[L,I]=sort(logL{:,:});

w0=X0(I(1:N),:);

h.IndividualParams=[];
tic
[x, p_x,accept,pCR] = dreamHParallel(w0,likelihood_fun,prior_fun,...
    size(w0,1),ceil(5e5/size(w0,1)), length(finalValues), 'BurnIn', 1e5,'StepSize',...
    1.38,'H', h);
toc

tic
[x2, p_x2,accept2,pCR2] = dreamHParallel(x(:,:,end)',likelihood_fun,prior_fun,...
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
plotMCMCDiagnostics(x, p_x,'name', {PI.par(:).name},'model', 'CIM')
getGelmanRubinStatistic(x_a,{PI.par(:).name})
%% Plotting results

postSamples =x(:,:,4e3:3e2:end);
logP_thinned = p_x(4e3:3e2:end,:);
plotMCMCDiagnostics(postSamples,logP_thinned,'name',...
    {PI.par(:).name},'model', 'CIM');

plotMCMCDiagnostics(postSamples([H.PopulationParams],:,:),logP_thinned,...
    'name', {PI.par([H.PopulationParams]).name},'model', 'PK model (TwoComp)');

postSamples=postSamples(:,:)';
%% 
% Population Parameters
plotBivariateMarginals_2((postSamples(:,[PI.H.PopulationParams PI.H.SigmaParams])),...
    'names',paramNames([PI.H.PopulationParams PI.H.SigmaParams]))
% Individual and population parameters
plotBivariateMarginals_2((postSamples(:, [PI.H.CellParams.Index PI.H.IndividualParams.Index])),...
    'names',paramNames([PI.H.CellParams.Index PI.H.IndividualParams.Index]))
% Sigma parameters
plotBivariateMarginals_2(exp(postSamples(:, PI.H.SigmaParams)),'names',...
    {PI.par(H.SigmaParams).name})

figure
hold on
for i=1:length(PI.data)
    histogram(exp(postSamples(:,PI.H.IndividualParams(1).Index(i))),...
        'FaceColor',PI.data(i).colors,'FaceAlpha', 0.5, 'Normalization', 'probability')
    
end
% histogram(exp(postSamples(:,PI.H.IndividualParams(1).EtaIndex)),...
%     'FaceColor',PI.data(i).colors,'FaceAlpha', 0.5, 'Normalization', 'probability')
legend({PI.data(:).Name},'interpreter', 'none')
ylabel('prob')
xlabel('Deviations wrt mean parameter')
title('Inter-individual variation in K_{CD8}')
%% Posterior predictions
simFun=@(x)getOutput(PI,@(p)sim(p,PI.tspan(end),PI.u,PI.tspan),x,...
    @(p)getPhi2(p,PI.H,length(PI.u),'initialValue',PI.x_0),PI.normIndx, PI.H);
tic
PI=getPosteriorPredictions(exp(postSamples),PI,simFun,observables);
toc
PI=getCredibleIntervals(PI,observables, exp(postSamples),H);
plotPosteriorPredictions(PI,observables)

%% Save results
save(strjoin({cd '/DREAM_MCMC_x_1.mat'},''), 'x')
save(strjoin({cd '/DREAM_MCMC_p_x_1.mat'},''), 'p_x')

load(strjoin({cd '/DREAM_MCMC_x_4.mat'},''))
load(strjoin({cd '/DREAM_MCMC_p_x_4.mat'},''))
