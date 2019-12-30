%% DREAM MCMC
finalValues = log([PI.par(:).finalValue]);
N = length(finalValues);

X0 =[ (finalValues); randn(100,length(finalValues))*0.02 + finalValues];
logL=rowfun(obj_fun,table(X0));
[L,I]=sort(logL{:,:});

w0=X0(I(1:N),:);

h.IndividualParams=[];
tic
[x, p_x,accept,pCR] = dream(w0,likelihood_fun,prior_fun,...
    size(w0,1),ceil(1e6/size(w0,1)), N, 'BurnIn', 2e5,'StepSize',...
    2.38,'H',h);
toc

tic
[x2, p_x2,accept2,pCR2] = dream(x(:,:,end)',likelihood_fun,prior_fun,...
    size(w0,1),ceil(1e6/size(w0,1)), N, 'BurnIn', 2e5,'StepSize',...
    2.38,'H', h);
toc


%% Concatenating chains
x_a=cat(3,x,x2);
logP=[p_x; p_x2];
x_mat=x(:,:)';

%% Diagnostics
plotMCMCDiagnostics(x, p_x,'name', {PI.par(:).name},'model', PI.model)
plotMCMCDiagnostics(x([PI.H.PopulationParams],:,:),p_x,...
    'name', {PI.par([PI.H.PopulationParams]).name},'model', PI.model);

%% Plotting results
indx = ceil(2e5/size(x,2)):6e2:size(x,3);
mean_px = mean(p_x(ceil(2e5/size(x,2)):end,:));
[mean_px_sorted,w_indx] =sort(mean_px);
[~,max_indx] = max(p_x(end,:));
w_indx = w_indx(end:-1:end-ceil(size(x,2)/2));
postSamples =x(:,w_indx,indx);
logP_thinned = p_x(indx,w_indx);
plotMCMCDiagnostics(postSamples,logP_thinned,'name',...
    {PI.par(:).name},'model', PI.model);

plotMCMCDiagnostics(postSamples([PI.H.PopulationParams],:,:),logP_thinned,...
    'name', {PI.par([PI.H.PopulationParams]).name},'model', PI.model);

plotMCMCDiagnostics(postSamples([PI.H.CellParams.Index],:,:),logP_thinned,...
    'name', {PI.par([PI.H.CellParams.Index]).name},'model', PI.model);

plotMCMCDiagnostics(postSamples([PI.H.SigmaParams],:,:),logP_thinned,...
    'name', {PI.par([PI.H.SigmaParams]).name},'model', PI.model);

postSamples=postSamples(:,:)';
logP_thinned=reshape(logP_thinned',1,[]);
%% 
% Population Parameters
plotBivariateMarginals_2((postSamples(:,[PI.H.PopulationParams])),...
    'names',paramNames([PI.H.PopulationParams]),'interpreter', 'tex')
% Individual and population parameters
plotBivariateMarginals_2((postSamples(:, [PI.H.CellParams.Index PI.H.IndividualParams.Index])),...
    'names',paramNames([PI.H.CellParams.Index PI.H.IndividualParams.Index]),'interpreter', 'tex')
% Sigma parameters
plotBivariateMarginals_2(exp(postSamples(:, PI.H.SigmaParams)),'names',...
    paramNames(PI.H.SigmaParams), 'interpreter', 'tex')

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
PI=getCredibleIntervals(PI,observables, exp(postSamples),PI.H);
plotPosteriorPredictions(PI,observables)

%% Posterior credible intervals
 PI=mcmcCI(PI, exp(postSamples), logP_thinned', 0.95);
 plotCI(PI, PI.model, 'name', paramNames, 'interpreter', 'tex')
%% Save results
save(strjoin({cd '/DREAM_MCMC_x_red2.mat'},''), 'x')
save(strjoin({cd '/DREAM_MCMC_p_x_red2.mat'},''), 'p_x')

load(strjoin({cd '/DREAM_MCMC_x_red2.mat'},''))
load(strjoin({cd '/DREAM_MCMC_p_x_red2.mat'},''))
