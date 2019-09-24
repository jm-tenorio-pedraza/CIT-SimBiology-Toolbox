%% DREAM MCMC

finalValues = [PI.par(:).finalValue];
X0 =[ log(finalValues); randn(100,length(finalValues))*0.1 + log(finalValues)];
logL=rowfun(obj_fun,table(X0));
[L,I]=sort(logL{:,:});

w0=X0(I(1:length(finalValues)*1),:);

tic
[x, p_x,accept] = dreamH(w0,likelihood_fun,prior_fun,...
    size(w0,1),ceil(1e6/size(w0,1)), length(finalValues), 'BurnIn', 1e5,'StepSize', 1.1,'H', H);
toc

tic
[x2, p_x2,accept2] = dreamH(x(:,:,end)',likelihood_fun,prior_fun,...
    size(w0,1),ceil(1e6/size(w0,1)), length(finalValues), 'BurnIn', 0, 'StepSize', 1.1,'H',H);
toc

tic
[x3, p_x3,accept3] = dreamH(x2(:,:,end)',likelihood_fun,prior_fun,...
    size(w0,1),ceil(1e6/size(w0,1)), length(finalValues), 'BurnIn', 2e5, 'StepSize', 0.8,'H',H);
toc

%% Concatenating chains
x_a=cat(3,x,x2);
logP=[p_x; p_x2];
x_mat=x_a(:,:)';

%% Diagnostics
plotMCMCDiagnostics(x_a, logP,'name', {PI.par(:).name})
plotMCMCDiagnostics(x, p_x,'name', {PI.par(:).name})

%% Plotting results

postSamples = x_mat(2e5:1e3:end,:);
plotMCMCDiagnostics(x,p_x,'name', {PI.par(:).name},'model', 'CIM');
plotBivariateMarginals_2((postSamples(:,H.PopulationParams)),'names',{PI.par(H.PopulationParams).name})
plotBivariateMarginals_2(exp(postSamples(:,[H.IndividualParams.EtaIndex H.IndividualParams.OmegaIndex])),...
    'names',{PI.par([H.IndividualParams.EtaIndex H.IndividualParams.OmegaIndex]).name})
plotBivariateMarginals_2(exp(postSamples(:,[H.IndividualParams.EtaIndex...
    H.IndividualParams.Index])),'names', {PI.par([H.IndividualParams.EtaIndex...
    H.IndividualParams.Index]).name})
plotBivariateMarginals_2(exp(postSamples(:, H.SigmaParams)),'names', {PI.par(H.SigmaParams).name})
plot(logP)
plot(lags,c)

%% Posterior predictions
simFun=@(x)getOutput(PI,@(p)sim(p,100,u,1:1:100),x,...
    @(p)getPhi2(p,H,length(u)), length(observables)-1:length(observables), 1:100);
tic
PI=getPosteriorPredictions(exp(postSamples),PI,simFun,observables);
toc
PI=getCredibleIntervals(PI,observables, exp(postSamples),H);
plotPosteriorPredictions(PI,observables)
plotGroupPosteriorPredictions(PI,observables(1),postSamples,H)

%% Save results
save('/Users/migueltenorio/Documents/MATLAB/SimBiology/CIM/output/PI/DREAM_MCMC_p.mat', 'x_a')
save('/Users/migueltenorio/Documents/MATLAB/SimBiology/CIM/output/PI/DREAM_MCMC_logP.mat', 'logP')