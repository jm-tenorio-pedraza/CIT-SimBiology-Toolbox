% DREAM MCMC

% X0 = postSample(1:30:end,[1:9 9:22 22:end]);
finalValues = [PI.par(:).finalValue];
X0 =[ log(finalValues); randn(100,length(finalValues))*0.1 + log(finalValues)];
logL=rowfun(obj_fun,table(X0));
[~,I]=sort(logL{:,:});

w0=X0(I(1:size(X0,2)*2),:);


tic
[x, p_x,accept] = dream(X0,@(p)(likelihood_fun(p)+prior_fun(p)),...
    size(X0,1),ceil(2e6/size(X0,1)), length(finalValues), 'BurnIn', 1e5,'StepSize', 1.1);
toc

tic
[x2, p_x2,accept2] = dream(x(:,:,end)',@(p)(likelihood_fun(p)+prior_fun(p)),...
    size(X0,1),3e3, 29, 'BurnIn', 1e4);
toc

tic
[x3, p_x3,accept3] = dream(x2(:,:,end)',@(p)(likelihood_fun(p)+prior_fun(p)),...
    size(X0,1),1e3, 29, 'BurnIn', 1e4,'StepSize', 1.19);
toc

tic
obj_fun(x_mat(end,:))
toc
%% Concatenating chains
x_a=cat(3,x,x2);
logP=[p_x; p_x2];
x_mat=x(:,:)';
%% Plotting results
plotTrace(x_a(:,:,2e4:end),PI)
figure
plotTrace(x,PI)
[c,lags,ess]=eacorr(x);

plotMCMCDiagnostics(x,p_x,PI,'model', 'Kosinsky et al. (2018)');
plotBivariateMarginals_2(exp(x_mat(1:1:end,H.PopulationParams)),PI)
plotBivariateMarginals_2(exp(x_mat(7e5:1e3:end,[H.IndividualParams.EtaIndex H.IndividualParams.OmegaIndex])),PI)
plotBivariateMarginals_2(exp(x_mat(7e5:1e3:end,[H.IndividualParams.EtaIndex...
    H.IndividualParams.Index])),PI,'names', {PI.par([H.IndividualParams.EtaIndex...
    H.IndividualParams.Index]).name})
plotBivariateMarginals_2(exp(x_mat(7e5:1e3:end, H.SigmaParams)),PI,'names', {PI.par(H.SigmaParams).name})
plot(logP)
plot(lags,c)

%% Posterior predictions
simFun=@(x)getOutput(PI,@(p)sim(p,100,u,1:1:100),x,...
    @(p)getPhi2(p,H,length(u)), 4:6);
tic
PI=getPosteriorPredictions(exp(x_mat(1e5:1000:end,:)),PI,simFun,observables);
toc
PI=getCredibleIntervals(PI,observables, exp(x_mat(1e5:1000:end,:)),H);
plotPosteriorPredictions(PI,observables)
plotGroupPosteriorPredictions(PI,observables(1),exp(x_mat(1e5:100:end,:)),H)

%% Save results
save('/Users/migueltenorio/Documents/MATLAB/SimBiology/Kosinsky/output/PI_kpro/DREAM_MCMC_kpro.mat', 'x_a')
save('/Users/migueltenorio/Documents/MATLAB/SimBiology/Kosinsky/output/PI_kpro/DREAM_MCMC_logP_kpro.mat', 'logP')