

%% External validation
PI_Morisada=getPIData('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Morisada.mat',...
    stateVar,{'Control' 'antiPD1'},observables,'zeroAction', 'omit','mergePhenotypes', true);
groups_Morisada = {'MOC1_antiPDL1';'MOC1_Control'};
[PI_Morisada.data(1:end).Group] = groups_Morisada{:,:};



% Posterior predictions
simFun=@(x)getOutput(PI,@(p)sim(p,PI_Morisada.tspan(end),u,PI_Morisada.tspan),x,...
    @(p)getPhi2(p,PI_Morisada.H,length(PI_Morisada.u),'initialValue',PI_Morisada.x_0),...
    PI_Morisada.normIndx, PI_Morisada.H);
tic
PI=getPosteriorPredictions(exp(postSamples),PI,simFun,observables);
toc
PI=getCredibleIntervals(PI,observables, exp(postSamples),H);
plotPosteriorPredictions(PI,observables)

