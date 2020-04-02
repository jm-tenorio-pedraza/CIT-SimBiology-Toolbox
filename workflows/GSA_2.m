%% GSA
time = 1:1:PI.tspan(end);
inputs = [PI.par(PI.H.PopulationParams).finalValue];

variant = [[PI.par([PI.H.CellParams.EtaIndex]).finalValue];reshape([PI.par([PI.H.CellParams.Index]).finalValue],[],length(PI.H.CellParams))];
eta = [[PI.par([PI.H.IndividualParams.EtaIndex]).finalValue];reshape([PI.par([PI.H.IndividualParams.Index]).finalValue],[],length(PI.H.IndividualParams))];

inputs = [repelem(inputs,size(PI.x_0,1),1) PI.x_0(:,1)];
inputs(:,[PI.H.CellParams.EtaIndex]) = inputs(:,[PI.H.CellParams.EtaIndex]).*(PI.H.CellIndx*variant);
inputs(:,[PI.H.IndividualParams.EtaIndex]) = inputs(:,[PI.H.IndividualParams.EtaIndex]).*eta;
group = [PI.data(:).Group];

% Get sensitivity matrix
sensmatrix = getSensitivities(inputs, PI,@(p)sim(p,PI.tspan(end),PI.u,time),...
    parameters, observables,time,'initialValue', true,'uniqueGroups',false);

% Get SVD
[U,S,V]=svd(sensmatrix,'econ');

%% plot singular values
plotSensitivities(S)
%% Get SHM
F = max(abs(V),[],2)'.*diag(S)'.*abs(U);
s=shmPlot2(F,group,time, observables,'tau',0.1);

%% Get PSS
pcs = V*S;
pcs = (pcs/max(max(abs(pcs))));
pc = plotPSS(pcs,4,paramNames(PI.H.PopulationParams),'threshold',-1);
%% Parameters
parameters_hat = cat(1,pc(:).p_hat);
parameters_hat = unique(parameters_hat,'stable');

%% Global PRCC SA
PI = globalSA(sim,PI,observables,'nsamples',1000,'time',1:3:PI.tspan(end),...
    'sigma', 1,'inputs', inputs(1,:));
plotPRCC(PI,parameters,observables(1),'time',1:3:PI.tspan(end))
%% Save results to cd
save(strjoin({cd 'parameters_hat_2.mat'},'/'), 'parameters_hat')
sensitivity = false;