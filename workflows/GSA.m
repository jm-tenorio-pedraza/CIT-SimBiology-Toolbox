%% GSA
time = 1:PI.tspan(end);
% inputs = [PI.par(PI.H.PopulationParams).finalValue];
inputs = [PI.par(PI.H.PopulationParams).MAP];
inputs = [repelem(inputs,size(PI.x_0,1),1) PI.x_0(:,1)];

if ~isempty(PI.H.CellParams(1).Index)
z = reshape([PI.par([PI.H.CellParams.Index]).finalValue],[],length(PI.H.CellParams));
inputs(:,[PI.H.CellParams.EtaIndex]) = inputs(:,[PI.H.CellParams.EtaIndex]).*(PI.H.CellIndx*z);

end
if ~isempty(PI.H.IndividualParams(1).Index)
    w = reshape([PI.par([PI.H.IndividualParams.Index]).finalValue],[],length(PI.H.IndividualParams));
 inputs(:,[PI.H.IndividualParams.EtaIndex]) = inputs(:,[PI.H.IndividualParams.EtaIndex]).*w;

end
 group = [PI.data(:).Name];
if ischar(group)
    group = {PI.data(:).Name};
end
% [group, indx] = unique(group, 'stable');
% u_subset = PI.u(indx,:);
% inputs = inputs(indx,:);
% Get sensitivity matrix
sensmatrix = getSensitivities(inputs, PI,@(p)sim(p,PI.tspan(end),PI.u,1:1:PI.tspan),...
    parameters,observables,time,'initialValue', true,'uniqueGroups',false);

% Get SVD
[U,S,V]=svd(sensmatrix,'econ');

%% plot singular values
plotSensitivities(S)
%% Get SHM
F = max(abs(V),[],2)'.*diag(S)'.*abs(U);
s=shmPlot2(F,group,time,observables,'tau',0.01);

%% Get PSS
pcs = V*S;
pcs = (pcs/max(max(abs(pcs))));
PI.pcs = pcs;

figure
pc = plotPSS(PI.pcs,6,PI.paramNames(PI.H.PopulationParams),'threshold',-1,'newFig', false);
%% Parameters
parameters_hat = cat(1,pc(:).p_hat);
parameters_hat = unique(parameters_hat,'stable');

%% Global PRCC SA
PI = globalSA2(sim,PI,observables,'samples',[exp(postSamples(:,[PI.H.PopulationParams])), ...
    repelem(PI.x_0,size(postSamples,1),1)],'nsamples',1e4,'time',1:1:PI.tspan(end),...
    'sigma', 1,'inputs', exp(finalValues(PI.H.PopulationParams)),'variation', 0.5);
PI = globalSA2(sim,PI,observables,'samples',[],'nsamples',1e2,'time',1:1:PI.tspan(end),...
    'sigma', 1,'inputs', exp(finalValues(PI.H.PopulationParams)),'variation', 0.5);

paramRanking = plotPRCC(PI,PI.paramNames(PI.H.PopulationParams),PI.observablesPlot,'time',...
    1:1:PI.tspan(end),'output', 'mean','kpi', 'sum');
parameters_hat1 =unique(paramRanking{1:8,:});

%% Joint results
unique([parameters_hat; parameters_hat1])
%% Save results to cd
save(strjoin({cd 'parameters_hat_2.mat'},'/'), 'parameters_hat')
sensitivity = false;

