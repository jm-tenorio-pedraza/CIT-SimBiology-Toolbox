%% Global sensitivity analysis

%% General setup for sensitivity evaluation
[sim,u]=initializePI(model,parameters,observables,PI,doses,'CT26');

% Hierarchical structure
H.PopulationParams=1:length(parameters);
H.IndividualParams=struct('name', [],'Index', [], 'EtaIndex', [], 'OmegaIndex', []);

% Generating PI
PI.par = getParamStruct2(sim,H,1,[],[],'Sigma', repelem(1,length(parameters),1));

% Parameter values to use
inputs = [PI.par(:).startValue];
time = {PI.data(:).dataTime};time = cell2mat(cellfun(@(x)x',time,'UniformOutput',false));
time = unique(time);
PI=getOutput(PI,@(p)sim(p,100,u,time),(inputs),...
    @(p)getPhi2(p,H,size(PI.data,1)), length(observables)-1:length(observables),time);
sensmatrix = getSensitivities(inputs, PI,@(x)sim(x,100, u, time), parameters, observables,time);
%% SVD of sensitivity matrix
[U, S, V] = svd(sensmatrix);
sigma = diag(S)/max(diag(S));
W = (V'*S(1:size(V,1),:));
pcs = W./max(max(W));

figure;
plot(sigma, '--d')
set(gca, 'YScale', 'log', 'YLim', [1e-5 1]);ylabel('Log(\sigma_i)');xlabel('Index i')
title('Singular values of sensitivity matrix')
%% SHM
D = U(:,1:size(V,1))*V'*S(1:size(V,1),:);
F = abs(U(:,1:size(V,1))).*diag(S)'.*(max(abs(V'),[],2)');
group = [PI.data(:).Group];

F_t = shmPlot2(F,group,time,observables,'tau', 0.01);

%% PSS
par_hat = plotPSS(pcs,3,parameters,'threshold', -1);

