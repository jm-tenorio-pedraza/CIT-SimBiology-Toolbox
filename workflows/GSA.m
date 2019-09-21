%% Global sensitivity analysis

% General setup for sensitivity evaluation
[sim,u]=initializePI(model,parameters,observables,PI,dose,'CT26');

%% Hierarchical structure
H.PopulationParams=1:length(parameters);
H.IndividualParams=struct('name', [],'Index', [], 'EtaIndex', [], 'OmegaIndex', []);

% Generating PI
PI.par = getParamStruct2(sim,H,1,[],[],'Sigma', repelem(1,length(parameters),1));

% Parameter values to use
inputs = [PI.par(:).startValue];
%time = {PI.data(:).dataTime};time = cell2mat(cellfun(@(x)x',time,'UniformOutput',false));
%time = unique(time);
time = 1:100;
PI=getOutput(PI,@(p)sim(p,100,u,time),(inputs),...
    @(p)getPhi2(p,H,size(PI.data,1)), [],time);
sensmatrix = getSensitivities(inputs, PI,@(x)sim(x,100, u, time), parameters, observables,time);
%% SVD of sensitivity matrix
[U, S, V] = svd(sensmatrix);
sigma = diag(S)/max(diag(S));
S_ij = (V*S(1:size(V,1),:));
pcs = S_ij./max(max(abs(S_ij)));
%% Singular values

figure;
s=plot(log10(sigma), '--d');
hold on
plot(s.XData,-1,'.-k')
% s.BaseValue = -3;
% set(gca, 'YScale', 'log', 'YLim', [1e-5 1]);
ylabel('Log10(\sigma_i/max(\sigma_i))');xlabel('Index i')
title('Singular values of sensitivity matrix')
%% SHM
F = abs(U(:,1:size(V,1))).*diag(S)'.*(max(abs(V'),[],2)');
group = [PI.data(:).Group];

F_t = shmPlot2(F,group,time,observables,'tau', 0.01);

%% PSS
par_hat = plotPSS(pcs,2,parameters,'threshold', -2);
%% Save output
save('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/Kosinsky/phat_Control.mat','par_hat')
save('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/CIM/phat_TV.mat','par_hat')
%% Compare parameter sensitivities between SAs
phat_Control = load('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/Kosinsky/phat_Control.mat','par_hat');
phat_TV = load('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/CIM/phat_TV.mat');
phat_TV = {phat_TV.par_hat(:).p_hat}';
phat_TV = cat(1,phat_TV{:,:});
phat_Control = {phat_Control.par_hat(:).p_hat}';
phat_Control = cat(1,phat_Control{:,:});
try
parameters_hat=unique([phat_Control; phat_TV]);
catch
    parameters_hat = phat_Control;
end
%%
save('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/Kosinsky/parameters_hat.mat', 'parameters_hat')
