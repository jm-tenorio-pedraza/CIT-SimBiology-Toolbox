%% GSA
time = 1:1:PI.tspan(end);
inputs = sim.Parameters.Value(1:end-1)';
% inputs = [PI.par(H.PopulationParams).finalValue];
inputs = [repelem(inputs,size(PI.x_0,1),1) PI.x_0(:,1)];
[group, indx] = unique([PI.data(:).Group], 'stable');
u_subset = PI.u(indx,:);
inputs = inputs(indx,:);
% Get sensitivity matrix
sensmatrix = getSensitivities(inputs, PI,@(p)sim(p,PI.tspan(end),u_subset,1:1:PI.tspan),...
    parameters, observables,time,'initialValue', true,'uniqueGroups',true);

% Get SVD
[U,S,V]=svd(sensmatrix,'econ');

%% plot singular values
s_ij = log10(diag(S)/max(diag(S)));
plot(s_ij, '-d')
hold on
plot(0:length(s_ij),ones(1,1+length(s_ij))*(-1),'-k')
title('Singular values of sensitivity matrix')
ylabel('Log_{10} \sigma_i / max \sigma_i')
%% Get SHM
F = max(abs(V),[],2)'.*diag(S)'.*abs(U);
s=shmPlot2(F,group,time, observables,'tau',0.1);

%% Get PSS
pcs = V*S;
pcs = (pcs/max(max(abs(pcs))));
pc = plotPSS(pcs,6,parameters(1:end-1),'threshold',-1);
%% Parameters
parameters_hat = cat(1,pc(:).p_hat);
parameters_hat = unique(parameters_hat,'stable');


%% Global PRCC SA
PI = globalSA(sim,PI,observables,'nsamples',50,'time',1:3:PI.tspan(end),'sigma', 1);
plotPRCC(PI,parameters,observables(1),'time',1:3:PI.tspan(end))
%% Save results to cd
save(strjoin({cd 'parameters_hat_2.mat'},'/'), 'parameters_hat')
sensitivity = false;