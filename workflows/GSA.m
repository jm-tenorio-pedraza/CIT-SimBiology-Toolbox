%% GSA
time = 1:100;

% inputs = sim.Parameters.Value';
inputs = [PI.par(H.PopulationParams).finalValue];

% Get sensitivity matrix
sensmatrix = getSensitivities(inputs, PI,@(p)sim(p,100,u,1:1:100),parameters, observables,time);

% Get SVD
[U,S,V]=svd(sensmatrix,'econ');

% plot singular values
s_ij = log10(diag(S)/max(diag(S)));
plot(s_ij, '-d')
hold on
plot(0:length(s_ij),ones(1,1+length(s_ij))*(-1),'-k')
title('Singular values of sensitivity matrix')
ylabel('Log_{10} \sigma_i / max \sigma_i')
%% Get SHM
F = max(abs(V),[],2)'.*diag(S)'.*abs(U);
s=shmPlot2(F,groups_subset,time, observables,'tau',0.01);

%% Get PSS
pcs = V*S;
pcs = (pcs/max(max(abs(pcs))));
pc = plotPSS(pcs,5,parameters,'threshold',-1);
%% Parameters

parameters_hat = cat(1,pc(:).p_hat);
parameters_hat = unique(parameters_hat,'stable');
%% Save results to cd
save(strjoin({cd 'parameters_hat.mat'},'/'))
sensitivity = false;