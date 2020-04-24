function [theta, H] = getTheta(PI,params, N_pop, N_indiv)

% Inputs:
%   - PI: Structure of PIs from the different optimizations to be used
%   - N_pop: Number of samples for the population-distributed parameters
%   - N_indiv: Number of samples for the individually-distributed
%   parameters per each population parameter
% Outputs:
%   - H: structure that contains the joint posterior samples, indexes from
%   all PIs
%   - theta: subset of posterior dimensions from the joint posterior to
%   consider
H.popIndx = [];
H.sigmaIndx = [];
H.cellIndx = [];
H.indivIndx = [];
theta = nan(N_pop*N_indiv, length(params));
for i=1:length(PI)
    H.par(i).popIndx = PI(i).PI.H.PopulationParams;
    H.par(i).etaIndx = [PI(i).PI.H.CellParams.EtaIndex PI(i).PI.H.IndivParams.EtaIndex];
    H.par(i).omegaIndx= [PI(i).PI.H.CellParams.OmegaIndx PI(i).PI.H.IndivParams.OmegaIndx];
    H.par(i).sigmaIndx = setdiff(PI(i).PI.H.SigmaParams,H.par(i).omegaIndx);
    H.par(i).name(i) = {PI(i).PI.par([H.par(i).popIndx H.par(i).sigmaIndx]).name};
    
    H.par(i).postSampleIndx = randsample(size(PI(i).PI.postSamples,1),N_pop, true);
    H.par(i).postSamples = PI(i).PI.postsamples(H.par(i).postSampleIndx, H.par(i).popIndx);
    H.par(i).omegaSamples = exp(PI(i).PI.postsamples(H.par(i).postSampleIndx, H.par(i).omegaIndx));
    H.par(i).sigmaSamples = exp(PI(i).PI.postsamples(H.par(i).postSampleIndx, H.par(i).sigmaIndx));
    
    eta = randn(N_pop*N_indiv, length(H.par(i).etaIndx).*repmat(H.par(i).omegaSamples),N_indiv,1);
    H.par(i).postSamples(:,H.par(i).etaIndx) = eta;
    H.data(i).variables = PI(i).PI.observables;
    
    [theta_indx, theta_order] = strcmp(params,  H.par(i).name(i));
    theta(:,theta_indx) = H.par(i).postSamples(:, theta_order);
    

end