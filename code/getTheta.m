function [theta, H] = getTheta(Meta,params, N_pop, N_indiv)

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
H.par=[];
theta = nan(N_pop*N_indiv, length(params));
for i=1:length(Meta)
    H.par(i).popIndx = Meta(i).Struct.PI.H.PopulationParams;
    H.par(i).etaIndx = [Meta(i).Struct.PI.H.CellParams.EtaIndex Meta(i).Struct.PI.H.IndividualParams.EtaIndex];
    try
            H.par(i).omegaIndx= [Meta(i).Struct.PI.H.CellParams.OmegaIndex Meta(i).Struct.PI.H.IndividualParams.OmegaIndex];
    catch
            H.par(i).omegaIndx= [Meta(i).Struct.PI.H.CellParams.OmegaIndex];

    end
    H.par(i).sigmaIndx = setdiff(Meta(i).Struct.PI.H.SigmaParams,H.par(i).omegaIndx);
    H.par(i).name = {Meta(i).Struct.PI.par([H.par(i).popIndx H.par(i).sigmaIndx]).name};
    
    H.par(i).postSampleIndx = randsample(size(Meta(i).Struct.PI.postSamples,1),N_pop, true);
    H.par(i).postSamples = repelem(Meta(i).Struct.PI.postSamples(H.par(i).postSampleIndx, H.par(i).popIndx),N_indiv,1);
    H.par(i).omegaSamples = repelem(exp(Meta(i).Struct.PI.postSamples(H.par(i).postSampleIndx, H.par(i).omegaIndx)),N_indiv,1);
    H.par(i).sigmaSamples = repelem(exp(Meta(i).Struct.PI.postSamples(H.par(i).postSampleIndx, H.par(i).sigmaIndx)),N_indiv,1);
    
    eta = randn(N_pop*N_indiv, length(H.par(i).etaIndx)).*(H.par(i).omegaSamples);
    H.par(i).postSamples(:,H.par(i).etaIndx) =  H.par(i).postSamples(:,H.par(i).etaIndx)+eta;
    H.par(i).variables = Meta(i).Struct.PI.observablesPlot;
    
    [theta_indx, theta_order] = ismember(params,  H.par(i).name);
    theta_order = theta_order(theta_order~=0);
    theta(:,theta_indx) = H.par(i).postSamples(:, theta_order);
    

end