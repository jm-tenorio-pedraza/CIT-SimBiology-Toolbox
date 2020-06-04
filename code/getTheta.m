function [theta, H] = getTheta(Meta,params, N_pop, N_cell,N_indiv)

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
N_tot = N_pop*N_cell*N_indiv;
theta = nan(N_tot, length(params));
for i=1:length(Meta)
    PI = Meta(i).Struct.PI;
    H.par(i).popIndx = PI.H.PopulationParams;                               % Indexes of the pop params
    H.par(i).etaIndx = [PI.H.CellParams.EtaIndex...                         % Indexes of the mean params of IEV+IIV params
        PI.H.IndividualParams.EtaIndex];
    N_cell_i  = length(PI.H.CellParams(1).Index);
    N_indiv_i = length(PI.H.IndividualParams(1).Index);
    
    p_cell_i = length([PI.H.CellParams]);
    p_indiv_i = length([PI.H.IndividualParams]);
    try
            H.par(i).omegaIndx= [PI.H.CellParams.OmegaIndex...              % Indexes of sigma params of IEV+IIV params
                PI.H.IndividualParams.OmegaIndex];
    catch
            H.par(i).omegaIndx= [PI.H.CellParams.OmegaIndex];

    end
   
    H.par(i).sigmaIndx = setdiff( PI.H.SigmaParams ,...
        H.par(i).omegaIndx );                                               % Identify indexes of error variance params
    H.par(i).name = {PI.par([H.par(i).popIndx ...
        H.par(i).sigmaIndx]).name};                                         % Set names of Pop and error var params
    
    H.par(i).popSampleIndx = randsample(size(PI.postSamples,1)...
        ,N_pop, true);                                                      % Random sample of population params
    H.par(i).postSamples = repelem(PI.postSamples(...
        H.par(i).popSampleIndx, H.par(i).popIndx),N_cell*N_indiv,1);               % Repeat the mean population params N_indiv times
    H.par(i).omegaSamples = repelem(exp(PI.postSamples(...                  % Repeat the omega population params N_indiv times
        H.par(i).popSampleIndx, H.par(i).omegaIndx)),N_cell*N_indiv,1);
    H.par(i).sigmaSamples = repelem(exp(PI.postSamples(...                  % Repeat the sigma population params N_indiv times
        H.par(i).popSampleIndx, H.par(i).sigmaIndx)),N_cell*N_indiv,1);
    try
    H.par(i).indivPostSamples = repelem(PI.postSamples( H.par(i).popSampleIndx,...
        [PI.H.IndividualParams(:).Index]), N_cell*N_indiv,1);
    catch
    end
    try
        H.par(i).cellPostSamples = repelem(PI.postSamples(H.par(i).popSampleIndx,...
            [PI.H.CellParams(:).Index]), N_cell*N_indiv,1);
    catch
    end
    %% Sample cell and individual parameters
    if N_cell_i~=0
        cellSampleIndx  = (randsample(1:size(PI.H.CellIndx,2),...
                N_tot,1));
        indivSampleIndx = nan(N_tot, 1);
        cellSampleIndxMat = nan(N_tot, N_cell_i);
        indivSampleIndxMat = nan(N_tot, N_indiv_i);

        for j=1:N_cell_i
            cellSampleIndx_j = ismember(cellSampleIndx,j);
            cellSampleIndxMat(:,j) = cellSampleIndx_j; 
            indivSampleIndx( cellSampleIndx_j,:) = randsample(find(PI.H.CellIndx(:,j)),...
                sum(cellSampleIndx_j),1);
            for k=1:N_indiv_i
                
            indivSampleIndxMat(cellSampleIndx_j,k) = ...
                ismember(indivSampleIndx( cellSampleIndx_j,:), k);
            end
        end
    elseif N_indiv_i ~=0
        indivSampleIndxMat = nan(N_tot, N_indiv_i);

        indivSampleIndx = randsample(1:length([PI.H.IndividualParams(:).Index]),N_tot,1);
        for k=1:N_indiv_i
                
            indivSampleIndxMat(:,k) = ...
                ismember(indivSampleIndx( :,:), k);
        end
    end
    
    try
        indivPostSamples = nan(N_tot, length(PI.H.IndividualParams));
                indx = 1:length(PI.H.IndividualParams(1).Index);

        for j=1:length(PI.H.IndividualParams)
            indivPostSamples=H.par(i).indivPostSamples(:,indx);

            indivPostSamples(:,j) =sum(indivPostSamples.*cellSampleIndxMat,2);
             indx = indx+length(PI.H.IndividualParams(1).Index);

        end
    catch
         indivPostSamples = [];

    end
    
    try
        cellPostSamples = nan(N_tot, length(PI.H.CellParams));
        indx = 1:length(PI.H.CellParams(1).Index);

        for j=1:length(PI.H.CellParams)
            cellPostSamples_j=H.par(i).cellPostSamples(:,indx);
            cellPostSamples(:,j) = sum(cellPostSamples_j.*cellSampleIndxMat,2);
            indx = indx+length(PI.H.CellParams(1).Index);
        end
    catch
                cellPostSamples = [];

    end
  
    H.par(i).postSamples(:,H.par(i).etaIndx) =  ...
         H.par(i).postSamples(:,H.par(i).etaIndx) +[cellPostSamples indivPostSamples];
    H.par(i).variables = Meta(i).Struct.PI.observablesPlot;
    
    [theta_indx, theta_order] = ismember(params,  H.par(i).name);           % Compare parameters against the required ones
    theta_order = theta_order(theta_order~=0);                              % Identify non-zero indexes
    theta(:,theta_indx) = H.par(i).postSamples(:, theta_order);             % Re-arrange samples to match desired output
    

end