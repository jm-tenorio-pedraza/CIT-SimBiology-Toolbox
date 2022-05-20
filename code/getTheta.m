function [theta, H] = getTheta(Meta,params, N_pop, N_cell,N_indiv,N_resp, varargin)
inputs = inputParser;
inputs.addParameter('covStructure', 'dependent')
inputs.parse(varargin{:})
inputs = inputs.Results;
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
N_tot = N_pop*N_cell*N_indiv*N_resp;
theta = nan(N_tot, length(params));
for i=1:length(Meta)
    PI = Meta(i).Struct.PI;
    H.par(i).popIndx = PI.H.PopulationParams;                               % Indexes of the pop params
    H.par(i).etaIndx = [PI.H.CellParams.EtaIndex...                         % Indexes of the mean params of IEV+IIV params
        PI.H.IndividualParams.EtaIndex PI.H.RespParams.EtaIndex];
    N_cell_i  = length(PI.H.CellParams(1).Index);
    N_indiv_i = length(PI.H.IndividualParams(1).Index);
    N_resp_i = length(PI.H.RespParams(1).Index);
    try
        H.par(i).omegaIndx= [PI.H.CellParams.OmegaIndex...              % Indexes of sigma params of IEV+IIV params
            PI.H.IndividualParams.OmegaIndex PI.H.RespParams.OmegaIndex];
    catch
        H.par(i).omegaIndx= [PI.H.CellParams.OmegaIndex];

    end
   
    H.par(i).sigmaIndx = setdiff( PI.H.SigmaParams ,...
        H.par(i).omegaIndx );                                               % Identify indexes of error variance params
    H.par(i).name = {PI.par([H.par(i).popIndx ...
        H.par(i).sigmaIndx]).name};                                         % Set names of Pop and error var params
    
    H.par(i).popSampleIndx = randsample(size(PI.postSamples,1),N_pop, true); % Random sample of population params
    H.par(i).postSamples = repelem(PI.postSamples(H.par(i).popSampleIndx,...
        H.par(i).popIndx),N_cell*N_indiv*N_resp,1);               % Repeat the mean population params N_indiv times
    H.par(i).omegaSamples = repelem(exp(PI.postSamples(...                  % Repeat the omega population params N_indiv times
        H.par(i).popSampleIndx, H.par(i).omegaIndx)),N_cell*N_indiv*N_resp,1);
    H.par(i).sigmaSamples = repelem(exp(PI.postSamples(...                  % Repeat the sigma population params N_indiv times
        H.par(i).popSampleIndx, H.par(i).sigmaIndx)),N_cell*N_indiv*N_resp,1);
    
    try
    H.par(i).indivPostSamples = repelem(PI.postSamples( H.par(i).popSampleIndx,...
        [PI.H.IndividualParams(:).Index]), N_cell*N_indiv*N_resp,1);
    catch
    end
    try
        H.par(i).cellPostSamples = repelem(PI.postSamples(H.par(i).popSampleIndx,...
            [PI.H.CellParams(:).Index]), N_cell*N_indiv*N_resp,1);
    catch
    end
    try
        H.par(i).respPostSamples = repelem(PI.postSamples(H.par(i).popSampleIndx,...
            [PI.H.RespParams(:).Index]), N_cell*N_indiv*N_resp,1);
    catch
    end
    %% Sample cell and individual parameters
    if strcmp(inputs.covStructure,'dependent')
    % If cell-specific parameters are present, sample them
    if N_cell_i~=0
        cellSampleIndx      = (randsample(1:size(PI.H.CellIndx,2),N_tot,1)); % Sample the indexes of the cell-specific estimates N_tot
        indivSampleIndx     = nan(N_tot, 1); 
        cellSampleIndxMat   = nan(N_tot, N_cell_i); % Matrix of logical indexes with each row representing a sample and each column the cell line
        indivSampleIndxMat  = nan(N_tot, N_indiv_i);
        respSampleIndxMat   = nan(N_tot, N_resp_i);
        for j=1:N_cell_i % For each unique cell type
            cellSampleIndx_j = ismember(cellSampleIndx,j); % Identify the samples corresponding to the j unique cell type
            cellSampleIndxMat(:,j) = cellSampleIndx_j;    
            indivSampleIndx(cellSampleIndx_j,:) = ... % For each individual select a random sample of the cell type when there are more than 1 estimates for 1 cell type
                randsample(find(PI.H.CellIndx(:,j)),...
                sum(cellSampleIndx_j),1);
            for k=1:N_indiv_i
                indivSampleIndxMat(cellSampleIndx_j,k) = ...
                    ismember(indivSampleIndx( cellSampleIndx_j,:), k); % Transform from matrix of numerical indexes to logical indexes
            end
        end
    else
        if N_indiv_i ~=0 % If there is no cell-specific estimates proceed to the individual estimates if these are present and do the same
            indivSampleIndxMat = nan(N_tot, N_indiv_i);
            indivSampleIndx = randsample(1:length([PI.H.IndividualParams(1).Index]),...
                N_tot,1);
            for k=1:N_indiv_i
                indivSampleIndxMat(:,k) = ...
                    ismember(indivSampleIndx( :,:), k);
            end
        elseif N_resp_i ~=0
            respSampleIndxMat = nan(N_tot, N_indiv_i);
            respSampleIndx = randsample(1:length([PI.H.RespParams(1).Index]),...
                N_tot,1);
            for k=1:N_indiv_i
                respSampleIndxMat(:,k) = ...
                    ismember(respSampleIndx( :,:), k);
            end
        end
    end
    
    if ~isempty(PI.H.IndividualParams(1).name) % If individual estimates are present
        indivPostSamples = nan(N_tot, length(PI.H.IndividualParams)); % Generate matrix of samples with num of columns equal to the number of invididually-estimated parameters
        indx = 1:length(PI.H.IndividualParams(1).Index); % Number of individual estimates in each parameter
        for j=1:length(PI.H.IndividualParams) % For each parameter
            indivPostSamples_j=H.par(i).indivPostSamples(:,indx); % Select the corresponding parameter samples 
            indivPostSamples(:,j) =sum(indivPostSamples_j.*indivSampleIndxMat,2); % In each row only 1 column != 0, sum just gets rid of them
            indx = indx+length(PI.H.IndividualParams(1).Index); % Update the index parameter after each parameter is sampled so the samples align with the actual parameter
        end
    else
         indivPostSamples = [];
    end
    
    if ~isempty(PI.H.CellParams(1).name)
        cellPostSamples = nan(N_tot, length(PI.H.CellParams));
        indx = 1:length(PI.H.CellParams(1).Index);
        for j=1:length(PI.H.CellParams)
            cellPostSamples_j=H.par(i).cellPostSamples(:,indx);
            cellPostSamples(:,j) = sum(cellPostSamples_j.*cellSampleIndxMat,2);
            indx = indx+length(PI.H.CellParams(1).Index);
        end
    else
        cellPostSamples = [];
    end
    
      if ~isempty(PI.H.RespParams(1).name)
        respPostSamples = nan(N_tot, length(PI.H.RespParams));
        indx = 1:length(PI.H.RespParams(1).Index);
        for j=1:length(PI.H.RespParams)
            respPostSamples_j=H.par(i).respPostSamples(:,indx);
            respPostSamples(:,j) = sum(respPostSamples_j.*respSampleIndxMat,2);
            indx = indx+length(PI.H.RespParams(1).Index);
        end
    else
        respPostSamples = [];
      end
    try
    H.par(i).postSamples(:,H.par(i).etaIndx) =  ...
        H.par(i).postSamples(:,H.par(i).etaIndx) +[cellPostSamples indivPostSamples respPostSamples];
    catch
    end
    else % If the lower hierarchy parameters are assumed to be independent from the ones in the higher level
        Z = randn(N_tot, length(H.par(i).etaIndx));
        try
        H.par(i).postSamples(:,H.par(i).etaIndx) = ...
            H.par(i).postSamples(:,H.par(i).etaIndx) + Z.*H.par(i).omegaSamples;
        catch 
        end
    end
     H.par(i).variables = Meta(i).Struct.PI.observablesPlot;
    
    [theta_indx, theta_order] = ismember(params,  H.par(i).name);           % Compare parameters against the required ones
    theta_order = theta_order(theta_order~=0);                              % Identify non-zero indexes
    theta(:,theta_indx) = H.par(i).postSamples(:, theta_order);             % Re-arrange samples to match desired output
end