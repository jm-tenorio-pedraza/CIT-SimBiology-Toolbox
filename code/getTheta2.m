function [theta, H] = getTheta2(Meta,params, N_tot, varargin)
inputs = inputParser;
inputs.addParameter('covStructure', 'dependent')
inputs.addParameter('fixedParams', 'uniform')

inputs.parse(varargin{:})
inputs = inputs.Results;
% Inputs:
%   - PI: Structure of PIs from the different optimizations to be used
%   - N_pop: Number of samples for the population-distributed parameters
%   - N_indiv: Number of samples for the individually-distributed
%           parameters per each population parameter
% Outputs:
%   - H: structure that contains the joint posterior samples, indexes from
%   all PIs
%   - theta: subset of posterior dimensions from the joint posterior to
%   consider
H=[];
% Number of PIs
n_pi = length(Meta);
% Allocate space for the posterior samples
theta = nan(N_tot, length(params));
for i=1:length(Meta)

    PI = Meta(i).Struct.PI;
    H.par(i).popIndx = PI.H.PopulationParams;                               % Indexes of the pop params
    H.par(i).etaIndx = [PI.H.CellParams.EtaIndex...                         % Indexes of the mean params of IEV+IIV params
        PI.H.IndividualParams.EtaIndex PI.H.RespParams.EtaIndex];
    H.par(i).omegaIndx= [PI.H.CellParams.OmegaIndex...              % Indexes of sigma params of IEV+IIV params
            PI.H.IndividualParams.OmegaIndex PI.H.RespParams.OmegaIndex];
    H.par(i).sigmaIndx = setdiff( PI.H.SigmaParams ,...
        H.par(i).omegaIndx );                                               % Identify indexes of error variance params
    H.par(i).name = {PI.par([H.par(i).popIndx ...
        H.par(i).sigmaIndx]).name};                                         % Set names of Pop and error var params
    
    H.par(i).popSampleIndx = randsample(size(PI.postSamples,1),N_tot, true); % Random sample of population params
    H.par(i).postSamples = (PI.postSamples(H.par(i).popSampleIndx,...
        H.par(i).popIndx));                                                     % Repeat the mean population params N_indiv times
    H.par(i).omegaSamples = (exp(PI.postSamples(...                  % Repeat the omega population params N_indiv times
        H.par(i).popSampleIndx, H.par(i).omegaIndx)));
    H.par(i).sigmaSamples = (exp(PI.postSamples(...                  % Repeat the sigma population params N_indiv times
        H.par(i).popSampleIndx, H.par(i).sigmaIndx)));
    postSamples=PI.postSamples(H.par(i).popSampleIndx,:);
    try
    H.par(i).indivPostSamples = (PI.postSamples( H.par(i).popSampleIndx,...
        [PI.H.IndividualParams(:).Index]));
    catch
    end
    try
        H.par(i).cellPostSamples = (PI.postSamples(H.par(i).popSampleIndx,...
            [PI.H.CellParams(:).Index]));
    catch
    end
    try
        H.par(i).respPostSamples = (PI.postSamples(H.par(i).popSampleIndx,...
            [PI.H.RespParams(:).Index]));
    catch
    end
    %% Sample cell and individual parameters
    if strcmp(inputs.covStructure,'dependent')
    
        if ~isempty(PI.H.IndividualParams(1).name) % If individual estimates are present
            n_indiv=length(PI.H.IndividualParams(1).Index);
            indivIndx=randsample(n_indiv,N_tot,true);
            indivIndx=cellfun(@(x)ismember(indivIndx,x), num2cell(1:n_indiv),'UniformOutput',false);
            indivIndx=cell2mat(indivIndx);
            indivPostSamples=arrayfun(@(x)sum(indivIndx.*postSamples(:,x.Index),2),PI.H.IndividualParams,'UniformOutput',false);
        else
            indivPostSamples = [];
        end
        
        if ~isempty(PI.H.CellParams(1).name)
            n_indiv=length(PI.H.CellParams(1).Index);
            indivIndx=randsample(n_indiv,N_tot,true);
            indivIndx=cellfun(@(x)ismember(indivIndx,x), num2cell(1:n_indiv),'UniformOutput',false);
            indivIndx=cell2mat(indivIndx);
            cellPostSamples=arrayfun(@(x)sum(indivIndx.*postSamples(:,x.Index),2),PI.H.CellParams,'UniformOutput',false);
        else
            cellPostSamples = [];
        end
        
        if ~isempty(PI.H.RespParams(1).name)
            n_indiv=length(PI.H.RespParams(1).Index);
            indivIndx=randsample(n_indiv,N_tot,true);
            indivIndx=cellfun(@(x)ismember(indivIndx,x), num2cell(1:n_indiv),'UniformOutput',false);
            indivIndx=cell2mat(indivIndx);
            respPostSamples=arrayfun(@(x)sum(indivIndx.*postSamples(:,x.Index),2),PI.H.RespParams,'UniformOutput',false);
            respPostSamples=respPostSamples{:,:};
        else
            respPostSamples = [];
        end
        %% Sample
        try
            H.par(i).postSamples(:,H.par(i).etaIndx) =  ...
                H.par(i).postSamples(:,H.par(i).etaIndx) + [cellPostSamples indivPostSamples respPostSamples];
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
    if strcmp(inputs.fixedParams,'uniform')
    % Sample uniformly for those parameters that were not estimated
    lb=Meta(i).Values(~theta_indx)*.9;
    ub=Meta(i).Values(~theta_indx)*1.1;
    
    theta2=rand(N_tot,sum(~theta_indx)).*(ub-lb)'+lb';
    elseif strcmp(inputs.fixedParams,'fixed')
        theta2 = repmat(Meta(i).Values(~theta_indx)',N_tot,1);
    end
    theta(:,~theta_indx)=log(theta2);
    H.par(i).theta=theta;
end
theta={H.par(1:n_pi).theta}';
theta=cell2mat(theta);
end