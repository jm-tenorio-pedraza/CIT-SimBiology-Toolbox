function par_struct = getParamStruct2(simFun,H,n_sim,sigmaParams,sigmaNames,varargin)
if nargin<5
    error('getParamStruct:toofewinputs','getParamStruct requires atleast 5 inputs.')
end
inputs=inputParser;
inputs.addParameter('Sigma',3);
inputs.addParameter('startSigma',repelem(.1,length(H.SigmaParams),1));
inputs.parse(varargin{:});
inputs=inputs.Results;

% Cell-specific params
cell_indx=arrayfun(@(x)(x.name),H.CellParams,'UniformOutput',false);        % Extract field names of individually-varying parameters
try
param_cell=cellfun(@(x)simFun.Parameters.Value(ismember(simFun.Parameters.Name,x)),...
    cell_indx,'UniformOutput',false);                                      % Compare names of parameters and extract the corresponding value
param_cell=reshape(repelem(cell2mat(param_cell),1,length(H.CellParams(1).Index))',[],1);
catch
    param_cell =[];
end

% Individual params
indiv_indx=arrayfun(@(x)(x.name),H.IndividualParams,'UniformOutput',false); % Extract field names of individually-varying parameters
try
param_indiv=cellfun(@(x)simFun.Parameters.Value(ismember(simFun.Parameters.Name,x)),...
    indiv_indx,'UniformOutput',false);                                      % Compare names of parameters and extract the corresponding value
param_indiv=reshape(repelem(cell2mat(param_indiv),1,n_sim)',[],1);
catch
    param_indiv =[];
end
% Population params
param=[simFun.Parameters.Value(H.PopulationParams)];       % 


% concatenating all params

    % Min            Max                Start           % Mu             % Sigma
p=[ param*1e-1       param*1e1           param           param 
    param_cell*.01   param_cell*100      ones(size(param_cell))      ones(size(param_cell))
    param_indiv*.01  param_indiv*100     ones(size(param_indiv))     ones(size(param_indiv))     
    sigmaParams*.01  sigmaParams*100     inputs.startSigma     sigmaParams];

if length(inputs.Sigma)<2
    % Sigma
    sigma=[ repelem(inputs.Sigma,length(param),1);
        repelem(inputs.Sigma,length(param_indiv),1);
        repelem(inputs.Sigma,length(sigmaParams),1)];
else
    sigma=inputs.Sigma;
end

% Parameter names
 paramNames=[simFun.Parameters.Name(H.PopulationParams);...
        repelem({(H.CellParams(1:end).name)}',length(H.CellParams(1).Index),1);
        repelem({(H.IndividualParams(1:end).name)}',length(H.IndividualParams(1).Index),1);...
        sigmaNames];
% condition = ~isempty(H.CellParams(1).name)+ ~isempty(H.IndividualParams(1).name);
% switch condition
%     case 0
%         paramNames=[simFun.Parameters.Name(H.PopulationParams);...
%         repelem({(H.CellParams(1:end).name)}',length(H.CellParams(1).Index),1);
%         repelem({(H.IndividualParams(1:end).name)}',length(H.IndividualParams(1).Index),1);...
%         sigmaNames];
%     case 1
%     if ~isempty(H.CellParams(1).name)
%     
%     paramNames([H.CellParams.EtaIndex])=arrayfun(@(x)strjoin({'eta' x.name},'_'),...
%         H.CellParams','UniformOutput', false);
%     paramNames([H.IndividualParams.EtaIndex])=arrayfun(@(x)strjoin({'eta' x.name},'_'),...
%         H.IndividualParams','UniformOutput', false);
%     elseif ~isempty(H.IndividualParams(1).name)
%          
%      paramNames=[simFun.Parameters.Name(H.PopulationParams);...
%         repelem({(H.CellParams(1:end).name)}',length(H.CellParams(1).Index),1);
%         repelem({(H.IndividualParams(1:end).name)}',length(H.IndividualParams(1).Index),1);...
%         sigmaNames];
%     paramNames([H.CellParams.EtaIndex])=arrayfun(@(x)strjoin({'eta' x.name},'_'),...
%         H.CellParams','UniformOutput', false);
%     else
%     paramNames([H.IndividualParams.EtaIndex])=arrayfun(@(x)strjoin({'eta' x.name},'_'),...
%         H.IndividualParams','UniformOutput', false);
%     end
%     
% catch
%     disp('Parameter array will be built using population params only')
%     paramNames=[simFun.Parameters.Name(setdiff(H.PopulationParams, [H.IndividualParams.EtaIndex]));...
%     sigmaNames];
% end
% Parameter structure
par_struct=struct('name', paramNames,'minValue',num2cell(p(:,1)), 'maxValue',...
    num2cell(p(:,2)),'startValue', num2cell(p(:,3)),'mu_prior',  num2cell(p(:,4)),'sigma_prior',  num2cell(sigma));
return