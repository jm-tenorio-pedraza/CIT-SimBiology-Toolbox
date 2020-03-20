function par_struct = getParamStruct3(simFun,H,n_sim,sigmaParams,sigmaNames,varargin)
if nargin<5
    error('getParamStruct:toofewinputs','getParamStruct requires atleast 5 inputs.')
end
inputs=inputParser;
inputs.addParameter('Sigma',3);
inputs.addParameter('startSigma',[repelem(1,length([H.CellParams.OmegaIndex]),1); ...
    repelem(1, length([H.IndividualParams.OmegaIndex]),1); repelem(.1,length(H.SigmaParams)...
    -length([H.CellParams.OmegaIndex])-length([H.IndividualParams.OmegaIndex]),1)]);
inputs.addParameter('ref', 'params')
inputs.addParameter('LB', [])
inputs.addParameter('UB', [])

inputs.parse(varargin{:});
inputs=inputs.Results;

% Cell-specific params
cell_indx=arrayfun(@(x)(x.name),H.CellParams,'UniformOutput',false);        % Extract field names of individually-varying parameters
try
param_cell=cellfun(@(x)simFun.Parameters.Value(ismember(simFun.Parameters.Name,x)),...
    cell_indx,'UniformOutput',false);                                      % Compare names of parameters and extract the corresponding value
param_cell=reshape(repelem(cell2mat(param_cell),1,length(H.CellParams(1).Index))',[],1);
n_cell = length(param_cell);
catch
    param_cell =[];
    n_cell = 0;
end

% Individual params
indiv_indx=arrayfun(@(x)(x.name),H.IndividualParams,'UniformOutput',false); % Extract field names of individually-varying parameters
try
param_indiv=cellfun(@(x)simFun.Parameters.Value(ismember(simFun.Parameters.Name,x)),...
    indiv_indx,'UniformOutput',false);                                      % Compare names of parameters and extract the corresponding value
param_indiv=reshape(repelem(cell2mat(param_indiv),1,n_sim)',[],1);
n_indiv = length(param_indiv);
catch
    n_indiv = 0;
    param_indiv =[];
end
% Population params
if strcmp(inputs.ref,'params')
    ref=[simFun.Parameters.Value(H.PopulationParams)];       % 
else
    ref = ones(length(H.PopulationParams),1);
end
param = [simFun.Parameters.Value(H.PopulationParams)];       % 

% concatenating all params
if ~isempty(inputs.LB)
        % Min                                   Max                                 Start                       % Mu           
p=[ inputs.LB                             inputs.UB                          param                       param 
    ones(size(param_cell))*1e-3           ones(size(param_cell))*1e2         ones(size(param_cell))      ones(size(param_cell))
    ones(size(param_indiv))*1e-3          ones(size(param_indiv))*1e2          ones(size(param_indiv))     ones(size(param_indiv))     
    sigmaParams*0.1                       sigmaParams*1000                   inputs.startSigma           sigmaParams];

else
    % Min                                   Max                                 Start                       % Mu           
p=[ 10.^(floor(log10(ref*1e-6)))          10.^(ceil(log10(ref*1e6)))         param                       param 
    ones(size(param_cell))*1e-3           ones(size(param_cell))*1e2         ones(size(param_cell))      ones(size(param_cell))
    ones(size(param_indiv))*1e-3          ones(size(param_indiv))*1e2        ones(size(param_indiv))     ones(size(param_indiv))     
    sigmaParams*0.1                       sigmaParams*1000                   inputs.startSigma           sigmaParams];
end
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

% Parameter structure
par_struct=struct('name', paramNames,'minValue',num2cell(p(:,1)), 'maxValue',...
    num2cell(p(:,2)),'startValue', num2cell(p(:,3)),'mu_prior',  num2cell(p(:,4)),'sigma_prior',  num2cell(sigma));
return