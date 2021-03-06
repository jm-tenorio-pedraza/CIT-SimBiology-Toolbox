function par_struct = getParamStruct(simFun,H,n_sim,sigmaParams,sigmaNames,varargin)
if nargin<5
    error('getParamStruct:toofewinputs','getParamStruct requires atleast 5 inputs.')
end
inputs=inputParser;
inputs.addParameter('Sigma',3,@isnumeric);
inputs.parse(varargin{:});
inputs=inputs.Results;


% Population params
param=simFun.Parameters.Value(H.PopulationParams);

% Individual params
indiv_indx=arrayfun(@(x)(x.name),H.IndividualParams,'UniformOutput',false);
param_indiv=cellfun(@(x)simFun.Parameters.Value(ismember(simFun.Parameters.Name,x)),indiv_indx,'UniformOutput',false);

param_indiv=reshape(repmat(cell2mat(param_indiv),1,n_sim),[],1);
% concatenating all params

    % Min            Max                Start           % Mu             % Sigma
p=[ param*0.01       param*100           param           param          
    param_indiv*.01  param_indiv*100     param_indiv     param_indiv     
    sigmaParams*.01  sigmaParams*100     sigmaParams     sigmaParams];
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
    repelem({(H.IndividualParams(1:end).name)},n_sim,1);...
    sigmaNames];

% Parameter structure
par_struct=struct('name', paramNames,'minValue',num2cell(p(:,1)), 'maxValue',...
    num2cell(p(:,2)),'startValue', num2cell(p(:,3)),'mu_prior',  num2cell(p(:,4)),'sigma_prior',  num2cell(sigma));
return