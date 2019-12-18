function paramNames = getParamNames(PI,sim,observables,varargin)

par = inputParser;
par.addParameter('all', true)
par.parse(varargin{:})
par = par.Results;

H= PI.H;
n_indiv = length(PI.H.IndividualParams);
n_cell = length(PI.H.CellParams);
n_sigma = length(PI.H.SigmaParams);
n_pop = length(PI.H.PopulationParams);
paramNames = repelem({'nan'}, n_pop+n_indiv+n_cell+n_sigma,1);

paramNames([H.PopulationParams]) = [sim.Parameters.Name(H.PopulationParams)];
try
paramNames_cell  = repelem(cellfun(@(x)strjoin({'\gamma' x}, '_'),...
    {H.CellParams(:).name},'UniformOutput', false),length(H.CellParams(1).Index));
psiNames_cell = cellfun(@(x)strjoin({'\psi' x}, '_'),...
    {H.CellParams(:).name},'UniformOutput', false);
catch
    paramNames_cell={};
    psiNames_cell = {};
end
try
paramNames_ind = repelem(cellfun(@(x)strjoin({'\eta' x}, '_'),...
    {H.IndividualParams(:).name},'UniformOutput', false),length(H.IndividualParams(1).Index));
omegaNames_ind = cellfun(@(x)strjoin({'\omega' x}, '_'),...
    {H.IndividualParams(:).name},'UniformOutput', false);
catch
    paramNames_ind = {};
    omegaNames_ind = {};
end
sigmaNames = cellfun(@(x)strjoin({'\sigma' x}, '_'),...
    observables,'UniformOutput', false);
 paramNames([H.IndividualParams(1:end).Index])= paramNames_ind;
 paramNames([H.CellParams(1:end).Index]) = paramNames_cell;
 paramNames([H.SigmaParams]) = [omegaNames_ind'; psiNames_cell'; sigmaNames'];
 
 for i=1:length(paramNames)
     underscore_indx = ismember(paramNames{i}, '_');
     if any(underscore_indx)
         paramNames_i =strrep(paramNames{i}, '_', '_{');
         paramNames_i(end+1:end+sum(underscore_indx)) = repelem('}', 1,sum(underscore_indx));
         paramNames(i) = {paramNames_i};
     else
         continue
     end
 end