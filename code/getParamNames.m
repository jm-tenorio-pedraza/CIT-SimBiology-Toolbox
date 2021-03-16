function paramNames = getParamNames(PI,sim,observables,varargin)

par = inputParser;
par.addParameter('all', true)
par.parse(varargin{:})
par = par.Results;

H= PI.H;
n_indiv = length([H.IndividualParams.EtaIndex]);
n_cell = length([H.CellParams.EtaIndex]);
n_resp = length([H.RespParams.EtaIndex]);
n_sigma = length(H.SigmaParams);
n_pop = length(H.PopulationParams);
paramNames = repelem({'nan'}, n_pop+n_indiv+n_cell+n_resp+n_sigma,1);

paramNames([H.PopulationParams]) = [sim.Parameters.Name(H.PopulationParams)];
try
paramNames_cell  = repelem(cellfun(@(x)strjoin({'w' x}, '_'),...
    {H.CellParams(:).name},'UniformOutput', false),length(H.CellParams(1).Index));
psiNames_cell = cellfun(@(x)strjoin({'\psi' x}, '_'),...
    {H.CellParams(:).name},'UniformOutput', false);
catch
    paramNames_cell={};
    psiNames_cell = {};
end
try
paramNames_ind = repelem(cellfun(@(x)strjoin({'z' x}, '_'),...
    {H.IndividualParams(:).name},'UniformOutput', false),length(H.IndividualParams(1).Index));
omegaNames_ind = cellfun(@(x)strjoin({'\omega' x}, '_'),...
    {H.IndividualParams(:).name},'UniformOutput', false);
catch
    paramNames_ind = {};
    omegaNames_ind = {};
end

try
paramNames_resp = repelem(cellfun(@(x)strjoin({'k' x}, '_'),...
    {H.RespParams(:).name},'UniformOutput', false),length(H.RespParams(1).Index));
nuNames_resp = cellfun(@(x)strjoin({'\nu' x}, '_'),...
    {H.RespParams(:).name},'UniformOutput', false);
catch
    paramNames_resp= {};
    nuNames_resp= {};
end
sigmaNames = cellfun(@(x)strjoin({'\sigma' x}, '_'),...
    observables,'UniformOutput', false);
 paramNames([H.IndividualParams(1:end).Index])= paramNames_ind;
 paramNames([H.CellParams(1:end).Index]) = paramNames_cell;
 paramNames([H.RespParams(1:end).Index]) = paramNames_resp;

 paramNames([H.SigmaParams]) = [psiNames_cell'; omegaNames_ind'; nuNames_resp'; sigmaNames'];
 
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