function PI = generatePI(tablesCell, rowNames,stateVar,varargin)
par = inputParser;
par.addParameter('groupNames', rowNames);
par.addParameter('varIndx', [])
par.parse(varargin{:});
par=par.Results;

PI.data =[];
if isempty(par.varIndx)
    
    varIndx = cellfun(@(x)true(ones(1,length(par.stateVar))),tablesCell, 'UniformOutput',true);
else
    varIndx = par.varIndx;
end
n_data = 0;
nCell = length(tablesCell);
for i=1:nCell
    table_i = tablesCell{i};
    values_conc = [table_i{:,2} table_i{:,4}]; % The concentration values are in columns 2 and 4
    sd_conc = [table_i{:,3} table_i{:,5}]; % The SD values are in columns 3 and 5
    
    dataValues = NaN(size(values_conc,1), size(varIndx,2));
    sdValues = NaN(size(sd_conc,1), size(varIndx,2));

    dataValues(:,varIndx(i,:)) = values_conc;
    sdValues(:,varIndx(i,:)) = sd_conc;
    
    PI.data(i).Name = rowNames(i);
    PI.data(i).dataValue = dataValues;
    PI.data(i).dataTime = [table_i{:,1}];
    PI.data(i).SD = sdValues;
    PI.data(i).Group = par.groupNames(i);
    
    n_data = n_data + sum(sum(~isnan(dataValues)));
end

tspan=(unique(cell2mat({PI.data(:).dataTime})));
PI.stateVar = stateVar;
PI.tspan = tspan;
PI.n_data = n_data;
return