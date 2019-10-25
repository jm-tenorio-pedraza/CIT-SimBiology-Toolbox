function[ PI,doses] = getDoses(PI,varargin)
p = inputParser;
p.addParameter('treatment', 'antiPDL1');
p.addParameter('normFactor', 0.022);
p.addParameter('t0', 0);

p.parse(varargin{:})
p=p.Results;

dose = {PI.data(:).Group};
dose = cellfun(@(x)str2double(strrep(x,'_mgkg', ''))*p.normFactor,dose,'UniformOutput',false);
[PI.data(1:end).dose] = dose{:,:};
doses = arrayfun(@(x) table(0, [x.dose],...
     [x.dose]/180, 'VariableNames', {'Time' 'Amount' 'Rate'}),...
    PI.data,'UniformOutput',false);
for i = 1:length(doses)
    doses{i}.Properties.VariableUnits= {'hour' 'milligram' 'microgram/second'};
end
return
