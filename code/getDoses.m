function[ PI,doses] = getDoses(PI,varargin)
p = inputParser;
p.addParameter('treatment', 'antiPDL1');
p.addParameter('MW', 1.5e5);
p.addParameter('t0', 0);
p.addParameter('doses', [])
p.addParameter('species', 'mouse')
p.parse(varargin{:})
p=p.Results;

if strcmp(p.species, 'mouse')
    BW = 0.022;
elseif strcmp(p.species, 'human')
    BW = 77;
end
if isempty(p.doses)
    dose = {PI.data(:).Group};
    dose = cellfun(@(x)str2double(strrep(x,'_mgkg', ''))*p.BW*1e-3/p.MW*1e6,dose,'UniformOutput',false);
else
    dose = repelem({'nan'}, size(p.doses,1),size(p.doses,2));

    for i=1:size(p.doses,2)
        dose(:,i) = num2cell(p.doses(:,i)*BW*1e-3/p.MW*1e6);
    end
end
doses = repelem({'nan'}, size(dose,1), size(dose,2));
for i=1:size(dose,2)
[PI.data(1:end).dose] = dose{:,i};

doses(:,i) = arrayfun(@(x) table(0, [x.dose],...
     [x.dose]/60, 'VariableNames', {'Time' 'Amount' 'Rate'}),...
    PI.data,'UniformOutput',false);
end

[PI.data(1:end).dose] = dose{:,1};
for i = 1:size(doses,1)
    for j=1:size(doses,2)
    doses{i,j}.Properties.VariableUnits= {'hour' 'microgram' 'microgram/second'};
    end
end
    
return
