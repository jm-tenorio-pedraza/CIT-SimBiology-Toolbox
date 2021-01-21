function[ PI,doses] = getDoses(PI,varargin)
p = inputParser;
p.addParameter('treatment', 'antiPDL1');
p.addParameter('MW', 1.5e5);
p.addParameter('t0', 0);
p.addParameter('doses', [])
p.addParameter('species', 'mouse')
p.addParameter('ApplicationRateConstant', 60)
p.addParameter('Schedule', [])
p.addParameter('BW', .022);
p.parse(varargin{:})
p=p.Results;

BW = p.BW;
MW = p.MW;

if isempty(p.doses)
    dose = {PI.data(:).Group};
    dose = cellfun(@(x)str2double(regexprep(x,'_mgkg\w*', ''))*BW*1e-3/MW*1e6,dose,'UniformOutput',false);
else

    dose = mat2cell(p.doses*BW*1e-3/p.MW*1e6,repelem(1,size(p.doses,1)),size(p.doses,2));
end
doses = repelem({'nan'}, size(dose,1), size(dose{1},2));
if isempty(p.Schedule)
    doseSchedule = repelem({0},size(dose,1),size(dose{1},2));
else
    doseSchedule = p.Schedule;
end
[PI.data(1:end).doseSchedule] = doseSchedule{:,:};
[PI.data(1:end).dose] = dose{:,:};

for i=1:size(dose{1},2) % For each druggable target
    doses(:,i) = arrayfun(@(x) table(x.doseSchedule(:,i), repelem([x.dose(:,i)], size(x.doseSchedule,1),1),...
        repelem([x.dose(:,i)]/60, size(x.doseSchedule,1),1), 'VariableNames', {'Time' 'Amount' 'Rate'}),...
        PI.data,'UniformOutput',false);
end

[PI.data(1:end).dose] = dose{:,:};
for i = 1:size(doses,1)
    for j=1:size(doses,2)
    doses{i,j}.Properties.VariableUnits= {'hour' 'micromole' 'micromole/second'};
    end
end
    
return
