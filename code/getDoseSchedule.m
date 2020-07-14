function u = getDoseSchedule(PI,doses,varargin)
par = inputParser;
par.addParameter('doseUnits', 'g')
par.addParameter('addVariants',false)
par.addParameter('parallel',false)
par.addParameter('startTime', zeros(1,length(doses)));
par.addParameter('freq', zeros(1,length(doses)));
par.addParameter('numDoses', zeros(1,length(doses)));
par.addParameter('dose', zeros(1,length(doses)));
par.parse(varargin{:})
par=par.Results;

[startTime, freq, numDoses,dose] = deal(par.startTime, par.freq, par.numDoses, par.dose);

antiPDL1=table([startTime(1):freq(1):(freq(1)*(numDoses(1)-1)+startTime(1))]',...
    repelem(dose(1),numDoses(1))', repelem(dose(1),numDoses(1))'/60,...
    'VariableNames',{'Time' 'Amount' 'Rate'});
antiCTLA4=table([startTime(2):freq(2):(freq(2)*(numDoses(2)-1)+startTime(2))]',...
    repelem(dose(2),numDoses(2))', repelem(dose(2),numDoses(2))'/60,...
    'VariableNames',{'Time' 'Amount' 'Rate'});
control=table([startTime(1):freq(1):(freq(1)*(numDoses(1)-1)+startTime(1))]',...
    zeros(numDoses(1),1), zeros(numDoses(1),1),...
    'VariableNames',{'Time' 'Amount' 'Rate'});
if strcmp(par.doseUnits, 'g')
% Create cell of doses
antiPDL1.Properties.VariableUnits={'day' 'milligram' 'milligram/second'};

antiCTLA4.Properties.VariableUnits={'day' 'milligram' 'milligram/second'};

control.Properties.VariableUnits={'day' 'milligram' 'milligram/second'};
else
antiPDL1.Properties.VariableUnits={'day' 'micromole' 'micromole/second'};

antiCTLA4.Properties.VariableUnits={'day' 'micromole' 'micromole/second'};

control.Properties.VariableUnits={'day' 'micromole' 'micromole/second'};
end    

u=repelem({NaN},size(PI.data,1),length(doses));
dose = cellfun(@(x) regexp(x,'_','split'),doses,'UniformOutput',false); 
dose = [dose{:,:}];
dose = dose(2:2:end);
group = [PI.data(:).Group]';
if ischar(group)
    group = {PI.data(:).Group}';
end
group = cellfun(@(x) regexp(x,'_', 'split'),group, 'UniformOutput',false);

for i=1:length(group)
    group_i = group{i,1}(1,2:end);
    [bool,indx] = ismember(dose, group_i);
    cellindx = 1:length(dose);
    for j = 1:length(dose)
        if bool(j)
            doseindx = cellindx(j);
                if doseindx==1
                    u(i,j)={antiPDL1};
                elseif doseindx==2
                    u(i,j)={antiCTLA4};
                end
        else
           u(i,j)={control};

        end
    end
end

u=u(:,1:length(dose));

ncol = ceil(sqrt(size(u,1)));
nrow = ceil(size(u,1)/ncol);

figure
for i=1:size(u,1)
    subplot(nrow,ncol,i)
    hold on
    for j=1:size(u,2)
        plot(u{i,j}.Time, u{i,j}.Amount, '*')
       
    end
    plot(u{i,j}.Time, zeros(length(u{i,j}.Time),1), '-k')

    legend(dose)
    title(PI.data(i).Group,'Interpreter', 'none')
%     ylim([0 1])
end
return

