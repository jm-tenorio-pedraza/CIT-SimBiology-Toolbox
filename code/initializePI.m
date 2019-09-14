function [sim,u]=initializePI(sim_model,parameters,observables,PI,doses,variant)

variants=getvariant(sim_model);

MOC1=variants(strcmp(get(variants,'Name'), variant));

% Create simFunction object
sim=createSimFunction(sim_model,parameters,observables, doses,MOC1,...
    'UseParallel', false);

% Create cell of doses
antiPDL1=table([7 12 17]', [0.0013 0.0013 0.0013]', [0 0 0]',...
    'VariableNames',{'Time' 'Amount' 'Rate'});
antiPDL1.Properties.VariableUnits={'day' 'micromole' 'micromole/second'};

antiCTLA4=table([7 12 17]', [0.00065 0.00065 0.00065]', [0 0 0]',...
    'VariableNames',{'Time' 'Amount' 'Rate'});
antiCTLA4.Properties.VariableUnits={'day' 'micromole' 'micromole/second'};

control=table([7 12 17]', [0 0 0]', [0 0 0]',...
    'VariableNames',{'Time' 'Amount' 'Rate'});
control.Properties.VariableUnits={'day' 'micromole' 'micromole/second'};
u=repelem({NaN},size(PI.data,1),length(doses));
dose = cellfun(@(x) regexp(x,'_','split'),doses,'UniformOutput',false); 
dose = [dose{:,:}];
dose = dose(2:2:end);
group = cellfun(@(x) regexp(x,'_', 'split'),[PI.data(:).Group]', 'UniformOutput',false);

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
return

