function [sim,u]=initializePI(sim_model,parameters,observables,PI)

variants=getvariant(sim_model);

MOC1=variants(strcmp(get(variants,'Name'), 'MOC1'));
% Define dose
control = sbiodose('rd');
control.TargetName = 'Dose_antiPDL1';
control.TimeUnits = 'day';
control.AmountUnits = 'micromole';
control.RateUnits = 'micromole/second';

% Create simFunction object
sim=createSimFunction(sim_model,parameters,observables, control,MOC1,...
    'UseParallel', false);

% Create cell of doses
antiPDL1_table=table([7 12 17]', [0.0013 0.0013 0.0013]', [0 0 0]',...
    'VariableNames',{'Time' 'Amount' 'Rate'});
antiPDL1_table.Properties.VariableUnits={'day' 'micromole' 'micromole/second'};

control_table=table([7 12 17]', [0 0 0]', [0 0 0]',...
    'VariableNames',{'Time' 'Amount' 'Rate'});
control_table.Properties.VariableUnits={'day' 'micromole' 'micromole/second'};
u=repelem({NaN},size(PI.data,1),1);
for i=1:size(PI.data,1)
    if strcmp(PI.data(i).Group, 'MOC1_Control')
        u(i,:)={control_table};
    elseif strcmp(PI.data(i).Group, 'MOC1_antiPDL1')
        u(i,:)={antiPDL1_table};
    end
end

