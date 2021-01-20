function conc_table = getIDperMl(conc_table, dose)

concentration = conc_table{:,2};
conc_norm = concentration/dose*100;
sd = conc_table{:,3};
sd_norm = sd/dose*100;
conc_table{:,4:5} = [conc_norm sd_norm];
conc_table.Properties.VariableNames(4:5) = {'Concentration [%ID/ml]' 'SD[%ID/ml]'};
return

