% Preprocess data for group-wise parameter estimation

load('PI_Clavijo.mat')
PI.data = PI.data';
[~,varindx]= (ismember(stateVar,PI.stateVar));
% varindx = varindx(varindx>0);
data_subset = arrayfun(@(x) x.dataValue(:,varindx), PI.data, 'UniformOutput', false)';

groups = unique({PI.data(:).Group});
data = [];
for i=1:length(groups)
    group_i = ismember({PI.data(:).Group}, groups(i));
    
    data_i = [data_subset(group_i)];
    time_i = {PI.data(group_i).dataTime};
    time = cellfun(@(x) length(x),time_i, 'UniformOutput',true);
            time = time_i{max(time)==time};
    mat_i = NaN(length(time),length(observables));

    for j=1:length(observables)
            datavar_j = cellfun(@(x) x(:,j),data_i, 'UniformOutput',false);
            
            % Matrix with each individual in the columns and the time
            % poitns in the rows
            mat_j = NaN(length(time),length(datavar_j));
            for k=1:length(data_i)
                mat_j(1:length(datavar_j{:,k}),k) = datavar_j{:,k};
            end
            mat_i(:,j) = mean(mat_j,2,'omitnan');
            
            
    end
        data(i).Name = groups(i);
        

    data(i).dataTime = time;
    data(i).dataValue = mat_i;
    data(i).Group = groups(i);
end
PI.data = data';

%% Subsetting to use only Control and antiPDL1

PI.data = PI.data(ismember([PI.data(:).Group], groups_subset));
PI.data(ismember([PI.data(:).Group], 'MOC1_Control_Mean')).Group = 'MOC1_Control';
