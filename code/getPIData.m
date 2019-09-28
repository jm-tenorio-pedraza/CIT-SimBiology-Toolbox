function PI=getPIData(data_ext,stateVar,groups_subset,observables)
% Preprocess data for group-wise parameter estimation

load(data_ext)
PI.data = PI.data';
[~,varindx]= (ismember(stateVar,PI.stateVar));
% varindx = varindx(varindx>0);
data_subset = arrayfun(@(x) x.dataValue(:,varindx), PI.data, 'UniformOutput', false)';

groups = unique({PI.data(:).Group});
data = [];
for i=1:length(groups)
    group_i = ismember({PI.data(:).Group}, groups(i));
    
    data_i = data_subset(group_i);
    time_i = {PI.data(group_i).dataTime};
    time = cellfun(@(x) length(x),time_i, 'UniformOutput',true);
            time = time_i{max(time)==time};
    mat_i_responders = NaN(length(time),length(observables));
    mat_i_progressors = NaN(length(time),length(observables));
    name_i = {};

    for j=1:length(observables)
            datavar_j = cellfun(@(x) x(:,j),data_i, 'UniformOutput',false);
            
            % Matrix with each individual in the columns and the time
            % poitns in the rows
            mat_j = NaN(length(time),length(datavar_j));

            for k=1:length(data_i)
                mat_j(1:length(datavar_j{:,k}),k) = datavar_j{:,k};
                if j==1
                    finalTV = datavar_j{:,k}(end)>0 | all( isnan(datavar_j{:,k}));
                    switch finalTV
                        case 1
                            name_i(1,k) = {'Progressor'};
                        case 0
                            name_i(1,k) = {'Responder'};
                    end
                end
            end
                mat_i_responders(:,j) = mean(mat_j(:, ismember(name_i, 'Responder')),2,'omitnan');
                mat_i_progressors(:,j) = mean(mat_j(:, ismember(name_i, 'Progressor')),2,'omitnan');
           
    end
    % Replacing zero values in TV for minimum tumor volume measured
    zeroIndx_responders = mat_i_responders(:,1)==0;
    zeroIndx_progressors = mat_i_progressors(:,1)==0;

    
    min_responders = min(mat_i_responders(~zeroIndx_responders,1));
    min_progressors = min(mat_i_progressors(~zeroIndx_progressors,1));
    mat_i_responders(zeroIndx_responders,1) = min_responders;
    mat_i_progressors(zeroIndx_progressors,1) = min_progressors;
    
    data(i).Name = strjoin({groups{i} 'Progressors'},'_');    
    data(i).dataTime = time;
    data(i).dataValue = mat_i_progressors;
    data(i).Group = groups(i);
    
    data(i+length(groups)).Name = strjoin({groups{i} 'Responders'},'_');    
    data(i+length(groups)).dataTime = time;
    data(i+length(groups)).dataValue = mat_i_responders;
    data(i+length(groups)).Group = groups(i);
    
    
    
    
end


nanIndx = arrayfun(@(x) all(all(isnan(x.dataValue))), data);

PI.data = data(:,~nanIndx)';

%% Subsetting to use only Control and antiPDL1

PI.data = PI.data(ismember([PI.data(:).Group], groups_subset));
try
PI.data(ismember([PI.data(:).Group], 'MOC1_Control_Mean')).Group = 'MOC1_Control';
catch
end
%% Control for data with non-admissible values (0)

zero_indx = arrayfun(@(x) sum(x.dataValue(:,1)==0,2),PI.data,'UniformOutput', false);   % Indexes of tumor volume equal to 0

[PI.data(1:end).zero_indx] = zero_indx{:,:};
dataValue = arrayfun(@(x) x.dataValue(~x.zero_indx,:), PI.data, 'UniformOutput', false);
dataTime = arrayfun(@(x) x.dataTime(~x.zero_indx,:), PI.data, 'UniformOutput', false);

[PI.data(1:end).dataValue] = dataValue{:,:};
[PI.data(1:end).dataTime] = dataTime{:,:};
PI.n_data=sum(cellfun(@(x)sum(sum(~isnan(x))),{PI.data.dataValue},'UniformOutput',true));

%% Plot data

figure;
hold on
arrayfun(@(x)plot(x.dataTime, x.dataValue(:,1),'*'),PI.data)
legend({PI.data(:).Name},'Interpreter', 'none')
set(gca,'YScale', 'log')
return
