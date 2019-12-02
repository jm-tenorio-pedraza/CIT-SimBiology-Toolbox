function PI=getPIData(data_ext,stateVar,groups_subset,observables,varargin)
% Preprocess data for group-wise parameter estimation
p=inputParser;
p.addParameter('zeroAction', 'input')
p.addParameter('mergePhenotypes', false)
p.addParameter('output', 'mean')

p.parse(varargin{:})
p=p.Results;
load(data_ext,'PI')

PI.data = PI.data';
[stateindx,varindx]= (ismember(stateVar, PI.stateVar));
varindx = varindx(varindx~=0);
% varindx = varindx(varindx>0);
data_subset = arrayfun(@(x) x.dataValue(:,varindx), PI.data, 'UniformOutput', false)';
groups = {PI.data(:).Group}';
dataTime = {PI.data(:).dataTime}';
names = {PI.data(:).Name}';
%% 
unique_groups = unique(groups);
data = [];
if strcmp(p.output, 'mean')
for i=1:length(unique_groups)
    group_i = ismember(groups, unique_groups(i));
    
    data_i = data_subset(group_i);
    time_i = {PI.data(group_i).dataTime};
    time = cellfun(@(x) length(x),time_i, 'UniformOutput',true);
    time = time_i{max(time)==time};
    mat_i_responders = NaN(length(time),length(varindx));
    mat_i_progressors = NaN(length(time),length(varindx));
    group_i = {};

    for j=1:length(varindx)
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
                            group_i(1,k) = {'Progressor'};
                        case 0
                            group_i(1,k) = {'Responder'};
                    end
                end
            end
                mat_i_responders(:,j) = mean(mat_j(:, ismember(group_i, 'Responder')),2,'omitnan');
                mat_i_progressors(:,j) = mean(mat_j(:, ismember(group_i, 'Progressor')),2,'omitnan');
           
    end
    % Replacing zero values in TV for minimum tumor volume measured
    if strcmp(p.zeroAction,'input')
    zeroIndx_responders = mat_i_responders(:,1)==0;
    zeroIndx_progressors = mat_i_progressors(:,1)==0;

    
    min_responders = min(mat_i_responders(~zeroIndx_responders,1));
    min_progressors = min(mat_i_progressors(~zeroIndx_progressors,1));
    mat_i_responders(zeroIndx_responders,1) = min_responders;
    mat_i_progressors(zeroIndx_progressors,1) = min_progressors;
    else
    end
    
    data(i).Name = strjoin({unique_groups{i} 'Progressors'},'_');    
    data(i).dataTime = time;
    data(i).dataValue = mat_i_progressors;
    data(i).Group = unique_groups(i);
    data(i).censoring = 'right';
    data(i+length(unique_groups)).Name = strjoin({unique_groups{i} 'Responders'},'_');    
    data(i+length(unique_groups)).dataTime = time;
    data(i+length(unique_groups)).dataValue = mat_i_responders;
    data(i+length(unique_groups)).Group = unique_groups(i);
    data(i+length(unique_groups)).censoring = 'left';
end
nanIndx = arrayfun(@(x) all(all(isnan(x.dataValue))), data);
PI.data = data(:,~nanIndx)';
PI.data = PI.data(ismember([PI.data(:).Group], groups_subset));

else
    data = data_subset';
    nanIndx = cellfun(@(x) (all(isnan(x))), data);
    PI.data = [];
    [PI.data(1:length(data(~nanIndx))).Name] = names{~nanIndx,:};
    [PI.data(1:length(data(~nanIndx))).dataTime] = dataTime{~nanIndx,:};
    [PI.data(1:length(data(~nanIndx))).dataValue] = data{~nanIndx,:};
    [PI.data(1:end).Group] = groups{~nanIndx,:};
    PI.data = PI.data(ismember({PI.data(:).Group}, groups_subset));

end



%% Change group identity

try
PI.data(ismember([PI.data(:).Group], 'MOC1_Control_Mean')).Group = {'MOC1_Control'};
catch
end

try
PI.data(ismember([PI.data(:).Group], 'MOC2_Control_Mean')).Group = {'MOC2_Control'};
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

%% Merge phenotypes
if p.mergePhenotypes
    name = ([PI.data(:).Group]);                                            % Extract names from the groups
    name_indx = cellfun(@(x)regexp(x,'Control'),name,'UniformOutput',false);% Identify those names with Control in it
    name_indx = cellfun(@(x)~isempty(x),name_indx);                         % Obtain indexes for those names corresponding to control groups
    group_i = unique(name(name_indx));                                        % Get unique group names that correspond to the datasets to merge
    data = PI.data;    
    PI.data = data(~name_indx);
    for i = 1:length(group_i)
        indx_i = (ismember(name, group_i(i)));                               % Get indexes of datasets to merge
         name_i = {data(indx_i).Name};

        data_i = data(indx_i);
        time_i = arrayfun(@(x) x.dataTime, data_i, 'UniformOutput', false);
        time_i = unique(cat(1,time_i{:,:}));
        dataset = nan(length(time_i), size(data_i(1).dataValue,2));
        for j=1:length(data_i)
            nan_indx = all(isnan(data_i(j).dataValue));
            time_indx = ismember(time_i, data_i(j).dataTime);
            dataset(time_indx,~nan_indx) = data_i(j).dataValue(:,~nan_indx);
        end
        PI.data(end+1).Name = name_i{1};
        PI.data(end).dataTime = time_i;
        PI.data(end).dataValue = dataset;
        PI.data(end).Group = group_i(i);
        PI.data(end).censoring = {data(find(indx_i,1)).censoring};
    end
end

%% Plot data
ncol = ceil(sqrt(length(varindx)));
nrow = ceil(length(varindx)/ncol);
figure;
colors = table2cell(table(linspecer(length(PI.data))));
[PI.data(1:end).colors] = colors{:,:};
for i = 1:length(varindx)
    subplot(nrow,ncol,i)
    hold on
    arrayfun(@(x)plot(x.dataTime, x.dataValue(:,i),'Color',x.colors,'Marker','*'),PI.data)
    legend({PI.data(:).Name},'Interpreter', 'none')
    indx = find(stateindx,i);
    title(stateVar(indx(i)))
    set(gca,'YScale', 'log')
end
PI.tspan = unique(cat(1,PI.data(:).dataTime));
return
