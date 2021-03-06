function PI=getPIData(data_ext,stateVar,groups_subset,observables,varargin)
% Preprocess data for group-wise parameter estimation
p=inputParser;
p.addParameter('zeroAction', 'input')
p.addParameter('mergePhenotypes', false)
p.addParameter('output', 'mean')
p.addParameter('maxIIV', false)
p.parse(varargin{:})
p=p.Results;

load(data_ext,'PI')                                                         % load data structure previously pre-processed
if size(PI.data,1)>size(PI.data,2)                                          % transpose data array for following processing
else
    PI.data = PI.data';
end
[~,varindx]= (ismember(stateVar,PI.stateVar));                              % identify and select variables of interest
data_subset = arrayfun(@(x) x.dataValue(:,varindx), PI.data,...
    'UniformOutput', false)';
groups = {PI.data(:).Group}';                                               % extract groups cell array
dataTime = {PI.data(:).dataTime}';                                          % extract vectors of measured time points cell array
names = {PI.data(:).Name}';                                                 % extract names cell array
%%
unique_groups = unique(groups);                                             % identify unique treatment groups
data = [];
if strcmp(p.output, 'mean')                                                 % processing for reduced mean model
    for i=1:length(unique_groups)                                               % For each unique treatment do
        group_i = ismember(groups, unique_groups(i));                           % Identify and select observations belonging to ith group
        data_i = data_subset(group_i);
        time_i = {PI.data(group_i).dataTime};
        
        time = cellfun(@(x) length(x),time_i, 'UniformOutput',true);            % Count how many time-resolved data points there are for each individual
        time = time_i{max(time)==time};                                         % Select the largest time end point
        mat_i_responders = NaN(length(time),length(observables));               % create NaN matrix with K rows and M columns for non-responders and responders
        mat_i_progressors = NaN(length(time),length(observables));
        group_i = {};
        
        for j=1:length(observables)                                             % For each jth variable do:
            datavar_j = cellfun(@(x) x(:,j),data_i, 'UniformOutput',false);     % Select jth the variable
            mat_j = NaN(length(time),length(datavar_j));                        % NaN mat with K rows and N columns (one for each individual)
            
            for k=1:length(data_i)                                              % For each individual do:
                t_k = length(datavar_j{:,k});                                   % determine the number of observed time points for the kth individual
                mat_j(1:t_k,k) = datavar_j{:,k};                                % assing the observed values in the first t_k rows of the kth column
                if j==1                                                         % For the 1st variable (i.e., tumor volume) check:
                    finalTV = datavar_j{:,k}(end)>0 |...                        % if the final observed time point is = 0 or all of the observations are NaN
                        all(isnan(datavar_j{:,k}));
                    switch finalTV
                        case 1
                            group_i(1,k) = {'Progressor'};                      % if true then assign the 'Progressor' label
                        case 0
                            group_i(1,k) = {'Responder'};                       % otherwise, assign the responder label
                    end
                end
            end
               mat_i_responders(:,j) = mean(mat_j(:, ismember(group_i,...          % assing time course to either the responders or non-responders matrix
                    'Responder')),2,'omitnan');
                mat_i_progressors(:,j) = mean(mat_j(:, ismember(group_i,...
                    'Progressor')),2,'omitnan');
            
            
        end
        
        if strcmp(p.zeroAction,'input')                                         % Check if any data correction for 0 values is to be taken
            zeroIndx_responders = mat_i_responders(:,1)==0;                     % identify rows with 0s in them
            zeroIndx_progressors = mat_i_progressors(:,1)==0;
            min_responders = min(mat_i_responders(~zeroIndx_responders,1));     % Calculate the min for each group
            min_progressors = min(mat_i_progressors(~zeroIndx_progressors,1));
            mat_i_responders(zeroIndx_responders,1) = min_responders;           % assing respective mins to each group where the zeros are
            mat_i_progressors(zeroIndx_progressors,1) = min_progressors;
        else
        end
        data(i).Name = strjoin({unique_groups{i} 'Progressors'},'_');           % First rows of data structure contain the progressors
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
    nanIndx = arrayfun(@(x) all(all(isnan(x.dataValue))), data);                % Identify and eliminate individuals where all obs are nan
    PI.data = data(:,~nanIndx)';
    PI.data = PI.data(ismember([PI.data(:).Group], groups_subset));             % Subset data to use only the input groups
    
else
    data = data_subset';                                                    % If no mean-model is to be used then do:
    nanIndx = cellfun(@(x) (all(isnan(x))), data);                          % Identify and eliminate individuals where all obs are nan
    PI.data = [];
    [PI.data(1:length(data(~nanIndx))).Name] = names{~nanIndx,:};           % subset all individual-level data
    [PI.data(1:length(data(~nanIndx))).dataTime] = dataTime{~nanIndx,:};
    [PI.data(1:length(data(~nanIndx))).dataValue] = data{~nanIndx,:};
    [PI.data(1:end).Group] = groups{~nanIndx,:};
    PI.data = PI.data(ismember({PI.data(:).Group}, groups_subset));
    if p.maxIIV
        
        
    end
    
end



%% Change group identity

try
    PI.data(ismember([PI.data(:).Group], 'MOC1_Control_Mean')).Group ...
        = {'MOC1_Control'};
catch
end

try
    PI.data(ismember([PI.data(:).Group], 'MOC2_Control_Mean')).Group ...
        = {'MOC2_Control'};
catch
end


%% Control for data with non-admissible values (x=0)

zero_indx = arrayfun(@(x) sum(x.dataValue(:,1)==0,2),PI.data,...            
    'UniformOutput', false);                                                % Identify indexes of tumor volume equal to 0
[PI.data(1:end).zero_indx] = zero_indx{:,:};                                % Add index array to data array
dataValue = arrayfun(@(x) x.dataValue(~x.zero_indx,:), PI.data,...          % Subset out the zero values and their respective time points
    'UniformOutput', false);
dataTime = arrayfun(@(x) x.dataTime(~x.zero_indx,:), PI.data,...
    'UniformOutput', false);        
[PI.data(1:end).dataValue] = dataValue{:,:};                                
[PI.data(1:end).dataTime] = dataTime{:,:};
PI.n_data=sum(cellfun(@(x)sum(sum(~isnan(x))),{PI.data.dataValue},...
    'UniformOutput',true));                                                 % Count number of non-nan data

%% Merge phenotypes
if p.mergePhenotypes
    name = ([PI.data(:).Group]);                                            % Extract names from the groups
    name_indx = cellfun(@(x)regexp(x,'Control'),name,'UniformOutput',false);% Identify  names with 'Control' in it
    name_indx = cellfun(@(x)~isempty(x),name_indx);                         % Obtain indexes for those names corresponding to control groups
    group_i = unique(name(name_indx));                                      % Get unique group names that correspond to the datasets to merge
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
ncol = ceil(sqrt(length(stateVar)));
nrow = ceil(length(stateVar)/ncol);
figure;
colors = table2cell(table(linspecer(length(PI.data))));
[PI.data(1:end).colors] = colors{:,:};
for i = 1:length(stateVar)
    subplot(nrow,ncol,i)
    hold on
    arrayfun(@(x)plot(x.dataTime, x.dataValue(:,i),'Color',...
        x.colors,'Marker','*'),PI.data)
    legend({PI.data(:).Name},'Interpreter', 'none')
    title(stateVar(i))
    set(gca,'YScale', 'log')
end
PI.tspan = unique(cat(1,PI.data(:).dataTime));
return
