function PI=getPIData2(data_ext,stateVar,groups_subset,observables,varargin)
% Preprocess data for group-wise parameter estimation
p=inputParser;
p.addParameter('zeroAction', 'input')
p.addParameter('mergePhenotypes', false)
p.addParameter('output', 'mean')
p.addParameter('maxIIV', false)
p.addParameter('f_training', 0.6)

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
%%
unique_groups = unique(groups);                                             % identify unique treatment groups
if strcmp(p.output, 'mean')                                                 % processing for reduced mean model
    data_cell = repelem({'nan'}, length(unique_groups),1);
else
    data_cell = repelem({'nan'}, length(groups),1);
end
data=[];
% [data(1:end).dataValue] = data_cell{:,:};

struct_indx = 1;
for i=1:length(unique_groups)                                               % For each unique treatment do
    group_i = ismember(groups, unique_groups(i));                           % Identify and select observations belonging to ith group
    data_i = data_subset(group_i);
    time_i = dataTime(group_i);
    
    time = cellfun(@(x) length(x),time_i, 'UniformOutput',true);            % Count how many time-resolved data points there are for each individual
    time = time_i{max(time)==time};                                         % Select the largest time end point
    mat_i_responders = NaN(length(time),length(observables));               % create NaN matrix with K rows and M columns for non-responders and responders
    mat_i_progressors = NaN(length(time),length(observables));
    n_indiv = length(data_i);
    
    for j=1:length(observables)                                             % For each jth variable do:
        datavar_j = cellfun(@(x) x(:,j),data_i, 'UniformOutput',false);     % Select jth the variable
        mat_j = NaN(length(time),n_indiv);                        % NaN mat with K rows and N columns (one for each individual)
        
        for k=1:n_indiv                                          % For each individual do:
            t_k = length(datavar_j{:,k});                                   % determine the number of observed time points for the kth individual
            mat_j(1:t_k,k) = datavar_j{:,k};                                % assing the observed values in the first t_k rows of the kth column
            if j==1                                                         % For the 1st variable (i.e., tumor volume) check:
                if k==1
                    phenotype_i = repelem({'nan'}, n_indiv,1);
                end
                finalTV = datavar_j{:,k}(end)>0.01 |...                        % if the final observed time point is = 0 or all of the observations are NaN
                    all(isnan(datavar_j{:,k}));
                switch finalTV
                    case 1
                        phenotype_i(k) = {'Progressor'};                      % if true then assign the 'Progressor' label
                    case 0
                        phenotype_i(k) = {'Responder'};                       % otherwise, assign the responder label
                end
            end
        end
        if strcmp(p.output, 'mean')                                                 % processing for reduced mean model
            
            mat_i_responders(:,j) = mean(mat_j(:, ismember(phenotype_i,...          % assing time course to either the responders or non-responders matrix
                {'Responder'})),2,'omitnan');
            mat_i_progressors(:,j) = mean(mat_j(:, ismember(phenotype_i,...
                {'Progressor'})),2,'omitnan');
        end
    end
    n_progressors = sum(ismember(phenotype_i,{'Progressor'}));
    n_responders = sum(ismember(phenotype_i,{'Responder'}));
    if strcmp(p.output, 'mean')                                                 % Assignment to structure
        
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
    else
        progressor_indx = ismember(phenotype_i,'Progressor');
        responder_indx = ismember(phenotype_i,'Responder');

        if p.maxIIV
            
            data_ij = data_i(progressor_indx);
            time_ij = time_i(progressor_indx);

            train_ij = getMaxVarData(data_ij,time_ij,floor(n_progressors*p.f_training));
            name_progressors = repelem({strjoin({unique_groups{i}...
                'Progressors'},'_')},n_progressors,1)  ;
            groups_ij = repelem(unique_groups(i),n_progressors,1);
            
            censoring_ij = repelem({'right'}, n_progressors,1);
            [data(struct_indx:struct_indx+n_progressors-1).Name] = name_progressors{:,:};         % First rows of data structure contain the progressors
            [data((struct_indx:struct_indx+n_progressors-1)).dataTime] = time_ij{:,:};
            [data((struct_indx:struct_indx+n_progressors-1)).dataValue] = data_ij{:,:};
            [data((struct_indx:struct_indx+n_progressors-1)).Group] = groups_ij{:,:};
            [data((struct_indx:struct_indx+n_progressors-1)).censoring] =censoring_ij{:,:};
            [data((struct_indx:struct_indx+n_progressors-1)).use] =train_ij{:,:};
            
            struct_indx = struct_indx + n_progressors;

            data_ij = data_i(responder_indx);
            time_ij = time_i(responder_indx);
            train_ij = getMaxVarData(data_ij,time_ij,floor(n_responders*p.f_training));
            name_responders = repelem({strjoin({unique_groups{i}...
                'Responders'},'_')},n_responders,1)  ;
            groups_ij = repelem(unique_groups(i),n_responders,1);
                        
            censoring_ij = repelem({'left'}, n_responders,1);
            [data(struct_indx:struct_indx+n_responders-1).Name] = name_responders{:,:};         % First rows of data structure contain the progressors
            [data((struct_indx:struct_indx+n_responders-1)).dataTime] = time_ij{:,:};
            [data((struct_indx:struct_indx+n_responders-1)).dataValue] = data_ij{:,:};
            [data((struct_indx:struct_indx+n_responders-1)).Group] = groups_ij{:,:};
            [data((struct_indx:struct_indx+n_responders-1)).censoring] =censoring_ij{:,:};
            [data((struct_indx:struct_indx+n_responders-1)).use] =train_ij{:,:};
            struct_indx = struct_indx+n_responders;

        else
            name_progressors = repelem({strjoin({unique_groups{i} 'Progressors'},'_')},n_progressors,1)  ;
            time_ij = time_i(progressor_indx);
            data_ij = data_i(progressor_indx);
            groups_ij = repelem(unique_groups(i),n_progressors,1);
            censoring_ij = repelem({'right'}, n_progressors,1);
            [data(struct_indx:struct_indx+n_progressors-1).Name] = name_progressors{:,:};         % First rows of data structure contain the progressors
            [data((struct_indx:struct_indx+n_progressors-1)).dataTime] = time_ij{:,:};
            [data((struct_indx:struct_indx+n_progressors-1)).dataValue] = data_ij{:,:};
            [data((struct_indx:struct_indx+n_progressors-1)).Group] = groups_ij{:,:};
            [data((struct_indx:struct_indx+n_progressors-1)).censoring] =censoring_ij{:,:};
            
            struct_indx = struct_indx + n_progressors;
            
            name_responders = repelem({strjoin({unique_groups{i} 'Responders'},'_')},n_responders,1)  ;
            time_ij = time_i(responder_indx);
            data_ij = data_i(responder_indx);
            groups_ij = repelem(unique_groups(i),n_responders,1);
            censoring_ij = repelem({'left'}, n_responders,1);
            
            [data(struct_indx:struct_indx+n_responders-1).Name] = name_responders{:,:};         % First rows of data structure contain the progressors
            [data((struct_indx:struct_indx+n_responders-1)).dataTime] = time_ij{:,:};
            [data((struct_indx:struct_indx+n_responders-1)).dataValue] = data_ij{:,:};
            [data((struct_indx:struct_indx+n_responders-1)).Group] = groups_ij{:,:};
            [data((struct_indx:struct_indx+n_responders-1)).censoring] = censoring_ij{:,:};
            
            struct_indx = struct_indx+n_responders;
        end
    end
end

PI.data = PI.data;
nanIndx = arrayfun(@(x) all(all(isnan(x.dataValue))), data);                % Identify and eliminate individuals where all obs are nan
PI.data = data(:,~nanIndx)';
try
    PI.data = PI.data(ismember({PI.data(:).Group}, groups_subset));             % Subset data to use only the input groups
catch
    PI.data = PI.data(ismember([PI.data(:).Group], groups_subset));             % Subset data to use only the input groups
end

%% Change group identity

try
    if p.maxIIV
    PI.data(ismember({PI.data(:).Group}, 'MOC1_Control_Mean')).use ...
        = {'train'};
    end
    PI.data(ismember([PI.data(:).Group], 'MOC1_Control_Mean')).Group ...
        = {'MOC1_Control'};
catch
end

try
    if p.maxIIV
    PI.data(ismember({PI.data(:).Group}, 'MOC2_Control_Mean')).use ...
        = {'train'};
    end
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
if strcmp(p.output, 'mean')
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
else
    ncol = ceil(sqrt(length(groups_subset)));
    nrow = ceil(length(groups_subset)/ncol);
    colors = table2cell(table(linspecer(length(PI.data))));
    [PI.data(1:end).colors] = colors{:,:};
    marker = repelem({'nan'},length(PI.data),1);
    marker(ismember([PI.data(1:end).use],{'test'}))= {'*'};
    marker(ismember([PI.data(1:end).use],{'train'}))= {'d'};
    [PI.data(1:end).marker] = marker{:,:};

    for j=1:length(stateVar)
        if j==1
        figure;

        for i = 1:length(groups_subset)
            group_i = ismember([PI.data(:).Group], groups_subset(i));
            data_i = PI.data(group_i);
            
            subplot(nrow,ncol,i)
            hold on
            h=arrayfun(@(x)plot(x.dataTime, x.dataValue(:,j),'Color',...
                x.colors,'Marker',x.marker),data_i);
            try
                legend(h,{data_i(1:end).use});
            catch
                 legend(h,[data_i(1:end).use]);

            end
            title(groups_subset(i),'Interpreter', 'none')
%             set(gca,'YScale', 'log')
        end 
        end
    end
end
tv = arrayfun(@(x) x.dataValue(~isnan(x.dataValue(:,1)),1),PI.data,'UniformOutput',false);
[PI.data(1:end).TV] = tv{:,:};
if p.maxIIV
    PI.test = PI.data(ismember([PI.data(1:end).use],{'test'}));
    PI.data = PI.data(ismember([PI.data(1:end).use],{'train'}));
end
PI.tspan = unique(cat(1,PI.data(:).dataTime));
PI.n_data=sum(cellfun(@(x)sum(sum(~isnan(x))),{PI.data.dataValue},...
    'UniformOutput',true));

if size(PI.data,1)>size(PI.data,2)
else
    PI.data=PI.data';
end
return
