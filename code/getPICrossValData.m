function MetaPI = getPICrossValData(data_ext,stateVar,groups_subset,observables,varargin)
% Preprocess data for group-wise parameter estimation
p=inputParser;
p.addParameter('zeroAction', 'input')
p.addParameter('mergePhenotypes', false)
p.addParameter('output', 'mean')
p.addParameter('crossvalidation', false)
p.addParameter('n_folds', 5)

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
struct_indx = 1;

MetaPI = [];
unique_groups = unique(groups);                                             % identify unique treatment groups
phen_group = repelem({'nan'},1, length(groups));
for i=1:length(groups)
    finalTV = data_subset{i}(end,1)<=0;
    switch finalTV
        case 1
            phen_group(1,i) = {'Responder'};                      % if true then assign the 'Progressor' label
        case 0
            phen_group(1,i) = {'Progressor'};                       % otherwise, assign the responder label
    end
end

if strcmp(p.output, 'mean')                                                 % processing for reduced mean model
for i=1:length(unique_groups)                                               % For each unique treatment do
    group_i = ismember(groups, unique_groups(i));                           % Identify and select observations belonging to ith group
    data_i = data_subset(group_i);
    time_i = {PI.data(group_i).dataTime};
    time = cellfun(@(x) length(x),time_i, 'UniformOutput',true);            % Count how many time-resolved data points there are for each individual
    time = time_i{max(time)==time};                                         % Select the largest time end point
    phen_group_i = phen_group(group_i);
    for j=1:length(observables)                                             % For each jth variable do:
        datavar_j = cellfun(@(x) x(:,j),data_i, 'UniformOutput',false);     % Select jth the variable
        mat_j = NaN(length(time),length(datavar_j));                        % NaN mat with K rows and N columns (one for each individual)

        for k=1:length(datavar_j)                                           % For each individual do:
            t_k = length(datavar_j{:,k});                                   % determine the number of observed time points for the kth individual
            mat_j(1:t_k,k) = datavar_j{:,k};                                % assing the observed values in the first t_k rows of the kth column
        end
        
        if j==1
            train_i_responders = NaN(length(time),1);     % create NaN matrix with K rows and M columns for non-responders and responders
            train_i_progressors = NaN(length(time),1);
            test_i_responders = NaN(length(time),1);      
            test_i_progressors = NaN(length(time),1);
            n_i_responders = sum(ismember(phen_group_i, 'Responder'));      % Count number of responders in ith group (N_i_r)
            n_i_progressors = sum(ismember(phen_group_i,'Progressor'));     % Count number of progressors in ith group (N_i_p)
            n_cv_r = floor(n_i_responders*1/p.n_folds);                      % Determine size of test set in responder group
            n_cv_p = floor(n_i_progressors*1/p.n_folds);                     % Determine size of test set in progressor group
            cv_r = randperm(n_i_responders);
            cv_p = randperm(n_i_progressors);                               % Sample test subjexts w/o replacement
            fold_indx_r = 1;
            fold_indx_p = 1;
            for l=1:p.n_folds
                test_indx_p = cv_p(fold_indx_p:fold_indx_p+n_cv_p-1);
                train_indx_p = setdiff(1:n_i_progressors,test_indx_p);
                mat_i_progressors = mat_j(:, ismember(phen_group_i,...
                    'Progressor'));
                train_i_progressors(:,j) = mean(mat_i_progressors(:,...      % Calculate mean over the training set
                    train_indx_p),2,'omitnan');
                test_i_progressors(:,j) = mean(mat_i_progressors(:,...      % Calculate mean over the training set
                    test_indx_p),2,'omitnan');
                fold_indx_p = fold_indx_p +n_cv_p;
                try
                    test_indx_r = cv_r(fold_indx_r:fold_indx_r+n_cv_r-1);
                    train_indx_r = setdiff(1:n_i_responders,test_indx_r);
                    mat_i_responders = (mat_j(:, ismember(phen_group_i,...          % assign time course to either the responders or non-responders matrix
                        'Responder')));
                    train_i_responders(:,j) = mean(mat_i_responders(:,...
                        train_indx_r),2,'omitnan');
                    test_i_responders(:,j) = mean(mat_i_responders(:,...
                        test_indx_p),2,'omitnan');
                    fold_indx_r = fold_indx_r +n_cv_r;
                catch
                end
                MetaPI(l).PI.data(i).Name = strjoin({unique_groups{i} 'Progressors'},'_');           % First rows of data structure contain the progressors
                MetaPI(l).PI.data(i).dataTime = time;
                MetaPI(l).PI.data(i).dataValue(:,1) = train_i_progressors;
                MetaPI(l).PI.data(i).testValue(:,1) = test_i_progressors;
                MetaPI(l).PI.data(i).Group = unique_groups(i);
                MetaPI(l).PI.data(i).censoring = 'right';
                
                MetaPI(l).PI.data(i+length(unique_groups)).Name = strjoin({unique_groups{i} 'Responders'},'_');
                MetaPI(l).PI.data(i+length(unique_groups)).dataTime = time;
                MetaPI(l).PI.data(i+length(unique_groups)).dataValue(:,1) = train_i_responders;
                MetaPI(l).PI.data(i+length(unique_groups)).testValue(:,1) = test_i_responders;
                MetaPI(l).PI.data(i+length(unique_groups)).Group = unique_groups(i);
                MetaPI(l).PI.data(i+length(unique_groups)).censoring = 'left';   
            end
        else
            for l=1:p.n_folds
                MetaPI(l).PI.data(i).Name = strjoin({unique_groups{i} 'Progressors'},'_');         % First rows of data structure contain the progressors
                MetaPI(l).PI.data(i).dataTime = time;
                MetaPI(l).PI.data(i).dataValue(:,j) = mean(mat_j,2);
                MetaPI(l).PI.data(i).testValue(:,j) = mean(mat_j,2);
                MetaPI(l).PI.data(i).Group = unique_groups(i);
                MetaPI(l).PI.data(i).censoring = 'right';    
                 MetaPI(l).PI.data(i+length(unique_groups)).dataValue(:,j) = mean(mat_j,2);
                MetaPI(l).PI.data(i+length(unique_groups)).testValue(:,j) = mean(mat_j,2);
            end  
        end
    end 
end

else
%     for i=1:length(groups)
%         
%         MetaPI(l).PI.data(i).Name = strjoin({groups{i} 'Progressors'},'_');           % First rows of data structure contain the progressors
%         MetaPI(l).PI.data(i).dataTime = time;
%         MetaPI(l).PI.data(i).dataValue(:,1) = train_i_progressors;
%         MetaPI(l).PI.data(i).testValue(:,1) = test_i_progressors;
%         MetaPI(l).PI.data(i).Group = unique_groups(i);
%         MetaPI(l).PI.data(i).censoring = 'right';
%         
%         MetaPI(l).PI.data(i+length(unique_groups)).Name = strjoin({unique_groups{i} 'Responders'},'_');
%         MetaPI(l).PI.data(i+length(unique_groups)).dataTime = time;
%         MetaPI(l).PI.data(i+length(unique_groups)).dataValue(:,1) = train_i_responders;
%         MetaPI(l).PI.data(i+length(unique_groups)).testValue(:,1) = test_i_responders;
%         MetaPI(l).PI.data(i+length(unique_groups)).Group = unique_groups(i);
%         MetaPI(l).PI.data(i+length(unique_groups)).censoring = 'left';
%     end
end

for l = 1:p.n_folds
    nanIndx = arrayfun(@(x) all(all(isnan(x.dataValue))), MetaPI(l).PI.data);                % Identify and eliminate individuals where all obs are nan
    MetaPI(l).PI.data = MetaPI(l).PI.data(:,~nanIndx)';
    MetaPI(l).PI.data = MetaPI(l).PI.data(ismember([MetaPI(l).PI.data(:).Group], groups_subset));             % Subset data to use only the input groups
    colors = table2cell(table(linspecer(length( MetaPI(l).PI.data))));
    [MetaPI(l).PI.data(1:end).colors] = colors{:,:};
    try
        MetaPI(l).PI.data(ismember([MetaPI(l).PI.data(:).Group], 'MOC1_Control_Mean')).Group ...
            = {'MOC1_Control'};
    catch
    end
    
    try
        MetaPI(l).PI.data(ismember([MetaPI(l).PI.data(:).Group], 'MOC2_Control_Mean')).Group ...
            = {'MOC2_Control'};
    catch
    end
    zero_indx = arrayfun(@(x) sum(x.dataValue(:,1)==0,2),MetaPI(l).PI.data,...
        'UniformOutput', false);                                                % Identify indexes of tumor volume equal to 0
    [MetaPI(l).PI.data(1:end).zero_indx] = zero_indx{:,:};                                % Add index array to data array
    dataValue = arrayfun(@(x) x.dataValue(~x.zero_indx,:), MetaPI(l).PI.data,...          % Subset out the zero values and their respective time points
        'UniformOutput', false);
    dataTime = arrayfun(@(x) x.dataTime(~x.zero_indx,:), MetaPI(l).PI.data,...
        'UniformOutput', false);
    [MetaPI(l).PI.data(1:end).dataValue] = dataValue{:,:};
    [MetaPI(l).PI.data(1:end).dataTime] = dataTime{:,:};
    MetaPI(l).PI.n_data=sum(cellfun(@(x)sum(sum(~isnan(x))),{MetaPI(l).PI.data.dataValue},...
        'UniformOutput',true));
end
%% Merge phenotypes
if p.mergePhenotypes
    for l=1:p.n_folds
        name = ([MetaPI(l).PI.data(:).Group]);                                            % Extract names from the groups
        name_indx = cellfun(@(x)regexp(x,'Control'),name,'UniformOutput',false);% Identify  names with 'Control' in it
        name_indx = cellfun(@(x)~isempty(x),name_indx);                         % Obtain indexes for those names corresponding to control groups
        group_i = unique(name(name_indx));                                      % Get unique group names that correspond to the datasets to merge
        data = MetaPI(l).PI.data;
        MetaPI(l).PI.data = data(~name_indx);
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
            MetaPI(l).PI.data(end+1).Name = name_i{1};
            MetaPI(l).PI.data(end).dataTime = time_i;
            MetaPI(l).PI.data(end).dataValue = dataset;
            MetaPI(l).PI.data(end).Group = group_i(i);
            MetaPI(l).PI.data(end).censoring = {data(find(indx_i,1)).censoring};
        end
    end
end

%% Plot data
ncol = ceil(sqrt(length(stateVar)));
nrow = ceil(length(stateVar)/ncol);
figure;

for l=1:p.n_folds
    colors = table2cell(table(linspecer(length(MetaPI(l).PI.data))));
    [MetaPI(l).PI.data(1:end).colors] = colors{:,:};
    for i = 1:length(stateVar)
        subplot(nrow,ncol,i)
        hold on
        arrayfun(@(x)plot(x.dataTime, x.dataValue(:,i),'Color',...
            x.colors,'Marker','*'),MetaPI(l).PI.data)
        legend({MetaPI(l).PI.data(:).Name},'Interpreter', 'none')
        title(stateVar(i))
        set(gca,'YScale', 'log')
    end
    MetaPI(l).PI.tspan = unique(cat(1,MetaPI(l).PI.data(:).dataTime));
end
return
