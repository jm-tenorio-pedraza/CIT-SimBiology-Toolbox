function [PI,doses] = getDataSets(dataset_file_ext,varargin)
p = inputParser;
p.addParameter('zeroAction','omit');
p.parse(varargin{:});
p=p.Results;

data = cellfun(@(x)readtable(x),dataset_file_ext,'UniformOutput',false);

varNames = cellfun(@(x) x.Properties.VariableNames,data,'UniformOutput',false);
varNames = cat(2,varNames{:});
varNames = unique(varNames,'stable');
varNames = varNames(~cellfun(@(x)strcmp(x(1:2),'SD'),varNames));            % Extract unique variable names from each data table of the datasets
varNames = setdiff(varNames, {'Time_hour_s__' 'Group' 'Dose_mg_kg_'});


groups = cellfun(@(x)unique(x.Group),data,'UniformOutput',false);
groups = cat(1,groups{:,:});                                                  % Extract the groups from each dataset

meanIndx = cellfun(@(x) cellfun(@(y)strcmp(y(1:2),'SD'),...
    x.Properties.VariableNames,'UniformOutput',true),data,...
    'UniformOutput', false);                                                % Extract column indexes where time points and mean data is stored

% meanIndx = cat(2,meanIndx{:,:});

nVar = length(varNames);                                                  % Set the number of columns that the common dataset will have
dataValue = cellfun(@(x)nan(size(x,1),nVar),data,'UniformOutput',false);    % Allocate space for each data matrix corresponding to each dataset
PI = [];                                                                    % Create PI structure
[PI.data(1:length(groups)).Name] = groups{:,:};                                 % Add one group entry to the PI structure for each considered condition in each dataset

[PI.data(1:length(groups)).Group] = groups{:,:};                                 % Add one group entry to the PI structure for each considered condition in each dataset

for i=1:length(data)
    data_i = data{i};
    [varIndx,dataIndx] = ismember(varNames,data_i.Properties.VariableNames);
    dataValue{i}(:,varIndx) = data_i{:,setdiff(dataIndx,0,'stable')};
end
indx=1;
% Add data values and time points
for i = 1:length(data)
    data_i = data{i};
    groups_i = unique(data_i.Group);
    for j = 1:length(groups_i)
        valueIndx = ismember(data_i.Group,groups_i(j));
        dataValue_i = dataValue{i}(valueIndx,:);
        dataTime_i = data_i{valueIndx,1};
        if strcmp(p.zeroAction,'omit')
            zeroIndx = ~any(dataValue_i==0,2);
        else
            zeroIndx = dataValue==dataValue_i;
        end
        PI.data(indx).dataTime = dataTime_i(zeroIndx,:);
        PI.data(indx).dataValue = dataValue_i(zeroIndx,:);
        indx = indx+1;
    end
end
     
% Add dosage

[PI,doses] = getDoses(PI);
doses = doses';
for i=1:length(doses)
    doses{i}.Properties.VariableUnits = {'hour' 'milligram' 'milligram/second'};
end

[PI.data(1:end).doses] = doses{:,:};

PI.stateVar = varNames;
PI.n_data=sum(cellfun(@(x)sum(sum(~isnan(x))),{PI.data.dataValue},'UniformOutput',true));
PI.data = PI.data';
PI.tspan = unique(cat(1,PI.data(:).dataTime));
return