function datasetSimBio=combineDatasets2(data_table_1,struct,varargin)
p = inputParser;
p.addParameter('dose',0,@isnumeric);
p.parse(varargin{:});
p = p.Results;

% Choose the table with the largest number of rows
name_fields = fieldnames(struct);
data_s = structfun(@(x) reshape(x{:,1:end},[],1),struct,'UniformOutput',false);

group = structFun(@(x)x.Properties.VariableNames,struct);
name_field_1 = char(name_fields(1));
group = repelem(group.(name_field_1)(1,2:end),4);
unique_groups = unique(group);

datasetSimBio=[];

% Create table with data from the struct 
data_table = struct2table(data_s);
data_table.Group = group';

% data_table.Dose = repelem(p.dose,size(data_table,1),1);

% Merge the two tables
for i = 1:length(group)
    group_i = unique_groups(i);                                             % Select group
    dataindx = ismember(group, group_i);                                    % Select indexes of new data corresponding to group
    tableindx = ismember(data_table_1.Group,group_i);                       % Select indexes of data corresponding to group
    t_merged = uniqe([data_table.Time data_table_1.Time]);                  % Time point vector merged
    dataTindx = ismember(t_merged, data_table.Time);                        % Select indexes of new data time points corresponding to merged vector
    tableTindx = ismember(t_merged, data_table_1.Time);                     % Select indexes of old data time points corresponding to merged vector
    datasetSimBio.Name(i) = group_i; 
    datasetSimBio.dataTime(i) = sort(t_merged(or(dataTindx, tableTindx)));  % Set merged time points as the first column of the 
    dataValue( = NaN(length(datasetSimBio.dataTime(i) ), size(data_table_1)
    [data_table_1(tableindx,2)...
        data_table(dataindx,2:end-1)];
    datasetSimBio(i,:) = datasetSimBio_i;
end

return
