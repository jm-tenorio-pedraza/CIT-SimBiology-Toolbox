function datasetSimBio=combineDatasets(data_table_1,data_table_2)
% Merges two data tables with different number of variables by adding rows
% or matching existing rows according to group variable
dim_1=size(data_table_1);
dim_2=size(data_table_2);

% Choose the table with the largest number of rows
if dim_1(1)>dim_2(1)
    datasetSimBio=data_table_1;
    table_j=data_table_2;
    
else
    datasetSimBio=data_table_2;
    table_j=data_table_1;
    dim_1=size(data_table_2);
    dim_2=size(data_table_1);
end

% Expand datasetSimBio in columns to add the relevant species
datasetSimBio(:,end+1:end+dim_2(2)-3)={NaN};
datasetSimBio.Properties.VariableNames=[datasetSimBio.Properties.VariableNames(1:dim_1(2)),...
    table_j.Properties.VariableNames(2:end-2)];

% Variables in main table that are not in secondary table
[a,b]=(ismember(datasetSimBio.Properties.VariableNames,table_j.Properties.VariableNames));
[c,d]=(ismember(table_j.Properties.VariableNames,datasetSimBio.Properties.VariableNames));
table_j=table_j(:,b(b~=0));
table_i=datasetSimBio(1:dim_2(1),:);
try
    table_i(:,a)=table_j;
catch
    
end
table_i(:,~a)={NaN};


% Check if table_i has combination of time points and group that match the
% ones in table_2
i_indx=and(ismember(datasetSimBio.Time, table_j.Time),ismember(datasetSimBio.Group, table_j.Group));
j_indx=and(ismember(table_j.Time, datasetSimBio.Time),ismember(table_j.Group, datasetSimBio.Group));

try
    datasetSimBio(i_indx,:)=table_i;
catch
    datasetSimBio(end+1:end+dim_2(1),:)=table_i;
end
