function dataset_SimBio = createDatasetSimBio(data_table, D, group_ID)
% data_table is a table with the time points in the first column and 
% all the measurements of the resonse variable in the rest of the columns
% D is structure with one field for each group:
% group_indx:   Column indexes of the group
% dosing_time:  Vector with the time points on which the dose was given
% doses:        Dose given, one for each dosing_time time point
% group_ID is a cell array with the column names in data_table

dataset_SimBio=array2table(NaN(numel(data_table(:,2:end)),3));
dataset_SimBio(:,4)={'NaN'};
main_indx=1;
for i=1:length(D)
dim_dosing_time=size(D(i).dosing_time);
dim_doses=size(D(i).doses);
dim_data_table=D(i).group_indx;
% Select time vector
data_time=data_table{:,1};
% Check dimension consistency
if dim_dosing_time(2)>dim_dosing_time(1)
    dosing_time=D(i).dosing_time';
else
    dosing_time=D(i).dosing_time;
end
if dim_doses(2)>dim_doses(1)
    doses=D(i).doses';
else
    doses=D(i).doses;
end
% Making vector of doses for each unique data point
data_time=sort(unique([data_time; dosing_time]));
dose_array=zeros(size(data_time));
dose_array(ismember(data_time,dosing_time),1)=doses;

dose_array=repmat(dose_array,length(dim_data_table),1);
% Expanding the group_id cell array
group_ID_i=repelem(group_ID(D(i).group_indx),...
    length(data_time));

% Expanding vector of responses to include dose time points
dataset_SimBio_i=NaN(length(data_time),length(dim_data_table));
for j=1:length(dim_data_table)
    data_indx=ismember(data_time,data_table{:,1});
    dataset_SimBio_i(data_indx,j)=data_table{:,dim_data_table(j)};
end
% Redimensionizing into 1 column
dataset_SimBio_i=reshape(dataset_SimBio_i,[],1);

% Input into main table
dim_i=size(dataset_SimBio_i);
% Time vector
dataset_SimBio(main_indx:main_indx+dim_i(1)-1,1)=...
    table(repmat(data_time,length(dim_data_table),1));
% Response variable
dataset_SimBio(main_indx:main_indx+dim_i(1)-1,2)=table(dataset_SimBio_i);
% Group
dataset_SimBio(main_indx:main_indx+dim_i(1)-1,4)=table([group_ID_i']);

% Dose
dataset_SimBio(main_indx:main_indx+dim_i(1)-1,3)=table(dose_array);
main_indx=main_indx+dim_i(1);
end