%% Morisada et al, 2017 data processing

clear all
addpath(genpath("/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/code"))
% addpath(genpath("/Users/migueltenorio/Documents/MATLAB/Parameter Identifications/aux functions"))

%% Fig 6
MOC1_treatment_names={'Control' 'antiPD1' '8Gyx2' '2Gyx10' 'antiPD1_8Gyx2' 'antiPD1_2Gyx2'};
MOC1_names=cellfun(@(y)cellfun(@(x)strjoin({'MOC1' y num2str(x)},'_'),num2cell(1:10),'UniformOutput', false),...
    MOC1_treatment_names,'UniformOutput',false);

MOC1_names=['Time' MOC1_names{:,:}];
Morisada_MOC1_Tumor_Fig6 = importPrismFile("/Users/migueltenorio/Documents/Data/CSV/Morisada_2017/MOC1_Fig6.csv",...
    [2 Inf], MOC1_names,0);

MC38_treatment_names={'Control' 'antiPD1' '8Gyx2' '2Gyx10' 'antiPD1_8Gyx2' 'antiPD1_2Gyx2'};
MC38_names=cellfun(@(y)cellfun(@(x)strjoin({'MC38' y num2str(x)},'_'),num2cell(1:10),'UniformOutput', false),...
    MC38_treatment_names,'UniformOutput',false);

MC38_names=['Time' MC38_names{:,:}];
Morisada_MC38_Tumor_Fig6 = importPrismFile("/Users/migueltenorio/Documents/Data/CSV/Morisada_2017/MC38_Fig6.csv",...
    [2 Inf], MC38_names,0);

%% Data array 
% Tumor volume table
TumorName=[reshape(repmat(Morisada_MOC1_Tumor_Fig6.Properties.VariableNames(2:end),...
    size(Morisada_MOC1_Tumor_Fig6,1),1),[],1);
    reshape(repmat(Morisada_MC38_Tumor_Fig6.Properties.VariableNames(2:end),...
    size(Morisada_MC38_Tumor_Fig6,1),1),[],1)];

Morisada_TumorTime = [repmat(Morisada_MOC1_Tumor_Fig6.Time,size(Morisada_MOC1_Tumor_Fig6,2)-1,1);
    repmat(Morisada_MC38_Tumor_Fig6.Time,size(Morisada_MC38_Tumor_Fig6,2)-1,1);];

MorisadaTumor_table=table(Morisada_TumorTime);
MorisadaTumor_table.Properties.VariableNames={'Time'};
MorisadaTumor_table(:,2)=table([reshape(Morisada_MOC1_Tumor_Fig6{:,2:end},[],1);
    reshape(Morisada_MC38_Tumor_Fig6{:,2:end},[],1)]);
MorisadaTumor_table.Properties.VariableNames(2)={'Tumor'};
MorisadaTumor_table(:,3)=reshape(TumorName,[],1);
MorisadaTumor_table.Properties.VariableNames(end)={'Name'};


time = unique([MorisadaTumor_table.Time]);
name = unique([MorisadaTumor_table.Name],'stable');

% Strucutre array
simBioData = [];
[simBioData(1:length(name)).Name] = name{:,:};
for i = 1:length(simBioData)
    name_i = name(i);
     
    tumorname_indx = ismember(MorisadaTumor_table.Name, name_i);
    tumortime_indx = ismember(time, MorisadaTumor_table.Time(tumorname_indx));
    
    array_i = NaN(length(time),1);
    
    array_i(tumortime_indx,1) = MorisadaTumor_table{tumorname_indx,2};
   
    
    nan_indx = all(isnan(array_i(:,1:end)),2);
    simBioData(i).dataTime = time(~nan_indx,1);
    simBioData(i).dataValue = array_i(~nan_indx,[1 2:2:end]);

end
dataIndx = arrayfun(@(x) ~isempty(x.dataValue), simBioData);
simBioData= simBioData(dataIndx);

group = cellfun(@(x) x(1:regexp(x,'_\d*$')-1),{simBioData.Name},'UniformOutput', false)';
empty_indx = cellfun(@(x) isempty(x), group, 'UniformOutput', true);
group(empty_indx) = name(empty_indx);

[simBioData(1:end).Group]=group{:,:};
PI.data = simBioData;
PI.stateVar = {'Tumor'};

%% Plotting results
groupNames = unique({PI.data(:).Group}, 'stable');
ncol = ceil(sqrt(length(groupNames)));
nrow = ceil(length(groupNames)/ncol);
for i=1:length(groupNames)
    subplot(2,6,i)
    hold on
    data_indx = ismember({PI.data(:).Group}, groupNames(i));
    data_i = PI.data(data_indx);
    arrayfun(@(x) plot(x.dataTime, x.dataValue),data_i)
    title(groupNames(i), 'interpreter', 'none')
end
%%
save('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Tumor_Morisada.mat', 'PI')