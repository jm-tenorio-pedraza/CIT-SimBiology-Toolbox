addpath(genpath("/Users/migueltenorio/Documents/MATLAB/SimBiology"))
addpath(genpath("/Users/migueltenorio/Documents/MATLAB/Parameter Identifications/aux functions"))
%% Import Prism files
% Clavijo et al, 2017
% MOC1
MOC1_treatment_names={'Control' 'antiPDL1' 'antiCTLA4' 'antiLy6G' 'antiLy6G_antiPDL1'...
    'antiCTLA4_antiLy6G' 'antiCTLA4_antiPDL1'};
MOC1_names=cellfun(@(y)cellfun(@(x)strjoin({'MOC1' y num2str(x)},'_'),num2cell(1:12),'UniformOutput', false),...
    MOC1_treatment_names,'UniformOutput',false);

MOC1_names=['Time' MOC1_names{:,:}];
Clavijo_MOC1_Tumor = importPrismFile("/Users/migueltenorio/Documents/Data/CSV/Clavijo_2017/MOC1 tumors_edit.csv",...
    [2 Inf], MOC1_names,0);
Clavijo_MOC1_Tumor = Clavijo_MOC1_Tumor(1:39,:);
Clavijo_Immune = readtable('/Users/migueltenorio/Documents/Data/Clavijo_2017_Immune_SimBio.xlsx');
% MOC2
MOC2_treatment_names={'Control' 'antiLy6G' 'antiPDL1' 'antiCTLA4'  'antiLy6G_antiPDL1'...
    'antiCTLA4_antiLy6G'};
MOC2_names=cellfun(@(y)cellfun(@(x)strjoin({'MOC2' y num2str(x)},'_'),num2cell(1:8),...
    'UniformOutput', false),MOC2_treatment_names,'UniformOutput',false);

MOC2_names=['Time' MOC2_names{:,:}];
Clavijo_MOC2_Tumor = importPrismFile("/Users/migueltenorio/Documents/Data/CSV/Clavijo_2017/MOC2 tumors grouped.csv",...
    [2 Inf], MOC2_names,0);

%% 
% Tumor volume table
TumorName=[reshape(repmat(Clavijo_MOC1_Tumor.Properties.VariableNames(2:end),...
    size(Clavijo_MOC1_Tumor,1),1),[],1); reshape(repmat(Clavijo_MOC2_Tumor.Properties.VariableNames(2:end),...
    size(Clavijo_MOC2_Tumor,1),1),[],1)];
Clavijo_TumorTime=[repmat(Clavijo_MOC1_Tumor.Time,size(Clavijo_MOC1_Tumor,2)-1,1);...
    repmat(Clavijo_MOC2_Tumor.Time,size(Clavijo_MOC2_Tumor,2)-1,1)];

ClavijoTumor_table=table(Clavijo_TumorTime);
ClavijoTumor_table.Properties.VariableNames={'Time'};
ClavijoTumor_table(:,2)=table([reshape(Clavijo_MOC1_Tumor{:,2:end},[],1); reshape(Clavijo_MOC2_Tumor{:,2:end},[],1)]);
ClavijoTumor_table.Properties.VariableNames(2)={'Tumor'};
ClavijoTumor_table(:,3)=reshape(TumorName,[],1);
ClavijoTumor_table.Properties.VariableNames(end)={'Name'};


time = unique([ClavijoTumor_table.Time; Clavijo_Immune.Time]);
name = unique([ClavijoTumor_table.Name; Clavijo_Immune.group]);
% Strucutre array
simBioData = [];
[simBioData(1:length(name)).Name] = name{:,:};
for i = 1:length(simBioData)
    name_i = name(i);
    immunename_indx = ismember(Clavijo_Immune.group, name_i);
    immunetime_indx = ismember(time, Clavijo_Immune.Time(immunename_indx));
    
    tumorname_indx = ismember(ClavijoTumor_table.Name, name_i);
    tumortime_indx = ismember(time, ClavijoTumor_table.Time(tumorname_indx));
    
    array_i = NaN(length(time),1+size(Clavijo_Immune,2)-3);
    try
       array_i(immunetime_indx,2:end) = Clavijo_Immune{immunename_indx,2:end-2};
    catch
       array_i(tumortime_indx,1) = ClavijoTumor_table{tumorname_indx,2};

    end
    try
       array_i(tumortime_indx,1) = ClavijoTumor_table{tumorname_indx,2};

    catch
       array_i(immunetime_indx,2:end) = Clavijo_Immune{immunename_indx,2:end-2};

    end
    nan_indx = all(isnan(array_i(:,1:end)),2);
    simBioData(i).dataTime = time(~nan_indx,1);
    simBioData(i).dataValue=array_i(~nan_indx,:);
end
dataIndx = arrayfun(@(x) ~isempty(x.dataValue), simBioData);
simBioData= simBioData(dataIndx);

group = cellfun(@(x) x(1:regexp(x,'_\d*$')-1),{simBioData.Name},'UniformOutput', false)';
empty_indx = cellfun(@(x) isempty(x), group, 'UniformOutput', true);
group(empty_indx) = name(empty_indx);

[simBioData(1:end).Group]=group{:,:};
PI.data = simBioData;
PI.stateVar = {'Tumor' 'CD8' 'CD107a' 'Treg' 'DC' 'GMDSC' 'Tumor_PDL1_Rel'...
    'Myeloid_PDL1_Rel' 'DC_Rel'  'GMDSC_Rel' 'CD8_logit' 'CD107a_logit'...
    'Treg_logit' 'DC_logit' 'GMDSC_logit'};
%% Save output
save('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/Clavijo_Data.mat', 'simBioData')
save('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat', 'PI')




