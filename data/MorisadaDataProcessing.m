%% 
Morisada_Immune=struct('CD8',Morisada_MOC1_CD8,'Treg', Morisada_MOC1_Treg,...
    'GMDSC' ,Morisada_MOC1_GMDSC,'DC',Morisada_MOC1_DC, 'Myeloid_PDL1_Rel',Morisada_MOC1_Myeloid_PDL1,...
    'Tumor_PDL1_Rel', Morisada_MOC1_Tumor_PDL1);

% Row names of data structure
ImmuneName=struct2cell(structfun(@(x) repmat(x.Properties.VariableNames(2:end),size(x{:,2:end},1),1),...
    Morisada_Immune,'UniformOutput', false));
ImmuneName=reshape(ImmuneName{1},[],1);
% Time vector
Morisada_ImmuneTime=struct2cell(structfun(@(x) repmat(x{:,1},size(x{:,2:end},2),1),...
    Morisada_Immune,'UniformOutput', false));
Morisada_ImmuneTime=cell2mat(Morisada_ImmuneTime(1));

% Variable matrix
Morisada_Immune=structfun(@(x) reshape(x{:,2:end},[],1),Morisada_Immune,'UniformOutput', false);
Morisada_Immune=struct2table(Morisada_Immune);
Morisada_Immune{:,1:4}=(Morisada_Immune{:,1:4}/1e4*100);
Morisada_Immune{:,5:6} = Morisada_Immune{:,5:6}./...
    repelem(Morisada_Immune{4:4:end,5:6},4,1);
Morisada_Immune{:,end+1:end+2} = Morisada_Immune{:,3:4}./...
    repelem(Morisada_Immune{4:4:end, 3:4},4,1);
Morisada_Immune.Properties.VariableNames(end-1:end) = {'GMDSC_Rel' 'DC_Rel'};
MorisadaImmune_table = table(Morisada_ImmuneTime,'VariableName', {'Time'});
MorisadaImmune_table(:,2:size(Morisada_Immune,2)+1)=Morisada_Immune;
MorisadaImmune_table.Properties.VariableNames(2:end)=Morisada_Immune.Properties.VariableNames;
MorisadaImmune_table(:,end+1)=ImmuneName;
MorisadaImmune_table.Properties.VariableNames(end)={'Name'};

% Tumor volume table
TumorName=repmat(Morisada_MOC1.Properties.VariableNames(2:end), size(Morisada_MOC1,1),1);
Morisada_TumorTime=repmat(Morisada_MOC1.Time,size(Morisada_MOC1,2)-1,1);
MorisadaTumor_table=table(Morisada_TumorTime);
MorisadaTumor_table.Properties.VariableNames={'Time'};
MorisadaTumor_table(:,2)=table(reshape(Morisada_MOC1{:,2:end},[],1));
MorisadaTumor_table.Properties.VariableNames(2)={'Tumor'};
MorisadaTumor_table(:,3)=reshape(TumorName,[],1);
MorisadaTumor_table.Properties.VariableNames(end)={'Name'};

% Compare tables and merge
% Elements in the larger array that are also in the smaller one
nameindx1=ismember(MorisadaTumor_table.Name,MorisadaImmune_table.Name);
[rowIndx1] = and(nameindx1,...
    ismember(MorisadaTumor_table.Time, MorisadaImmune_table.Time));
% Elements in the smaller array that are not in the larger one
nameindx2 = ismember(MorisadaImmune_table.Name, MorisadaTumor_table.Name);
[rowIndx2] = and(nameindx2,...
    ~ismember(MorisadaImmune_table.Time, MorisadaTumor_table.Time));

time = unique([MorisadaTumor_table.Time; MorisadaImmune_table.Time]);
name = unique([MorisadaTumor_table.Name; MorisadaImmune_table.Name]);
% Strucutre array
simBioData = [];
[simBioData(1:length(name)).Name] = name{:,:};
for i = 1:length(simBioData)
    name_i = name(i);
    immunename_indx = ismember(MorisadaImmune_table.Name, name_i);
    immunetime_indx = ismember(time, MorisadaImmune_table.Time);
    
    tumorname_indx = ismember(MorisadaTumor_table.Name, name_i);
    tumortime_indx = ismember(time, MorisadaTumor_table.Time);
    
    array_i = NaN(length(time),1+size(MorisadaImmune_table,2)-2);
    try
       array_i(immunetime_indx,2:end) = MorisadaImmune_table{immunename_indx,2:end-1};
    catch
    end
    array_i(tumortime_indx,1) = MorisadaTumor_table{tumorname_indx,2};
    nan_indx = all(isnan(array_i(:,1:end)),2);
    simBioData(i).dataTime = time(~nan_indx,1);
    simBioData(i).dataValue=array_i(~nan_indx,:);
end

group = cellfun(@(x) x(1:regexp(x,'_\d*$')-1),{simBioData.Name},'UniformOutput', false)';
[simBioData(1:end).Group]=group{:,:};
PI.data = simBioData;
PI.stateVar = {'Tumor' 'CD8' 'Treg' 'GMDSC' 'DC' 'Myeloid_PDL1_Rel' 'Tumor_PDL1_Rel' 'GMDSC_Rel' 'DC_Rel'};

%% Adding logit-transforms
stateVar = {'CD8' 'Treg' 'GMDSC' 'DC'};
[arrayindx,stateindx] = ismember(stateVar,PI.stateVar);
stateindx = stateindx(stateindx~=0);
logit_transforms = arrayfun(@(x) x.dataValue(:,stateindx)./(100-x.dataValue(:,stateindx)),...
    PI.data,'UniformOutput', false);
for i=1:length(logit_transforms)
    PI.data(i).dataValue(:,length(PI.stateVar)+1:length(PI.stateVar)+length(stateindx)) = logit_transforms{1,i};
end
PI.stateVar = [PI.stateVar {'CD8_logit' 'Treg_logit' 'GMDSC_logit' 'DC_logit'}];
%%

% Save output
save('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/Morisada_Data.mat', 'simBioData')
save('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Morisada.mat', 'PI')