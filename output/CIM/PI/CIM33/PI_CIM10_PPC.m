
%% PKPD Human efficacy general setup
% Search paths
clear all
warning off
if ispc
    addpath(genpath('\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox'))
    cd('\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\output')
    out = sbioloadproject('\Users\jmten\OneDrive\Dokumente\GitHub\sbio-projects\CIM_10_1_PKPD.sbproj');
    data_ext = {'\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\PI_Clavijo_2.mat'};
    data_ext1 = {'\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\PI_Morisada_3.mat'};
else
    addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
    cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output')
    out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/sbio-projects/CIM_10_1_PK.sbproj');
    data_ext = {'/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo_2.mat'};
    data_ext1={'/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Morisada_3.mat'};
end
%%
PI_MOC1 = load(strjoin({cd 'CIM/PI/CIM33/PI_CIM10_MOC1_PKPD_1_8_kinCD8.mat'},'/'));
PI_MOC2 = load(strjoin({cd 'CIM/PI/CIM33/PI_CIM10_MOC2_PKPD_1_8.mat'},'/'));
PI_MC38 = load(strjoin({cd 'CIM/PI/CIM33/PI_CIM10_MC38_PKPD_1_8.mat'},'/'));
PI_MOC1=PI_MOC1.PI;
PI_MOC2=PI_MOC2.PI;
PI_MC38=PI_MC38.PI;
PI = [];
names = {PI_MOC1.output(1:end).Name PI_MOC2.output(1:end).Name PI_MC38.output(1:end).Name};
TV = {PI_MOC1.output(1:end).TV PI_MOC2.output(1:end).TV PI_MC38.output(1:end).TV};
CD8 = {PI_MOC1.output(1:end).CD8 PI_MOC2.output(1:end).CD8 PI_MC38.output(1:end).CD8};
Treg = {PI_MOC1.output(1:end).Treg PI_MOC2.output(1:end).Treg PI_MC38.output(1:end).Treg};
MDSC = {PI_MOC1.output(1:end).MDSC PI_MOC2.output(1:end).MDSC PI_MC38.output(1:end).MDSC};
DC = {PI_MOC1.output(1:end).DCm PI_MC38.output(1:end).DCm};
CD8_E = {PI_MOC1.output(1:end).CD8_E  PI_MC38.output(1:end).CD8_E};
PDL1_I = {PI_MOC1.output(1:end).PDL1_I  PI_MC38.output(1:end).PDL1_I};
PDL1_T = {PI_MOC1.output(1:end).PDL1_T  PI_MC38.output(1:end).PDL1_T};

PI.output=[];
[PI.output(1:length(names)).Name]=names{:,:};
[PI.output(1:length(names)).TV]=TV{:,:};
[PI.output(1:length(names)).CD8]=CD8{:,:};
[PI.output(1:length(names)).Treg]=Treg{:,:};
[PI.output(1:length(names)).MDSC]=MDSC{:,:};

[PI.output([1:9 14:15]).DCm]=DC{:,:};
[PI.output([1:9 14:15]).CD8_E]=CD8_E{:,:};
[PI.output([1:9 14:15]).PDL1_I]=PDL1_I{:,:};
[PI.output([1:9 14:15]).PDL1_T]=PDL1_T{:,:};
%% Add treatments, cell and response factors
treatments = cellfun(@(x)strrep(x,'MOC1_',''),names,'UniformOutput',false);
treatments = cellfun(@(x)strrep(x,'MOC2_',''),treatments,'UniformOutput',false);
treatments = cellfun(@(x)strrep(x,'MC38_',''),treatments,'UniformOutput',false);
treatments = cellfun(@(x)strrep(x,'_Progressor',''),treatments,'UniformOutput',false);
treatments = cellfun(@(x)strrep(x,'_Responder',''),treatments,'UniformOutput',false);
treatments = cellfun(@(x)strrep(x,'_Immune',''),treatments,'UniformOutput',false);
treatments = cellfun(@(x)strrep(x,'_Mean',''),treatments,'UniformOutput',false);

cell=cellfun(@(x) strsplit(x,'_'),names,'UniformOutput',false);
cell=cellfun(@(x) x{1},cell,'UniformOutput',false);
response=cellfun(@(x) strsplit(x,'_'),names,'UniformOutput',false);
response=cellfun(@(x) x{end},response,'UniformOutput',false);

[PI.output(1:end).Treatment]=treatments{:,:};
[PI.output(1:end).Cell]=cell{:,:};
[PI.output(1:end).Response]=response{:,:};

%% Add CI tables

TV = {PI_MOC1.CI(1:end).TV PI_MOC2.CI(1:end).TV PI_MC38.CI(1:end).TV};
CD8 = {PI_MOC1.CI(1:end).CD8 PI_MOC2.CI(1:end).CD8 PI_MC38.CI(1:end).CD8};
Treg = {PI_MOC1.CI(1:end).Treg PI_MOC2.CI(1:end).Treg PI_MC38.CI(1:end).Treg};
MDSC = {PI_MOC1.CI(1:end).MDSC PI_MOC2.CI(1:end).MDSC PI_MC38.CI(1:end).MDSC};
DC = {PI_MOC1.CI(1:end).DCm PI_MC38.CI(1:end).DCm};
CD8_E = {PI_MOC1.CI(1:end).CD8_E  PI_MC38.CI(1:end).CD8_E};
PDL1_I = {PI_MOC1.CI(1:end).PDL1_I  PI_MC38.CI(1:end).PDL1_I};
PDL1_T = {PI_MOC1.CI(1:end).PDL1_T  PI_MC38.CI(1:end).PDL1_T};

PI.CI=[];
[PI.CI(1:length(names)).Name]=names{:,:};
[PI.CI(1:length(names)).TV]=TV{:,:};
[PI.CI(1:length(names)).CD8]=CD8{:,:};
[PI.CI(1:length(names)).Treg]=Treg{:,:};
[PI.CI(1:length(names)).MDSC]=MDSC{:,:};

[PI.CI([1:9 14:15]).DCm]=DC{:,:};
[PI.CI([1:9 14:15]).CD8_E]=CD8_E{:,:};
[PI.CI([1:9 14:15]).PDL1_I]=PDL1_I{:,:};
[PI.CI([1:9 14:15]).PDL1_T]=PDL1_T{:,:};

[PI.CI(1:end).Treatment]=treatments{:,:};
[PI.CI(1:end).Cell]=cell{:,:};
[PI.CI(1:end).Response]=response{:,:};

%% Add individual data
PI.IndivData = [PI_MOC1.IndivData; PI_MOC2.IndivData; PI_MC38.IndivData];

observables = [repelem({PI_MOC1.observablesFields},length(PI_MOC1.IndivData),1);...
    repelem({PI_MOC2.observablesFields},length(PI_MOC2.IndivData),1);...
    repelem({PI_MC38.observablesFields},length(PI_MC38.IndivData),1)];
[PI.IndivData(1:end).observables] = observables{:,:};
PI.tspan = 1:100;
PI.variableUnits = PI_MOC1.variableUnits;
%% plot
plotPosteriorPredictions2(PI,{'TV' 'CD8' 'Treg' 'MDSC' 'DCm' 'CD8_E' 'PDL1_I' 'PDL1_T'},'treatments',{'Control','antiPDL1','antiCTLA4','antiCTLA4_antiPDL1'})

h=gcf;
h.Children.Children(2).XLabel.String = 'Time[days]';
h.Children.Children(2).XLabel.FontWeight = 'bold';
h.Children.Children(2).XLabel.FontSize = 12;
h.Children.Children(2).XTick = 0:5:35;
h.Children.Children(2).XTickLabel =cellfun(@(x)num2str(x),num2cell(0:5:35),'UniformOutput',false);

h.Children.Children(6).XLabel.String = 'Time[days]';
h.Children.Children(6).XLabel.FontWeight = 'bold';
h.Children.Children(6).XLabel.FontSize = 12;
h.Children.Children(6).XTick = 0:5:30;
h.Children.Children(6).XTickLabel =cellfun(@(x)num2str(x),num2cell(0:5:30),'UniformOutput',false);
