if ispc
    wd='\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\';
else
    wd='/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox\';
end


%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: C:\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\Sun_2019_Fig7f.csv
%
% Auto-generated by MATLAB on 16-Sep-2023 05:40:20

%% Setup the Import Options and import the data Sun_2019_Fig7
opts = delimitedTextImportOptions("NumVariables", 10);

% Specify range and delimiter
opts.DataLines = [3, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["MOC1_Control", "VarName2", "MOC1_1", "VarName4", "MOC1_2", "VarName6", "MOC1_5", "VarName8", "MOC1_10", "VarName10"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
Sun2019Fig7f = readtable("C:\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\Sun_2019_Fig7f.csv", opts);
%% Setup the Import Options and import the data Sun_2019_fig1
opts = delimitedTextImportOptions("NumVariables", 8);

% Specify range and delimiter
opts.DataLines = [3, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["MOC1_Control", "VarName2", "MOC1_TIL_10", "VarName4", "MOC1_TIL_GMDSC", "VarName6", "MOC1_GMDSC", "VarName8"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
Sun2019Fig1c = readtable("C:\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\Sun_2019_Fig1c.csv", opts);





%%
SI=[];
dataNames={'Sun2019_Fig7f_Control'; "Sun2019_Fig7f_Tumor:TIL_1:1"; "Sun2019_Fig7f_Tumor:TIL_1:2";...
    "Sun2019_Fig7f_Tumor:TIL_1:5"; "Sun2019_Fig7F_Tumor:TIL_1:10"; 'Sun2019_Fig1c_Control';...
    'Sun2019_Fig1c_Tumor:TIL_10:1'; 'Sun2019_Fig1c_Tumor:TIL:GMDSC_1:10:30'; 'Sun2019_Fig1c_Tumor:GMDSC_1:30'};
[SI.data(1:9).Name]=dataNames{:,:};

%% Dataset-specific data
% Sun 2019, Fig7f
dataTimes = round(Sun2019Fig7f{:,1:2:end}*24);
dataTimes=mat2cell(dataTimes', repelem(1,size(dataTimes',1)));
[SI.data(1:5).dataTime]=dataTimes{:,:};

dataValue=Sun2019Fig7f{:,2:2:end};
dataValue(1,:)=repelem(1,5);
dataValue=mat2cell(dataValue',repelem(1,size(dataValue',1)));
dataValue=cellfun(@(x)x',dataValue,'UniformOutput', false);
[SI.data(1:5).dataValue]=dataValue{:,:};


SD=NaN(5,size(Sun2019Fig7f,1));
SD=mat2cell(SD',repelem(1,size(SD',1)));
[SI.data(1:5).SD]=SD{:,:};

% Sun 2019, Fig1c
dataTimes = round(Sun2019Fig1c{:,1:2:end}*24);
dataTimes=mat2cell(dataTimes', repelem(1,size(dataTimes',1)));
[SI.data(6:9).dataTime]=dataTimes{:,:};

dataValue=Sun2019Fig1c{:,2:2:end};
dataValue(1,:)=repelem(1,4);
dataValue=mat2cell(dataValue',repelem(1,size(dataValue',1)));
dataValue=cellfun(@(x)x',dataValue,'UniformOutput', false);
[SI.data(6:9).dataValue]=dataValue{:,:};


SD=NaN(4,size(Sun2019Fig1c,1));
SD=mat2cell(SD',repelem(1,size(SD',1)));
[SI.data(6:9).SD]=SD{:,:};

%% Common dataset
groupNames={'MOC1_Sun2019_Fig7f' 'MOC1_Sun2019_Fig7f' 'MOC1_Sun2019_Fig7f'...
    'MOC1_Sun2019_Fig7f' 'MOC1_Sun2019_Fig7f' 'MOC1_Sun2019_Fig1c' 'MOC1_Sun2019_Fig1c'...
    'MOC1_Sun2019_Fig1c' 'MOC1_Sun2019_Fig1c'};
[SI.data(1:9).Group]=expNames{:,:};

expNames= {'MOC1_TIL_1:0' 'MOC1_TIL_1:1' 'MOC1_TIL_1:2' 'MOC1_TIL_1:5' 'MOC1_TIL_1:10' ....
     'MOC1_TIL_1:0'  'MOC1_TIL_1:10'  'MOC1_TIL_GMDSC_1:10:30'  'MOC1_GMDSC_1:30'};
[SI.data(1:9).Experiment]=expNames{:,:};

cell={'MOC1' 'MOC1' 'MOC1' 'MOC1' 'MOC1' 'MOC1' 'MOC1' 'MOC1' 'MOC1'};
[SI.data(1:9).Cell]=cell{:,:};

responseGroup={'Responder' 'Responder' 'Responder' 'Responder' 'Responder' ...
     'Responder'  'Responder'  'Responder'  'Responder'};
[SI.data(1:9).Response]=responseGroup{:,:};


datacount= {7 7 7 7 7 9 9 9 9 };
[SI.data(1:5).Count]=datacount{:,:};

% Add initial values for the three species considered
%    T_0    CD8_0   GMDSC_0
x0=[.01     0       0; ...
    .01     .01     0; ...
    .01     .02     0; ...
    .01     .05     0; ...
    .01     .1      0; ...
    .01     0       0; ...
    .01     .1      0; ...
    .01     .1      .3; ...
    .01     0       .3;];
SI.x_0=x0;
%% Clear temporary variables
clear opts
%% Save structure
save(strjoin({wd 'data/PI_CIM10_MOC1_InVitro.mat'},''), 'SI')