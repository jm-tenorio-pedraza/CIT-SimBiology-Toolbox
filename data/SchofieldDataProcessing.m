% Search paths
clear all
warning off
if ispc
    addpath(genpath('\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox'))
    cd('\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data')

else
    addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
    cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data')

end
%% 
heading = {'Time' 'Concentration [µg/ml]' 'SD[µg/ml]'};

Schonfied_2021_Fig2a = importCSVfile('/Users/migueltenorio/Documents/Data/CSV/PK/Schofield_2021/Schofield_2021_Fig2a.csv',...
    {'antiPDL1_10mgkg_C1_Time' 'antiPDL1_10mgkg_C1_Value'}, {'double' 'double'}, [3 Inf]);
antiPDL1_10mgkg_SD = postProcessPKTable(Schonfied_2021_Fig2a(1:16,1:2), 'roundDigits',4,...
    'VarNames', heading);

%%
Schonfied_2021_Fig2bi = importCSVfile('/Users/migueltenorio/Documents/Data/CSV/PK/Schofield_2021/Schofield_2021_Fig2bi.csv',...
    {'antiCTLA4_10mgkg_C1_Time' 'antiCTLA4_10mgkg_C1_Value' 'antiCTLA4_1mgkg_C1_Time'...
    'antiCTLA4_1mgkg_C1_Value'}, {'double' 'double' 'double' 'double'}, [3 Inf]);
antiCTLA4_10mgkg_MD = postProcessPKTable(Schonfied_2021_Fig2bi(1:38,1:2), 'roundDigits', 4,...
    'VarNames', heading);
antiCTLA4_1mgkg_MD = postProcessPKTable(Schonfied_2021_Fig2bi(1:28,3:4),'roundDigits', 4,...
    'VarNames', heading);

%%
Schonfied_2021_Fig2bii = importCSVfile('/Users/migueltenorio/Documents/Data/CSV/PK/Schofield_2021/Schofield_2021_Fig2bii.csv',...
    {'antiPDL1_10mgkg_C1_Time' 'antiPDL1_10mgkg_C1_Value' 'antiPDL1_1mgkg_C1_Time'...
    'antiPDL1_1mgkg_C1_Value'}, {'double' 'double' 'double' 'double'}, [3 Inf]);
antiPDL1_10mgkg_MD = postProcessPKTable(Schonfied_2021_Fig2bii(1:32,1:2), 'roundDigits',4,...
    'VarNames', heading);
antiPDL1_1mgkg_MD = postProcessPKTable(Schonfied_2021_Fig2bii(1:38,3:4),'roundDigits', 4,...
    'VarNames', heading);

%% Generate Fig 2a
errorbar(antiPDL1_10mgkg_SD{:,1}, antiPDL1_10mgkg_SD{:,2}, antiPDL1_10mgkg_SD{:,3})
set(gca, 'YScale', 'log')
ylim([1e0, 1e4])
xlim([-25, 400])
ylabel('Concentration [µg/ml]')
xlabel('Time [hr]')
grid on
title('anti-mouse PD-L1 clone 10F.9G2')

%% Generate Fig 2b
figure
subplot(2,1,1)
hold on
errorbar(antiCTLA4_10mgkg_MD{:,1}, antiCTLA4_10mgkg_MD{:,2}, antiCTLA4_10mgkg_MD{:,3})
errorbar(antiCTLA4_1mgkg_MD{:,1}, antiCTLA4_1mgkg_MD{:,2}, antiCTLA4_1mgkg_MD{:,3})
set(gca, 'YScale', 'log')
ylim([1e-3, 1e4])
xlim([-25, 525])
grid on
title('anti-mouse CTLA-4 clone 9D9')

subplot(2,1,2)
hold on
errorbar(antiPDL1_10mgkg_MD{:,1}, antiPDL1_10mgkg_MD{:,2}, antiPDL1_10mgkg_MD{:,3})
errorbar(antiPDL1_1mgkg_MD{:,1}, antiPDL1_1mgkg_MD{:,2}, antiPDL1_1mgkg_MD{:,3})
set(gca, 'YScale','log')
ylim([1e-3, 1e4])
xlim([-25, 525])
ylabel('Concentration [µg/ml]')
xlabel('Time [hr]')
grid on
title('anti-mouse PD-L1 clone 10F.9G2')


%% Add dose-normalizing columns
antiPDL1_10mgkg_SD = getIDperMl(antiPDL1_10mgkg_SD, 200);
antiPDL1_10mgkg_MD = getIDperMl(antiPDL1_10mgkg_MD, 200);
antiPDL1_1mgkg_MD = getIDperMl(antiPDL1_1mgkg_MD, 20);
antiCTLA4_10mgkg_MD = getIDperMl(antiCTLA4_10mgkg_MD, 200);
antiCTLA4_1mgkg_MD = getIDperMl(antiCTLA4_1mgkg_MD, 20);
%% Generate PI.data
expCell = {antiPDL1_10mgkg_SD; antiPDL1_10mgkg_MD; antiPDL1_1mgkg_MD;...
    antiCTLA4_10mgkg_MD; antiCTLA4_10mgkg_MD};
groupNames = {'10_mgkg_SD'; '10_mgkg_MD'; '1_mgkg_MD'; '10_mgkg_MD'; '1_mgkg_MD'};
varNames = {'C1_antiPDL1_10mgkg_SD' 'C1_antiPDL1_10mgkg_MD' 'C1_antiPDL1_1mgkg_MD'...
    'C1_antiCTLA4_10mgkg_MD' 'C1_antiCTLA4_1mgkg_MD'};
stateVar = {'Blood.antiPDL1' 'antiPDL1_ID_g_Blood_free' 'Blood.antiCTLA4' 'antiCTLA4_ID_g_Blood_free'};
PI = generatePI(expCell, varNames, stateVar, 'groupNames', groupNames, 'varIndx',...
    logical([1 1 0 0; 1 1 0 0; 1 1 0 0; 0 0 1 1; 0 0 1 1]));

%% Generate PI doses
MDSchedule = [0 0; 71 71; 167 167; 239 239; 335 335; 406 406;];
doseSchedule = {[0 0]; MDSchedule; MDSchedule; MDSchedule; MDSchedule};
doses = [10 0; 10 0; 1 0; 0 10; 0 1];
[PI, PI.u] = getDoses(PI,'doses', doses, 'Schedule', doseSchedule);
PI.variableUnits = {'µg/ml' '%ID/ml' 'µg/ml' '%ID/ml'};

%% Save
save(strjoin({cd 'PI_Schofield.mat'},'/'),'PI')