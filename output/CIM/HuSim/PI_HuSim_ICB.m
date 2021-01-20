%% PK general setup
% Search paths
clear all
warning off
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output')

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/sbio-projects/CIM_4.sbproj');
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(cs, 'MaximumWallClock', 2.5)
MODEL = 'CIM 4';
variants = get(model, 'variants');
sbioaccelerate(model, cs)
 %% Setting up parameters, data and simulations

parameters = {'Blood'; 'Tumor'; 'Peripheral'; 'CL_antiPDL1'; 'Q23'; 'Q12';...
    'PDL1_Tumor' ;'kdeg_PDL1';'kin_CD8';'K_IFNg';'KDE_MDSC';'K_MDSC'; ...
    'kin_MDSC';'kin_TIC';'kpro_Tumor'; 'kpro_Tumor_Linear';'kill_CD8'; 'S_L';...
    'S_R'; 'K_CTLA4';'K_PDL1'; 'kill_Treg'};
parameters = [parameters; 'T_0'];                
% Define outputs
groups_subset = {'MOC1_Control', 'MOC1_antiPDL1', 'MOC1_antiCTLA4', 'MOC1_antiCTLA4_antiPDL1'};
observables={'TV' 'CD8' 'Tumor.antiPDL1' 'Tumor.antiCTLA4' 'IAR' 'kcoi' 'kco'};
stateVar={'Tumor'};
doses = {'Blood.Dose_antiPDL1' 'Blood.Dose_antiCTLA4'};

PI=getPIData2('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat',...
    stateVar,groups_subset,observables,'zeroAction', 'omit','mergePhenotypes',...
    false,'output', 'mean','maxIIV', true);
PI.data = PI.data([1 4 2 3]);
PI.variableUnits={'Volume [mL]'};
PI.normIndx = [];
PI.model = 'Human simulations';
PI.observablesPlot={'TV' 'CD8' 'Tumor_antiPDL1' 'Tumor_antiCTLA4' 'IAR' 'kcoi' 'kco'};
PI.tspan = 1:3:365;

% Get simulation function
[sim,PI.u]=initializePI(model,parameters,observables,PI,doses, 'Human','doseUnits', 'mole');
PI = getDoseRegimen(PI);

%% Hierarchical model simulation
PI.H = getHierarchicalStruct(parameters(1:end-1),PI,'n_sigma', length(observables),...
    'rand_indx', [14:16],'cell_indx',[], 'n_indiv', length(PI.u),'CellField', 'Name');
% Generating PI
SigmaNames = getVarNames(PI, observables);
[beta, sigma_prior] = getVarValues([0.1 0.1 .001], [0.1, 0.1 0.001], [1 1 1], PI);
lb = [1e-3  1e-3    1e-4    1e-4    1e-4    1e0    1e0  1e0    1e-3   1e-3    1e-4 1e-3   1e-4    1e-3    1e-2   1e-2     1e-3 1e-7     1e-7   1e0   1e3    1e-4];
ub = [1e1   1e1     1e2     1e2     1e1     1e6    2e6  1e6   1e2     1e2     1e1 1e1    1e2     1e2     1e1    1e2  1e3 1        1      1e4   1e6     1e3];
PI.par = getParamStruct2(sim,PI.H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior,'LB', lb', 'UB', ub');
try
    finalValues =log([PI.par(:).finalValue]);
catch
    finalValues =log([PI.par(:).startValue]);

end

PI= assignPrior(PI);
% Log-ikelihood function
likelihood_fun=@(p)likelihood(exp(p),sim,PI,'censoring',false);
prior_fun_MCMC=@(p)getPriorPDFMCMC2(exp(p),PI);
paramNames = getParamNames(PI,sim, observables);
%% Load PIs and posterior parameter samples
PI_PK = load(strjoin({cd 'PK_mAb_ThreeComp/PI/ThreeComp_4/PI_PK_ThreeComp4_4_TMDD_11.mat'},'/'));
PI_CIM = load(strjoin({cd 'CIM/PI/CIM21/PI_CIM21_Control_14.mat'},'/'));
PI_ICB = load(strjoin({cd 'CIM/PI/CIM21/PI_CIM4_ICB_21_14.mat'},'/'));
Meta=[];
Meta(1).Struct = PI_PK;
Meta(2).Struct = PI_CIM;
Meta(3).Struct = PI_ICB;
%% Select which dimensions to sample from 
N_pop = 500;
N_indiv = 5;
N_cell = 2;
[theta,H] = getTheta(Meta, parameters, N_pop, N_cell, N_indiv, 'covStructure', 'independent');

% Indexes of parameter uncertainties to propagate
indx_un = find(ismember(parameters, {'Blood' 'Tumor' 'Peripheral' 'CL_antiPDL1'...
    'Q23' 'Q12' 'kpro_Tumor'}));
priorMu = log([PI.par(indx_un).startValue]);

Theta = theta;
delta = theta(:,indx_un)-mean(theta(:,indx_un));
Theta(:,indx_un) =priorMu+delta;

TV_range = [pi*(4/3)*(.5^3) pi*(4/3)*(5^3)]; % SLD_0: 0.5 to 5 cm
T_0 = log((rand(size(theta,1),1)*(TV_range(2)-TV_range(1))... % Initial T_0 (1e3 cells/µliter)
    +TV_range(1))/1.53e-9*1e-6);
Theta(:,end) = T_0;
%% Scaling parameters
%%               1 2 3    4   5   6     7   8       9 101112    13141516           17   18  19202122    23  
scalingExp1 =   [0 0 0    0   0   0     0   0       0 0 0 0     0 0 0 0            0    0   0 0 0 0     0]; % none
scalingExp2 =   [0 0 0    0   0   0     0   0       0 0 0 1     0 0 0 1            0    1   1 0 0 0     0]; % Volume parameters
scalingExp3 =   [0 0 0    0   0   0     0   -.25    0 0 0 1     0 0 0 1         -.25    1   1 0 0 -.25  0]; % Volume, killing rate, PDL1 turnover rate
scalingExp4 =   [0 0 0    0   0   0     0   0       0 0 0 1     0 0 0 0            0    1   1 0 0 0     0]; % Volume params w/o kpro_Tumor_Linear
scalingExp5 =   [0 0 0    0   0   0     0   -.25    0 0 0 1     0 0 0 0         -.25    1   1 0 0 -.25  0]; % Volume params and killing rate w/o kpro_Tumor_Linear
scalingExp6 =   [0 0 0    0   0   0     0   0       0 0 0 1     0 0 0 .75          0    1   1 0 0 0     0]; % Volume and killing rate params

scalingExp8 =   [0 0 0    0   0   0     0   0       0 0 0 .9    0 0 0 0            0   .9  .9 0 0 0     0]; % Volume params w/o kpro_Tumor_Linear
scalingExp9 =   [0 0 0    0   0   0     0   -.25    0 0 0 .9    0 0 0 0         -.25   .9  .9 0 0 -.25  0]; % Volume params and killing rate w/o kpro_Tumor_Linear
scalingExp10 =  [0 0 0    0   0   0     0   -.25    0 0 0 .9    0 0 0 .65       -.25   .9  .9 0 0 -.25  0]; % Volume params and killing rate w/o kpro_Tumor_Linear
scalingExp11 =  [0 0 0    0   0   0     0   -.25    0 0 0 .8    0 0 0 .55       -.25   .8  .8 0 0 -.25  0]; % Volume params and killing rate w/o kpro_Tumor_Linear
scalingExp12 =  [0 0 0    0   0   0     0   -.25    0 0 0 .8    0 0 0 .55       0      .8  .8 0 0 0     0]; % Volume params and killing rate w/o kpro_Tumor_Linear
scalingExp13 =  [0 0 0    0   0   0     0   -.25    0 0 0 1     0 0 0 1         0      1   1  0 0 0     0]; % Volume params and killing rate w/o kpro_Tumor_Linear

scalingFactor1 = (77/.022).^(scalingExp1);
scalingFactor2 = (77/.022).^(scalingExp2);
scalingFactor3 = (77/.022).^(scalingExp3);
scalingFactor4 = (77/.022).^(scalingExp4);
scalingFactor5 = (77/.022).^(scalingExp5);
scalingFactor6 = (77/.022).^(scalingExp6);
scalingFactor7 = (77/.022).^(scalingExp7);
scalingFactor8 = (77/.022).^(scalingExp8);
scalingFactor9 = (77/.022).^(scalingExp9);
scalingFactor10 = (77/.022).^(scalingExp10);
scalingFactor11 = (77/.022).^(scalingExp11);
scalingFactor12 = (77/.022).^(scalingExp12);
scalingFactor13 = (77/.022).^(scalingExp13);

Theta1 = Theta + log(scalingFactor1);
Theta2 = Theta + log(scalingFactor2);
Theta3 = Theta + log(scalingFactor3);
Theta4 = Theta + log(scalingFactor4);
Theta5 = Theta + log(scalingFactor5);
Theta6 = Theta + log(scalingFactor6);
Theta7 = [theta + log(scalingFactor7)];
Theta7(:,[indx_un(end) end]) = [priorMu(end) + delta(:,(end)) T_0];
Theta8 = Theta + log(scalingFactor8);
Theta9 = Theta + log(scalingFactor9);
Theta10 = Theta + log(scalingFactor10);
Theta11 = Theta + log(scalingFactor11);
Theta12 = Theta + log(scalingFactor12);
Theta13 = Theta + log(scalingFactor13);

table( exp(mean((Theta1))'),exp(mean((Theta2))'),exp(mean((Theta3))'),...
    exp(mean((Theta4))'),exp(mean((Theta5))'),exp(mean((Theta6))'),...
    exp(mean((Theta7))'),exp(mean((Theta8))'),exp(mean((Theta9))'),...
    exp(mean((Theta10))'),exp(mean((Theta11))'),exp(mean((Theta12))'),...
    exp(mean((Theta13))'), exp(mean((theta))'),'VariableNames', ...
    {'Human1' 'Human2' 'Human3' 'Human4' 'Human5' 'Human6' 'Human7' 'Human8'...
    'Human9' 'Human10' 'Human11' 'Human12' 'Human13' 'Mouse'}, 'RowNames', parameters)
sigma = nan(N_pop*N_indiv*N_cell,2);

sigma_indx = randsample(size(PI_CIM.PI.postSamples,1),N_pop*N_indiv*N_cell,true);
sigma(:,2) = PI_CIM.PI.postSamples(sigma_indx,end-6);
sigma_indx = randsample(size(PI_ICB.PI.postSamples,1),N_pop*N_indiv*N_cell,true);
sigma(:,1) = PI_ICB.PI.postSamples(sigma_indx,end);
sigma = exp(sigma);
clearvars beta delta doses finalValues groups_subset lb MODEL PI_CIM PI_ICB PI_PK...
    sigma_indx sigmga_prior ub variants scalingExp1 scalingExp2 scalingExp3 ...
    scalingExp4 scalingExp5 scalingExp6 scalingExp7 scalingExp8 scalingExp9 ...
    scalingExp10 scalingExp11 
%% Plot bivariate marginals
plotBivariateMarginals_2(exp(Theta3(:,:)),...
       'names',parameters,'interpreter', 'tex')
plotCorrMat(Theta2, parameters)

%% Posterior predictions
simFun=@(x)getOutput2(PI,@(p)sim(p,PI.tspan(end),PI.u,PI.tspan),x,...
    {},PI.normIndx,'Output', 'data');
dataOutput = simFun(exp(mean(Theta3(:,:))));
[PI.data(1:end).simValue] = dataOutput{:,:};
figure
hold on
arrayfun(@(x)plot(PI.tspan, x.simValue(:,1), 'Color', x.colors), PI.data, 'UniformOutput', false)
title('IAR')

tic
PI1=getPosteriorPredictions2(exp(Theta1(:,:)),PI,simFun,PI.observablesPlot);
toc

tic
PI2=getPosteriorPredictions2(exp(Theta2(1:end,:)),PI,simFun,PI.observablesPlot);
toc

tic
PI3=getPosteriorPredictions2(exp(Theta3(1:end,:)),PI,simFun,PI.observablesPlot);
toc

tic
PI4=getPosteriorPredictions2(exp(Theta4(1:end,:)),PI,simFun,PI.observablesPlot);
toc

tic
PI5=getPosteriorPredictions2(exp(Theta5(1:end,:)),PI,simFun,PI.observablesPlot);
toc

tic
PI6=getPosteriorPredictions2(exp(Theta6(1:end,:)),PI,simFun,PI.observablesPlot);
toc

tic
PI7=getPosteriorPredictions2(exp(Theta7(1:end,:)),PI,simFun,PI.observablesPlot);
toc

tic
PI8=getPosteriorPredictions2(exp(Theta8(1:end,:)),PI,simFun,PI.observablesPlot);
toc

tic
PI9=getPosteriorPredictions2(exp(Theta9(1:end,:)),PI,simFun,PI.observablesPlot);
toc

tic
PI10=getPosteriorPredictions2(exp(Theta10(1:end,:)),PI,simFun,PI.observablesPlot);
toc

tic
PI11=getPosteriorPredictions2(exp(Theta11(1:end,:)),PI,simFun,PI.observablesPlot);
toc


PI12=getPosteriorPredictions2(exp(Theta12(1:end,:)),PI,simFun,PI.observablesPlot);

PI13=getPosteriorPredictions2(exp(Theta13(1:end,:)),PI,simFun,PI.observablesPlot);

%% Plot posterior predictions
treatments = {'Control' 'antiPDL1' 'antiCTLA4' 'antiPDL1_antiCTLA4'};

% Calculate SLD and response %
PI1 = getSLD(PI1, sigma, 'TV_0', exp(Theta1(:,end))*.00153);
[PI1, ORR1] = getORR(PI1, treatments,'TV_0', exp(Theta1(:,end))*.00153);
PI2 = getSLD(PI2, sigma, 'TV_0', exp(Theta2(:,end))*.00153);
[PI2, ORR2] = getORR(PI2, treatments,'TV_0', exp(Theta2(:,end))*.00153);
PI3 = getSLD(PI3, sigma, 'TV_0', exp(Theta3(:,end))*.00153);
[PI3, ORR3] = getORR(PI3, treatments,'TV_0', exp(Theta3(:,end))*.00153);
PI4 = getSLD(PI4, sigma, 'TV_0', exp(Theta4(:,end))*.00153);
[PI4, ORR4] = getORR(PI4, treatments,'TV_0', exp(Theta4(:,end))*.00153);
PI5 = getSLD(PI5, sigma, 'TV_0', exp(Theta5(:,end))*.00153);
[PI5, ORR5] = getORR(PI5, treatments,'TV_0', exp(Theta6(:,end))*.00153);
PI6 = getSLD(PI6, sigma, 'TV_0', exp(Theta6(:,end))*.00153);
[PI6, ORR6] = getORR(PI6, treatments,'TV_0', exp(Theta6(:,end))*.00153);
PI7 = getSLD(PI7, sigma, 'TV_0', exp(Theta7(:,end))*.00153);
[PI7, ORR7] = getORR(PI7, treatments,'TV_0', exp(Theta7(:,end))*.00153);
PI8 = getSLD(PI8, sigma, 'TV_0', exp(Theta8(:,end))*.00153);
[PI8, ORR8] = getORR(PI8, treatments,'TV_0', exp(Theta8(:,end))*.00153);
PI9 = getSLD(PI9, sigma, 'TV_0', exp(Theta9(:,end))*.00153);
[PI9, ORR9] = getORR(PI9, treatments,'TV_0', exp(Theta9(:,end))*.00153);
PI10 = getSLD(PI10, sigma, 'TV_0', exp(Theta10(:,end))*.00153);
[PI10, ORR10] = getORR(PI10, treatments,'TV_0', exp(Theta10(:,end))*.00153);
PI11 = getSLD(PI11, sigma, 'TV_0', exp(Theta11(:,end))*.00153);
[PI11, ORR11] = getORR(PI11, treatments,'TV_0', exp(Theta11(:,end))*.00153);
PI12 = getSLD(PI12, sigma, 'TV_0', exp(Theta12(:,end))*.00153);
[PI12, ORR12] = getORR(PI12, treatments,'TV_0', exp(Theta12(:,end))*.00153);
PI13 = getSLD(PI13, sigma, 'TV_0', exp(Theta13(:,end))*.00153);
[PI13, ORR13] = getORR(PI13, treatments,'TV_0', exp(Theta13(:,end))*.00153);

%%
% Plot SLD
plotORR(PI1, treatments,'output', {'SLD'  })
plotORR(PI2, treatments,'output', {'SLD'  })
plotORR(PI3, treatments,'output', {'SLD' })
plotORR(PI4, treatments,'output', {'SLD' })
plotORR(PI5, treatments,'output', {'SLD' })
plotORR(PI6, treatments,'output', {'SLD' })
plotORR(PI7, treatments,'output', {'SLD' })
plotORR(PI8, treatments,'output', {'SLD' })
plotORR(PI9, treatments,'output', {'SLD' })
plotORR(PI10, treatments,'output', {'SLD' })
plotORR(PI11, treatments,'output', {'SLD' })
plotORR(PI12, treatments,'output', {'SLD' })
plotORR(PI13, treatments,'output', {'SLD' })

% Plot individual time courses
plotORR(PI2, treatments, 'output', {'SLD'},'lines', 'individual')
plotORR(PI4, treatments, 'output', {'SLD'},'lines', 'individual')
plotORR(PI5, treatments, 'output', {'SLD'},'lines', 'individual')

%% PFS
% Calculate SLD and response %
[PI1, response1] = getPFS(PI1, treatments);
[PI1,~, ~] = getSurvivalTime(PI1, treatments, Theta2);

[PI2, response2] = getPFS(PI2, treatments);
[PI2,T2, censor2] = getSurvivalTime(PI2, treatments);

[PI3, response3] = getPFS(PI3, treatments);
[PI3,~, ~] = getSurvivalTime(PI3, treatments, Theta2);

[PI4, response4] = getPFS(PI4, treatments);
[PI4,~, ~] = getSurvivalTime(PI4, treatments, Theta5);

[PI5, response5] = getPFS(PI5, treatments);
[PI5,~, ~] = getSurvivalTime(PI5, treatments, Theta5);

[PI6, response6] = getPFS(PI6, treatments);
[PI6,~, ~] = getSurvivalTime(PI6, treatments, Theta6);

[PI7, response7] = getPFS(PI7, treatments);
[PI7,~, ~] = getSurvivalTime(PI7, treatments, Theta7);

[PI8, response8] = getPFS(PI8, treatments);
[PI8,~, ~] = getSurvivalTime(PI8, treatments, Theta8);

[PI9, response9] = getPFS(PI9, treatments);
[PI9,T9, censor9] = getSurvivalTime(PI9, treatments);

[PI10, response10] = getPFS(PI10, treatments);
[PI10,~, ~] = getSurvivalTime(PI10, treatments);

[PI11, response11] = getPFS(PI11, treatments);
[PI11,T11, censor11] = getSurvivalTime(PI11, treatments);

[PI12, response12] = getPFS(PI12, treatments);
[PI12,~, ~] = getSurvivalTime(PI12, treatments, Theta9);

[PI13, response13] = getPFS(PI13, treatments);
[PI13,~, ~] = getSurvivalTime(PI13, treatments, Theta9);


%% Survival functions
plotSurvivalFunction(PI1,364,treatments)
plotSurvivalFunction(PI2,364,treatments)
plotSurvivalFunction(PI3,364,treatments)
plotSurvivalFunction(PI4,364,treatments)
plotSurvivalFunction(PI5,364,treatments)
plotSurvivalFunction(PI6,364,treatments)
plotSurvivalFunction(PI7,364,treatments)
plotSurvivalFunction(PI8,364,treatments)
plotSurvivalFunction(PI9,364,treatments)
grid on
plotSurvivalFunction(PI10,364,treatments)
plotSurvivalFunction(PI11,364,treatments)
plotSurvivalFunction(PI12,364,treatments)
plotSurvivalFunction(PI13,364,treatments)

%% Analysing parameter-output relations
corrInOut = plotInputToOutput(Theta2, {'SLD'}, PI2, treatments, parameters);
corrInOut2 = plotInputToOutput(PI1.Theta, {'SLD'}, PI1, treatments, parameters,'plotOutput', true);
figure
hold on
arrayfun(@(x)plot((x.CD8(:,end)),(x.TV(:,end)), '+'), PI2.output)
%% Comparing ORR
ORR = {ORR1 ORR2 ORR3 ORR4 ORR5 ORR6 ORR7 ORR8 ORR9 ORR10 ORR11 ORR12 ORR13};
parametrizationNames =   {'Par1' 'Par2' 'Par3' 'Par4' 'Par5' 'Par6' 'Par7' 'Par8' 'Par9' 'Par10'...
        'Par11' 'Par12' 'Par13'};
controlORR = array2table(cell2mat(cellfun(@(x) x{:,2}, ORR,...
    'UniformOutput', false)), 'RowNames', {'PD' 'SD' 'PR' 'CR'}, 'VariableNames',...
  parametrizationNames);
antiPDL1ORR = array2table(cell2mat(cellfun(@(x) x{:,3}, ORR, ...
    'UniformOutput', false)),'RowNames', {'PD' 'SD' 'PR' 'CR'}, 'VariableNames',...
  parametrizationNames );
antiCTLA4ORR = array2table(cell2mat(cellfun(@(x) x{:,4}, ORR,...
    'UniformOutput', false)),'RowNames', {'PD' 'SD' 'PR' 'CR'}, 'VariableNames',...
  parametrizationNames) ;
antiPDL1_antiCTLA4ORR = array2table(cell2mat(cellfun(@(x) x{:,5}, ORR, ...
    'UniformOutput', false)),'RowNames', {'PD' 'SD' 'PR' 'CR'}, 'VariableNames',...
  parametrizationNames);


%% 
pfs2 = getMedianPFS(PI2,T2, censor2,treatments);
pfs9 = getMedianPFS(PI9,T9, censor9,treatments);
pfs11 = getMedianPFS(PI11,T11, censor11,treatments);

%% Plot predictions of ORR against realized values
ORR_Clinical = table({'PD' 'SD' 'PR' 'CR'}', 'VariableNames', {'Response'});
ORR_Clinical(:, end+1) = {nan(1,1)};
ORR_Clinical(:, end+1) = {64.6  6.2     9.2 0}';
ORR_Clinical(:, end+1) = {69.8  0       1.6 0}';
ORR_Clinical(:, end+1) = {64.3  5.4     7.8 0}';
ORR_Clinical.Properties.VariableNames(2:end) = treatments;
ORR_Preclinical = ORR9;
figure
hold on
markers = {'d' 'o' 's' '*'};
for i=1:length(treatments)
    indx = (i-1)*2+2;
    CI = cellfun(@(x)str2double(x), ORR_Preclinical{:, indx+1});

    for j=1:4
    h = errorbar(ORR_Clinical{j,i+1}, ORR_Preclinical{j, indx}*100,...
        (CI(j,1)-ORR_Preclinical{j, indx})*100,(CI(j,2)-ORR_Preclinical{j, indx})*100);
    try
    h.Color = PI1.data(i).colors;
    h.LineStyle = 'none';
    h.Marker = markers{j};
    h.MarkerEdgeColor = PI1.data(i).colors;
    h.MarkerFaceColor = PI1.data(i).colors;
    h.MarkerSize = 12;
    catch
    end
    end
end
ax = gca;
plot(0:10:100, 0:10:100, '-k')
grid on
legend(ax.Children(3:4:end),treatments(end:-1:1),'interpreter', 'none')
title('RECIST 1.1. Classification comparison between simulations and realized clinical values')
xlabel('Clinical estimates [% of patients]')
ylabel('Clinical simulation [% of simulations]')


%% MEdian PFS
PFS_Clinical = table(treatments');
PFS_Clinical(:,end+1) = {nan 1.9 2 1.9}';
PFS_Clinical(:,end+1) = {nan 1.8 1.8 1.9}';
PFS_Clinical(:,end+1) = {nan 2.8 2 2.1}';
figure
hold on
for i=1:4

h=errorbar(PFS_Clinical{i,2}, pfs9{i,2}, PFS_Clinical{i,2}-PFS_Clinical{i,3},...
    PFS_Clinical{i,2}-PFS_Clinical{i,4}, pfs9{i,2}-pfs9{i,3}, pfs9{i,2}-pfs9{i,4});
h.MarkerFaceColor = [PI.data(i).colors];

end

%% Save output
for i=1:13
    num_i = num2str(i);
    save(strjoin({cd '/CIM/HuSim/CIM21/PI' num_i '.mat'},''),strjoin({'PI', num_i}, ''))
end
