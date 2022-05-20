%% PK general setup
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
%% Load project 
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-11);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-9);
set(cs, 'MaximumWallClock', 2.5)
MODEL = 'CIM 10';
variants = get(model, 'variants');
set(cs, 'Time','day')
set(cs, 'StopTime',365*3)
%% Setting up parameters, data and simulations

parameters = {'Blood'; 'Peripheral';'CL_antiPDL1';'Q12';...
    'f_L'; 'V_vs';'kin_CD8'; 'kin_Treg'; 'kin_DC';'kin_MDSC';'kin_CD45';...
    'kpro_Tumor'; 'kpro_Tumor_Linear';'kill_CD8'; ...
    'K_MDSC';'K_IFNg';'K_CTLA4'; 'K_PDL1'; 'KDE_MDSC';'KDE_Treg';
    'V_intra';'PDL1_Tumor_ss';'PDL1_Immune_ss';'kel_CD45';'kel_DC';'kel_Effector';'kel_MDSC';...
    'kel_Naive';'kel_Treg';'T_0'};
% Define outputs
groups_subset = {'MOC1_Control', 'MOC1_antiPDL1', 'MOC1_antiCTLA4', 'MOC1_antiCTLA4_antiPDL1'};
observables={'Tumor' 'TumorInter.Ag' 'TumorInter.Cancer' 'TumorInter.CD8_N'...
    'TumorInter.CD8_E'  'TumorInter.CD4Foxp3' 'TumorInter.DC'...
    'TumorInter.GMDSC' 'TumorInter.PDL1_Tumor' 'TumorInter.PDL1_Immune'...
    'TumorInter.CTLA4_Treg' 'TumorInter.CTLA4_CD8'  'TumorInter.antiPDL1' 'TumorInter.antiCTLA4'};
stateVar={'Tumor'};
doses = {'Blood.antiPDL1' 'Blood.antiCTLA4'};
%%
PI=getPIData4(data_ext, stateVar,groups_subset,'output', 'mean',...
    'responseGrouping', true, 'kineticGrouping', false, 'zeroHandling',...
    'replace','logTransform',false);
PI.data = PI.data([1 6 2 4]);
cellgroup=repelem({'MOC1'},4);
[PI.data(1:end).Cell]=cellgroup{:,:};
PI.data(1).Group='MOC1_Control';
PI.variableUnits={'Volume [mL]'};
PI.normIndx = [];
PI.model = 'Human simulations';
PI.observablesPlot={'TV' 'Ag' 'CancerCells' 'CD8_N' 'CD8_E' 'CD4Foxp3' 'DC'...
    'GMDSC' 'PDL1_Tumor' 'PDL1_Immune' 'CTLA4_Treg' 'CTLA4_CD8' 'TumorInter_antiPDL1'...
    'TumorInter_antiCTLA4'};
PI.tspan = 1:14:(365*3);

% Get simulation function
[sim,PI.u]=initializePI(model,parameters,observables,PI,doses, 'Human','doseUnits', 'mole');
PI = getDoseRegimen(PI);

%% Hierarchical model simulation
PI.H = getHierarchicalStruct3(parameters(1:end-1),PI,'n_sigma', length(observables),...
    'rand_indx', [] , 'cell_indx',[], 'resp_indx', [],'n_indiv', length(PI.u),...
    'CellField','Cell');% Generating PI
SigmaNames = getVarNames(PI, observables);
[beta, sigma_prior] = getVarValues2([.1 .1 .1 .1], ...
    [.6 .6 .6 .6], [1 1 1 1], PI);
PI.par = getParamStruct3(sim,PI.H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior,'startSigma',...
    ones(length(PI.H.SigmaParams), 1)*.6, 'ref', 'ones');
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
PI_PK = load(strjoin({cd 'PK_mAb_ThreeComp/PI/CIM10/PI_CIM10_PK_1.mat'},'/'));
PI_MOC1 = load(strjoin({cd 'CIM/PI/CIM33/PI_CIM10_MOC1_PKPD_1_8_kinCD8.mat'},'/'));
PI_MOC2 = load(strjoin({cd 'CIM/PI/CIM33/PI_CIM10_MOC2_PKPD_1_8.mat'},'/'));
PI_MC38 = load(strjoin({cd 'CIM/PI/CIM33/PI_CIM10_MC38_PKPD_1_8.mat'},'/'));

Meta=[];
Meta(1).Struct = PI_PK;
Meta(1).Struct = PI_MC38;
Meta(2).Struct = PI_MOC1;
Meta(3).Struct = PI_MOC2;
values=repelem({sim.Parameters.Value},length(Meta),1);
[Meta(1:end).Values]=values{:,:};
%% Select which dimensions to sample from 
N_pop =5e2;
[theta1,~] = getTheta2(Meta, parameters, N_pop, 'covStructure', 'independent','fixedParams','fixed');
[theta2,~] = getTheta2(Meta, parameters, N_pop, 'covStructure', 'independent');
% Indexes of parameter uncertainties to propagate
PK_params={'Blood' 'Peripheral' 'CL_antiPDL1'...
    'Q12' 'f_L' 'V_vs'};
PD_params = {'kpro_Tumor'};
indx_pk = find(ismember(parameters,PK_params));
indx_pd = find(ismember(parameters, PD_params));
priorMu_pk = log([PI.par(indx_pk).startValue]);
priorMu_pd = log([PI.par(indx_pd).startValue]);

Theta = theta1;
% PK parameters
[pk_indx, pk_order]=ismember({PI_PK.PI.par(:).name},PK_params);
pk_order=pk_order(~pk_order==0);
theta_pk = PI_PK.PI.postSamples(:,ismember({PI_PK.PI.par(:).name},PK_params));
sampleIndx_pk=randsample(size(theta_pk,1),N_pop*3, true);
delta_pk = theta_pk(sampleIndx_pk,:)-mean(theta_pk);
Theta(:,indx_pk) =priorMu_pk+delta_pk(:,pk_order);

% PD parameters
[pd_indx, pd_order]=ismember(parameters,PD_params);
pd_order=pd_order(~pd_order==0);
theta_pd = Theta(:,pd_indx);
delta_pd = theta_pd-mean(theta_pd);
Theta(:,indx_pd) =priorMu_pd + delta_pd(:,pd_order);

TV_range = [pi*(4/3)*(.5^3) pi*(4/3)*(5^3)]; % SLD_0: 0.5 to 5 cm
T_0 = log((rand(size(Theta,1),1)*(TV_range(2)-TV_range(1))... % Initial T_0 (1e6 cells)
    +TV_range(1))/1.53e-3.*exp(Theta(:,ismember(parameters,'V_intra'))));
Theta(:,ismember(parameters,'T_0')) = T_0;
%% Set linear tumor growth to 0
kproTumorLinear_indx=ismember(parameters,'kpro_Tumor_Linear');
Theta(:,kproTumorLinear_indx)=-log(1000);

%% Scaling parameters
a=1;    b=.9;  c=.8;   d=.75;  e=.65; f=.55; g=-.25;
%%               1 2 3    4   5   6     7   8       9 101112    13141516           17   18  19202122    23  24  25  26  27  28  29  30  
scalingExp1 =   [0 0 0    0   0   0     0   0       0 0 0 0     0 0 0 0            0    0   0 0 0 0     0   0   0   0   0   0   0   0]; % none
scalingExp2 =   [0 0 0    0   0   0     d   d       d d d 0     0 0 0 0            0    0   0 0 0 0     0   0   0   0   0   0   0   0]; % Cell infiltration at .75
scalingExp3 =   [0 0 0    0   0   0     b   b       b b b 0     0 0 0 0            0    0   0 0 0 0     0   0   0   0   0   0   0   0]; % Cell infiltration at .9
scalingExp4 =   [0 0 0    0   0   0     a   a       a a a 0     0 0 0 0            0    0   0 0 0 0     0   0   0   0   0   0   0   0]; % Cell infiltration at 1
scalingExp5 =   [0 0 0    0   0   0     d   d       d d d 0     0 0 d 0            0    0   0 0 0 0     0   0   0   0   0   0   0   0]; % Cell infiltration at .75 and killing rate at .75
scalingExp6 =   [0 0 0    0   0   0     b   b       b b b 0     0 0 b 0            0    0   0 0 0 0     0   0   0   0   0   0   0   0]; % Cell infiltration at .9 and killing rate at .9
scalingExp7 =   [0 0 0    0   0   0     a   a       a a a 0     0 0 a 0            0    0   0 0 0 0     0   0   0   0   0   0   0   0]; % Cell infiltration at 1 and killing rate at 1
scalingExp8 =   [0 0 0    0   0   0     d   d       d d d 0     0 0 a a            0    0   a a 0 0     0   0   0   0   0   0   0   0]; % Cell infiltration at .75 and cell amounts at 1 
scalingExp9 =   [0 0 0    0   0   0     b   b       b b b 0     0 0 a a            0    0   a a 0 0     0   0   0   0   0   0   0   0]; % Cell infiltration at .9 and cell amounts at 1
scalingExp10=   [0 0 0    0   0   0     a   a       a a a 0     0 0 a a            0    0   a a 0 0     0   0   0   0   0   0   0   0]; % Cell infiltration at 1 and cell amounts at 1
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

Theta1 = Theta + log(scalingFactor1);
Theta2 = Theta + log(scalingFactor2);
Theta3 = Theta + log(scalingFactor3);
Theta4 = Theta + log(scalingFactor4);
Theta5 = Theta + log(scalingFactor5);
Theta6 = Theta + log(scalingFactor6);
Theta7 = Theta + log(scalingFactor7);
Theta8 = Theta + log(scalingFactor8);
Theta9 = Theta + log(scalingFactor9);
Theta10= Theta + log(scalingFactor10);

HumanParams=table( exp(mean((Theta1))'),exp(mean((Theta2))'),exp(mean((Theta3))'),...
    exp(mean((Theta4))'),exp(mean((Theta5))'),exp(mean((Theta6))'),...
    exp(mean((Theta7))'), exp(mean((theta1)))',sim.Parameters.Units,'VariableNames', ...
    {'Human1' 'Human2' 'Human3' 'Human4' 'Human5' 'Human6' 'Human7'  'Mouse','Units'}, 'RowNames', parameters)
clearvars beta delta doses finalValues groups_subset lb MODEL PI_CIM PI_ICB PI_PK...
    sigma_indx sigmga_prior ub variants scalingExp1 scalingExp2 scalingExp3 ...
    scalingExp4 scalingExp5 scalingExp6 scalingExp7 scalingExp8 scalingExp9 ...
    scalingExp10 scalingExp11 
%% Plot bivariate marginals
plotBivariateMarginals_2(exp(Theta3(:,:)),...
       'names',parameters,'interpreter', 'tex')
plotCorrMat(Theta2, parameters)
% Plot all samples
plotBivariateMarginals_2((theta1(:,:)),...
       'names',parameters,'interpreter', 'tex')
plotCorrMat(theta1, parameters)
%Plot 1st cell line samples
plotBivariateMarginals_2((theta1(1:N_pop,:)),...
       'names',parameters,'interpreter', 'tex')

%% Posterior predictions
simFun=@(x)getOutput2(PI,@(p)sim(p,PI.tspan(end),PI.u,PI.tspan),x,...
    {},PI.normIndx,'Output', 'data');
dataOutput = simFun(exp(mean(Theta1(:,:))));
[PI.data(1:end).simValue] = dataOutput{:,:};
figure
hold on
arrayfun(@(x)plot(PI.tspan, x.simValue(:,1)), PI.data, 'UniformOutput', false)
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
%% plot Each species
figure
nrow=ceil(sqrt(length(PI.observablesPlot)));
ncol=ceil(length(PI.observablesPlot)/nrow);
colors=linspecer(length(PI.data));
for i=1:length(PI.observablesPlot)
    output_i=PI.observablesPlot{i};
    subplot(nrow,ncol,i)
    hold on
for j=1:length(PI.data)
    plot(PI4.tspan,mean(PI4.output(j).(output_i),1,'omitnan'),'Color',colors(j,:));
end
title(output_i)
end
legend([PI.data(:).Group])
%% Plot posterior predictions
treatments = {'Control' 'antiPDL1' 'antiCTLA4' 'antiPDL1_antiCTLA4'};

% Calculate SLD and response %
PI1 = getSLD(PI1, sigma, 'TV_0', exp(Theta1(:,end))*.00153/.375);
[PI1, ORR1] = getORR(PI1, treatments,'TV_0', exp(Theta1(:,end))*.00153/.375);
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
[PI1,~, ~] = getSurvivalTime(PI1, treatments,'N',size(Theta1,1));

[PI2, response2] = getPFS(PI2, treatments);
[PI2,T2, censor2] = getSurvivalTime(PI2, treatments);

[PI3, response3] = getPFS(PI3, treatments);
[PI3,~, ~] = getSurvivalTime(PI3, treatments);

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
plotSurvivalFunction(PI1,PI.tspan(end),treatments)
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
