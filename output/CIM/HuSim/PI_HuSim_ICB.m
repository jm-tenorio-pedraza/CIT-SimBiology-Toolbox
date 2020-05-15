%% PK general setup

% Search paths
clear all
warning off
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/PK_mAb_ThreeComp/HuSim')
sensitivity = false;

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/CIM_4.sbproj');
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-15);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-12);
set(cs, 'MaximumWallClock', 0.25)
MODEL = 'CIM 4';
variants = get(model, 'variants');

%% Setting up parameters, data and simulations

parameters = {'Blood'; 'Tumor'; 'Peripheral'; 'CL_antiPDL1'; 'Q23'; 'Q12';...
    'kdeg_PDL1';'kin_CD8';'K_IFNg';'KDE_MDSC';'K_MDSC'; ...
    'kin_MDSC';'kin_TIC';'kpro_Tumor'; 'kpro_Tumor_Linear';'kill_CD8'; 'S_L';...
    'S_R'; 'K_CTLA4';'K_PDL1'; 'kill_Treg'};
parameters = [parameters; 'T_0'];
% Define outputs
groups_subset = {'MOC1_Control', 'MOC1_antiPDL1', 'MOC1_antiCTLA4', 'MOC1_antiCTLA4_antiPDL1'};
observables={'TV' 'CD8' 'Tumor.antiPDL1' 'Tumor.antiCTLA4'};
stateVar={'Tumor'};
doses = {'Blood.Dose_antiPDL1' 'Blood.Dose_antiCTLA4'};

PI=getPIData2('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat',...
    stateVar,groups_subset,observables,'zeroAction', 'omit','mergePhenotypes',...
    false,'output', 'mean','maxIIV', true);
PI.data = PI.data([1 4 2 3]);
PI.variableUnits={'Volume [mL]'};
PI.normIndx = [];
PI.model = 'Human simulations';
PI.observablesPlot={'TV' 'CD8' 'Tumor_antiPDL1' 'Tumor_antiCTLA4'};
PI.tspan = 1:3:365;

% Get simulation function
[sim,PI.u]=initializePI(model,parameters,observables,PI,doses, 'Human','doseUnits', 'mole');
% Get initial values
%% Doses (Durvalumab + Tremelimumab)
n_doses_antiPDL1_mono = 27;
n_doses_antiCTLA4_mono = 9;
dosing_times_antiPDL1_mono = 0:14:365;
dosing_times_antiCTLA4_mono = [0:28:28*6 (28*6+12*7):12*7:365];
dosing_times_antiPDL1_combi = [0:28:28*3 (28*3+7*2):14:365];
dosing_times_antiCTLA4_combi = 0:28:28*3;

maxDosingFreq = dosing_times_antiPDL1_mono; % Maximum dosing frequency
nDoses_max = length(maxDosingFreq);

doseScalingFactor = 77*1e-3/1.5e5*1e6; % kg*mg/kg*g/mg*mole/g*µmol/mol

control = table(maxDosingFreq', repelem(0,nDoses_max,1),...
    repelem(0/60, nDoses_max,1),'VariableNames',{'Time' 'Amount' 'Rate'});

antiPDL1_mono = table(dosing_times_antiPDL1_mono', repelem(10*doseScalingFactor,n_doses_antiPDL1_mono,1),...
    repelem(10*doseScalingFactor/60, n_doses_antiPDL1_mono,1),'VariableNames',{'Time' 'Amount' 'Rate'});

antiCTLA4_mono = table(dosing_times_antiCTLA4_mono',...
    repelem(10*doseScalingFactor,n_doses_antiCTLA4_mono,1),...
    repelem(10*doseScalingFactor/60, n_doses_antiCTLA4_mono,1),...
    'VariableNames',{'Time' 'Amount' 'Rate'});

antiPDL1_combi = table(dosing_times_antiPDL1_combi',...
    [repelem(20*doseScalingFactor, 4,1); repelem(10*doseScalingFactor, length(dosing_times_antiPDL1_combi)-4,1)],...
    [repelem(20*doseScalingFactor/60, 4,1); repelem(10*doseScalingFactor/60, length(dosing_times_antiPDL1_combi)-4,1)],...
    'VariableNames',{'Time' 'Amount' 'Rate'});

antiCTLA4_combi = table(dosing_times_antiCTLA4_combi',...
    repelem(10*doseScalingFactor,4,1), repelem(10*doseScalingFactor/60,4,1),...
    'VariableNames',{'Time' 'Amount' 'Rate'});

control.Properties.VariableUnits = {'day' 'micromole' 'micromole/minute'};
antiPDL1_mono.Properties.VariableUnits = {'day' 'micromole' 'micromole/minute'};
antiCTLA4_mono.Properties.VariableUnits = {'day' 'micromole' 'micromole/minute'};
antiPDL1_combi.Properties.VariableUnits = {'day' 'micromole' 'micromole/minute'};
antiCTLA4_combi.Properties.VariableUnits = {'day' 'micromole' 'micromole/minute'};

PI.u = {};
PI.u(1,1:2) = {control control};
PI.u(2,1:2) = {antiPDL1_mono control};
PI.u(3, 1:2) = {control antiCTLA4_mono};
PI.u(4, 1:2) = {antiPDL1_combi antiCTLA4_combi};
%% Hierarchical model simulation
PI.H = getHierarchicalStruct(parameters(1:end-1),PI,'n_sigma', length(observables),...
    'rand_indx', [14:16],'cell_indx',[2 4 12 13], 'n_indiv', length(PI.u),'CellField', 'Name');
% Generating PI
SigmaNames = getVarNames(PI, observables);
[beta, sigma_prior] = getVarValues([0.1 0.1 .001], [0.1, 0.1 0.001], [1 1 1], PI);
lb = [1e-3  1e-3    1e-4    1e-4    1e-4    1e0     1e0    1e-3   1e-3    1e-4 1e-3   1e-4    1e-3    1e-2   1e-2     1e-3 1e-7     1e-7   1e0   1e3    1e-4];
ub = [1e1   1e1     1e2     1e2     1e1     1e6    1e6   1e2     1e2     1e1 1e1    1e2     1e2     1e1    1e2  1e3 1        1      1e4   1e6     1e3];
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
%% Objective function

% Obj function
obj_fun=@(x)(likelihood_fun(x)*(-1)+prior_fun_MCMC(x)*(-1));
tic
obj_fun((finalValues))
toc


%% Posterior samples
PI_Preclinical=load('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/PK_mAb_ThreeComp/PI/ThreeComp_4/PI_PK_ThreeComp4_5_TMDD_0.mat');
postSamples_PopParams = (PI_Preclinical.PI.postSamples(:,[1:8])) ;
omega = (PI_Preclinical.PI.postSamples(:,12));
sigma = PI_Preclinical.PI.postSamples(:,14);
N_pop = 500;
N_indiv = 2;
randIndx = randsample(size(postSamples_PopParams,1),N_pop, true);

Z_indiv = (randn(N_indiv*N_pop, 2)).*repmat(omega(randIndx), N_indiv,2);
% Z = Z_indiv*PI.H.CellIndx';
scalingExp = [0.9 0.9 0.9 0.9 0.9 0.9 0 0];
scalingFactor = (77/.022).^(scalingExp);
Theta1 = [repmat(postSamples_PopParams(randIndx,:), N_indiv, 1) Z_indiv...
    repmat(omega(randIndx), N_indiv,1) repmat(sigma(randIndx), N_indiv,1) ];

Theta1(:,1:8) = Theta1(:,1:8) + log(scalingFactor);
%  plotBivariateMarginals_2(exp(Theta1(:,1:7)),...
%        'names',parameters([1:7]),'interpreter', 'tex')
% plotBivariateMarginals_2(Theta(:,8:18),...
%       'names',{PI.par([8:18]).name},'interpreter', 'tex')
Delta = repmat(postSamples_PopParams(randIndx,1:4), N_indiv, 1)- mean(repmat(postSamples_PopParams(randIndx,1:4), N_indiv, 1));
Theta2 = [repmat(log([PI.par(1:4).startValue]), N_pop*N_indiv, 1)+Delta Theta1(:,5:end)];
plotBivariateMarginals_2((Theta2(:,1:7)),...
       'names',parameters([1:7]),'interpreter', 'tex')
%% Posterior predictions
simFun=@(x)getOutput(PI,@(p)sim(p,PI.tspan(end),PI.u,PI.tspan),x,...
    @(p)getPhi2(p,PI.H,length(PI.u),'initialValue',PI.x_0),PI.normIndx, PI.H);
tic
PI=getPosteriorPredictions(exp(Theta1),PI,simFun,PI.observablesPlot);
toc
tic
PI2 = getPosteriorPredictions(exp(Theta2), PI,simFun, PI.observablesPlot);
toc
PI=getCredibleIntervals(PI,PI.observablesPlot, (Theta1),PI.H);
plotPosteriorPredictions(PI,PI.observablesPlot,'output','indiv')

PI2 = getCredibleIntervals(PI2, PI2.observablesPlot, Theta2, PI.H);
PI2_subset = PI2;
PI2_subset.data = PI2.data(5:12);
plotPosteriorPredictions(PI2_subset, PI.observablesPlot, 'output', 'indiv')
PI2.postSamples = Theta2;
PI.postSamples = Theta1;
%% Save
save(strjoin({cd 'PI_PK_TwoComp4_HuSim_Theta1'}, '/'), 'PI')
save(strjoin({cd 'PI_PK_TwoComp4_HuSim_Theta2'}, '/'), 'PI2')
