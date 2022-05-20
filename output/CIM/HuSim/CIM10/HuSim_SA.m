%% General setup for Kosinsky et al with k_pro as the parameter to vary for each individual
%% Search paths
clear all
warning off
opengl('save', 'software')
%%
if ispc
    addpath(genpath('\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox'))
    cd('\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\output')
    out = sbioloadproject('\Users\jmten\OneDrive\Dokumente\GitHub\sbio-projects\CIM_10_1_PKPD.sbproj');
    data_ext = {'\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\PI_Clavijo_2.mat'};
    data_ext1 = {'\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\PI_Morisada_3.mat'};
else
    addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
    cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output')
    out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/sbio-projects/CIM_10_1_PKPD.sbproj');
    data_ext = {'/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo_2.mat'};
    data_ext1={'/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Morisada_3.mat'};
end
%% Load project 
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(cs, 'MaximumWallClock', 2.5)
variants = get(model, 'variants');
set(cs, 'Time','day')
set(cs, 'StopTime',365*2)
%% Rules
% Initial tumor volume
model.Rules(21,1).Active=true;
% Initial tumor cell number
model.Rules(1,1).Active=true;
% Initial PDL1_Tumor
model.Rules(18,1).Active=true;
% Initial CTLA4_CD8
model.Rules(19,1).Active=true;
% Initial CTLA4_Treg
model.Rules(23,1).Active=true;
% Initial PDL1_Immune
model.Rules(33,1).Active=true;
% Carrying capacity
model.Rules(54,1).Active=true;
% kdep_Treg
model.Reactions(54,1).Active=false;
% saturation_antiCTLA4
model.Rules(55,1).Active = false;

%% Parameter setup
parameters = {'Blood'; 'Peripheral';'CL_antiPDL1';'Q12';...
    'f_L'; 'V_vs';'kin_CD8'; 'kin_Treg'; 'kin_DC';'kin_MDSC';'kin_CD45';...
    'kpro_Tumor'; 'kpro_Tumor_Linear';'kill_CD8'; ...
    'K_MDSC';'K_IFNg';'K_CTLA4'; 'K_PDL1'; 'KDE_MDSC';'KDE_Treg';
    'V_intra';'PDL1_Tumor_ss';'PDL1_Immune_ss';'kel_CD45';'kel_DC';'kel_Effector';'kel_MDSC';...
    'kel_Naive';'kel_Treg';'kpro_CD8';'kdif';'kdeg_CTLA4';'kdeg_PDL1';'T_0'};
% Define outputs% Define outputs
groups_subset = {'MOC1_Control', 'MOC1_antiPDL1', 'MOC1_antiCTLA4', 'MOC1_antiCTLA4_antiPDL1'};
observables={'Tumor'   'TumorInter.Cancer'  'TumorInter.CD8_E' 'TumorInter.CD8_N'...
    'TumorInter.DC' 'TumorInter.GMDSC'...
    'TumorInter.CD4Foxp3' 'TumorInter.Ag'  'TumorInter.CD45'...
    'TumorInter.CTLA4_CD8' 'TumorInter.CTLA4_Treg'...
    'TumorInter.PDL1_Tumor' 'TumorInter.PDL1_Immune'...
    'TumorInter.antiPDL1' 'TumorInter.antiCTLA4'};

stateVar={'Tumor'};
doses = {'Blood.Dose_antiPDL1' 'Blood.Dose_antiCTLA4'};
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

PI.observablesPlot={'Tumor' 'Cancer' 'CD8_E' 'CD8_N' 'DC' 'GMDSC'...
    'CD4Foxp3' 'Ag' 'CD45' 'CTLA4_CD8' 'CTLA4_Treg' 'PDL1_Tumor'...
    'PDL1_Immune', 'Tumor_Inter_antiPDL1', 'Tumor_Inter_antiCTLA4'};
PI.observablesPlot={'Tumor' 'Cancer' 'CD8_E' 'CD8_N' 'DC' 'GMDSC'...
    'CD4Foxp3' 'Ag' 'CD45' 'CTLA4_{CD8}' 'CTLA4_{Treg}' 'PDL1_{Tumor}'...
    'PDL1_{Immune}', 'antiPDL1_{TumorInter}', 'antiCTLA4_{TumorInter}'};
PI.tspan = 0:1:(365*2);
% Get simulation function
[sim,PI.u]=initializePI(model,parameters,observables,PI,doses, 'Human_InitialValues','doseUnits', 'mole');
PI = getDoseRegimen(PI);

%% Hierarchical model simulation
PI.H = getHierarchicalStruct3(parameters(1:end),PI,'n_sigma', length(observables),...
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
PI.paramNames  =paramNames;
%% Load PIs and posterior parameter samples
PI_PK = load(strjoin({cd 'PK_mAb_ThreeComp/PI/CIM10/PI_CIM10_PK_1.mat'},'/'));
PI_MOC1 = load(strjoin({cd 'CIM/PI/CIM33/PI_CIM10_MOC1_PKPD_1_8_kinCD8.mat'},'/'));

Meta=[];
Meta(1).Struct = PI_PK;
Meta(1).Struct = PI_MOC1;
values=repelem({sim.Parameters.Value},length(Meta),1);
[Meta(1:end).Values]=values{:,:};

%% Select which dimensions to sample from 
N_pop =1e3;
[theta1,~] = getTheta2(Meta, parameters, N_pop, 'covStructure', 'dependent','fixedParams','uniform');
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
sampleIndx_pk=randsample(size(theta_pk,1),N_pop, true);
delta_pk = theta_pk(sampleIndx_pk,:)-mean(theta_pk);
Theta(:,indx_pk) =priorMu_pk+delta_pk(:,pk_order);

% PD parameters
[pd_indx, pd_order]=ismember(parameters,PD_params);
pd_order=pd_order(~pd_order==0);
theta_pd = Theta(:,pd_indx);
delta_pd = theta_pd-mean(theta_pd);
Theta(:,indx_pd) =priorMu_pd + delta_pd(:,pd_order);
histogram(exp(Theta(:,indx_pd)))

TV_range = [.5 1].^3*(4/3)*pi;
TV_0 = log((rand(size(Theta,1),1)*(TV_range(2)-TV_range(1))... % Initial T_0 (1e6 cells)
    +TV_range(1)));
% histogram(exp(TV_0))
T_0=log(exp(TV_0)/.00153.*exp(Theta(:,ismember(parameters,'V_intra'))));
Theta(:,ismember(parameters,'T_0')) = T_0;
histogram(exp(T_0))

table(parameters,exp(mean(Theta))')
%% 
plotBivariateMarginals_2(exp(Theta(:,:)),...
       'names',parameters,'interpreter', 'tex')
plotCorrMat(Theta, parameters)
plotCorrMat(exp(Theta), parameters)

%%
PI = globalSA2(sim,PI,observables,'samples',exp(Theta(1:end,:)),'nsamples',1e3,'time',1:14:PI.tspan(end),...
    'sigma', 1,'inputs', exp(mean(Theta(1:end,:))),'variation', 0.5);
paramRanking = plotPRCC3(PI,PI.paramNames(PI.H.PopulationParams),PI.observablesPlot,'time',...
    1:14:PI.tspan(end),'output', 'mean','kpi', 'sum','nparams',5,'ncol',4);
parameters_hat1 =unique(paramRanking{1:8,:});
%% Save
save('GSA_CIM10_Human_PKPD_1_8_kinCD8_MOC1.mat', 'PI')

