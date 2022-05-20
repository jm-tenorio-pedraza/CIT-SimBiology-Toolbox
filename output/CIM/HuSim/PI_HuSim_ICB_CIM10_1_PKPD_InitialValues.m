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
%% Load project 
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-11);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-9);
set(cs, 'MaximumWallClock', 2.5)
variants = get(model, 'variants');
set(cs, 'Time','day')
set(cs, 'StopTime',365*3)
model.Rules(21,1).Active=false;
model.Rules(1,1).Active=false;
model.Rules(18,1).Active=false;
model.Rules(19,1).Active=false;
model.Rules(23,1).Active=false;
model.Rules(33,1).Active=false;

%% Setting up parameters, data and simulations
parameters = {'Blood'; 'Peripheral';'CL_antiPDL1';'Q12';...
    'f_L'; 'V_vs';'kin_CD8'; 'kin_Treg'; 'kin_DC';'kin_MDSC';'kin_CD45';...
    'kpro_Tumor'; 'kpro_Tumor_Linear';'kill_CD8'; ...
    'K_MDSC';'K_IFNg';'K_CTLA4'; 'K_PDL1'; 'KDE_MDSC';'KDE_Treg';
    'V_intra';'PDL1_Tumor_ss';'PDL1_Immune_ss';'kel_CD45';'kel_DC';'kel_Effector';'kel_MDSC';...
    'kel_Naive';'kel_Treg'};
% Define outputs
groups_subset = {'MOC1_Control', 'MOC1_antiPDL1', 'MOC1_antiCTLA4', 'MOC1_antiCTLA4_antiPDL1'};
observables={'Tumor'   'TumorInter.Cancer'  'TumorInter.CD8_E' 'TumorInter.CD8_N'...
    'TumorInter.DC' 'TumorInter.GMDSC'...
    'TumorInter.CD4Foxp3' 'TumorInter.Ag'  'TumorInter.CD45'...
    'TumorInter.CTLA4_CD8' 'TumorInter.CTLA4_Treg'...
    'TumorInter.PDL1_Tumor' 'TumorInter.PDL1_Immune'...
    'TumorInter.antiPDL1' 'TumorInter.antiCTLA4'};

parameters=[parameters; observables(1:end-2)'];
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
% PI.observablesPlot={'Tumor'  'CTLA4_CD8' 'PDL1_Tumor' 'PDL1_Immune'...
%     'CTLA4_Treg' 'Cancer' 'CD8_N' 'DC' 'GMDSC'...
%     'CD4Foxp3' 'Ag' 'CD8_E' 'CD45' 'Tumor_Inter_antiPDL1' 'Tumor_Inter_antiCTLA4'};
PI.observablesPlot={'Tumor' 'Cancer' 'CD8_E' 'CD8_N' 'DC' 'GMDSC'...
    'CD4Foxp3' 'Ag' 'CD45' 'CTLA4_CD8' 'CTLA4_Treg' 'PDL1_Tumor'...
    'PDL1_Immune', 'Tumor_Inter_antiPDL1', 'Tumor_Inter_antiCTLA4'};

PI.tspan = 0:14:(365*3);

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

%% Load PIs with initial values and posterior samples
%% Save output
for i=1:16
    num_i = num2str(i);
    load(strjoin({cd '/CIM/HuSim/CIM10/PI' num_i '.mat'},''),strjoin({'PI', num_i}, ''))
end

%% Extendend theta
HumanPI=[];
CancerIndx=ismember(observables(1:end-2),'TumorInter.Cancer');
TumorIndx=ismember(observables(1:end-2),'Tumor');

HumanPI(1).Theta=[PI1.theta(:,1:end-1) log(PI1.initialValues)];
HumanPI(2).Theta=[PI2.theta(:,1:end-1) log(PI2.initialValues)];
HumanPI(3).Theta=[PI3.theta(:,1:end-1) log(PI3.initialValues)];
HumanPI(4).Theta=[PI4.theta(:,1:end-1) log(PI4.initialValues)];
HumanPI(5).Theta=[PI5.theta(:,1:end-1) log(PI5.initialValues)];
HumanPI(6).Theta=[PI6.theta(:,1:end-1) log(PI6.initialValues)];
HumanPI(7).Theta=[PI7.theta(:,1:end-1) log(PI7.initialValues)];
HumanPI(8).Theta=[PI8.theta(:,1:end-1) log(PI8.initialValues)];
HumanPI(9).Theta=[PI9.theta(:,1:end-1) log(PI9.initialValues)];
HumanPI(10).Theta=[PI10.theta(:,1:end-1) log(PI10.initialValues)];
HumanPI(11).Theta=[PI11.theta(:,1:end-1) log(PI11.initialValues)];
HumanPI(12).Theta=[PI12.theta(:,1:end-1) log(PI12.initialValues)];
HumanPI(13).Theta=[PI13.theta(:,1:end-1) log(PI13.initialValues)];
HumanPI(14).Theta=[PI14.theta(:,1:end-1) log(PI14.initialValues)];
HumanPI(15).Theta=[PI15.theta(:,1:end-1) log(PI15.initialValues)];
HumanPI(16).Theta=[PI16.theta(:,1:end-1) log(PI16.initialValues)];


%% Save output
save(strjoin({cd '/CIM/HuSim/CIM10/HumanPI_initialValues.mat'},''),'HumanPI')