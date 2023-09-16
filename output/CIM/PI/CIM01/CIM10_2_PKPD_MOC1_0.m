%% General setup
%% Search paths
clear all
warning off
opengl('save', 'software')
%%
if ispc
    wd='\Users\jmten\OneDrive\Dokumente\GitHub\';
else
    wd='/Users/migueltenorio/Documents/GitHub/';
end
    addpath(genpath(strjoin({wd 'CIT-SimBiology-Toolbox'}, '')))
    cd(strjoin({wd 'CIT-SimBiology-Toolbox/output'},''))
    out=sbioloadproject(strjoin({wd 'sbio-projects/CIM_10_1_PKPD.sbproj'},''));
    data_ext=strjoin({wd 'CIT-SimBiology-Toolbox/data/PI_Clavijo_2.mat'},'');
%% Load project 
% Extract model
model=out.m1;
variants = getvariant(model);
initialStruct = struct('name', {'MOC1';'MOC2';'MC38'}, 'initialValue', {5; 0.1; 0.1},...
    'variant', {variants(1); variants(3); variants(4)});

cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(cs, 'MaximumWallClock', .5)
set(cs, 'Time','day')
set(cs, 'StopTime',100)

    %% Rule setup
% Initial tumor volume
model.Rules(21,1).Active=true;
% Initial tumor cell number
model.Rules(1,1).Active=true;
% Initial PDL1_Tumor
model.Rules(18,1).Active=true;
% Initial CTLA4_CD8
model.Rules(19,1).Active=true;
% Initial CTLA4_Treg
model.Rules(22,1).Active=true;
% Initial PDL1_Immune
model.Rules(32,1).Active=true;
% Carrying capacity
model.Rules(50,1).Active=false;
% Observers
model.Rules(2,1).Active=true;
model.Rules(3,1).Active=true;
model.Rules(4,1).Active=true;
model.Rules(5,1).Active=true;
model.Rules(6,1).Active=true;
model.Rules(9,1).Active=true;
model.Rules(10,1).Active=true;
model.Rules(12,1).Active=true;
% saturation_antiCTLA4
model.Rules(51,1).Active = false;

%% Parameter setup
parameters=arrayfun(@(x)x.Name,model.Parameters,'UniformOutput',false);
constantIndx=arrayfun(@(x)x.ConstantValue==1,model.Parameters,'UniformOutput',true);
zeroIndx=arrayfun(@(x)x.Value~=0,model.Parameters,'UniformOutput',true);
dimIndx=arrayfun(@(x)~ismember({x.ValueUnits},'dimensionless'),model.Parameters,'UniformOutput',true);
parameters=parameters(and(and(constantIndx,zeroIndx),dimIndx));
excludedParams={'cells_per_mole' 'T_0' 'vol_Tumor' 'Tumor_0' 'kdep_max' 'antiLy6G_dose1'...
    'antiLy6G_lag' 'antiLy6G_on' 'antiLy6G_off' 'ID_antiCTLA4_Blood_free_min' ...
    'ID_PDL1_Blood_free_min' 'kdeg_antiPDL1_ADA' 'ADA_onset' 'BW' 'cell' 'V_int' 'GMDSC_antiLy6G' 'IARkill_w_PDL1'...
    'kon_antiCTLA4' 'kon_antiPDL1' 'kon_mFcRn_mIgG2b' 'kon_mFcRn_rIgG2b' ...
    'k12' 'k21'  'EC50_antiCTLA4' 'TotalBloodVol'};
exclIndx=ismember(parameters,excludedParams);
parameters=parameters(~exclIndx);

parameters=[sort(parameters(1:end-1));'T_0'];

% Define outputs% Define outputs
groups_subset = {'MOC1_Control', 'MOC1_Control_Mean'};
observables={'Tumor'  'CD8'  'Treg' 'DCm'...
    'MDSC' 'CD8_E' 'PDL1_T' 'PDL1_I'};
stateVar={'Tumor'  'CD8' 'Treg' 'DC'...
    'GMDSC' 'CD107a_Rel' 'Tumor_PDL1_Rel' 'Myeloid_PDL1_Rel'};
doses = {'Blood.Dose_antiPDL1' 'Blood.Dose_antiCTLA4'};
%% Load PI data structure for each treatment group
PI=getPIData4({data_ext}, stateVar,groups_subset,'output', 'mean',...
    'responseGrouping', true, 'kineticGrouping', false, 'zeroHandling',...
    'replace','logTransform',false);
PI.n_data=sum(cellfun(@(x)sum(sum(~isnan(x))), {PI.data(:).dataValue}));
cellgroup=repelem({'MOC1'},2);
[PI.data(1:end).Cell]=cellgroup{:,:};
PI.data(1).Group='MOC1_Control';
PI.data(2).Group='MOC1_Control';
%%
PI.variableUnits={'Volume [mL]' 'Percentage [%]' 'Percentage [%]'  'Percentage [%]' ...
     'Percentage [%]'   'Relative units []' ...
    'Relative units []' 'Relative units []'};
PI.observablesFields = {'TV'  'CD8'  'Treg' 'DCm'...
    'MDSC' 'CD8_E' 'PDL1_T' 'PDL1_I'};
PI.normIndx = 6:8;
PI.model = 'MOC1 Calibration';
PI.observablesPlot={'TV' 'CD8' 'Treg' 'DCm'...
    'MDSC' 'CD107a' 'PDL1_T' 'PDL1_I'};
plotData(PI, PI.observablesPlot, 'responseGrouping', true, 'kineticGrouping', false)
% Get initial values
[PI.x_0, PI.variants] = getInitialValues({PI.data(:).Group},...
    initialStruct);
%% Get simulation function
[sim,PI.u]=initializePI(model,parameters,observables,PI,doses, 'PKPD_SA_MOC1_0','doseUnits', 'mole');

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

% Log-ikelihood function
PI = assignPrior(PI,'sigmaDist', 'JP');
% Log-ikelihood function
likelihood_fun=@(p)likelihood(exp(p),sim,PI,'censoring',false,...
    'logTransform',true,'errorModel','additive','constantVar',.01,'indivData',false);
prior_fun=@(p)getPriorPDFMCMC2(exp(p),PI);
paramNames = getParamNames(PI,sim, observables);
PI.paramNames  =paramNames;
clearvars beta lb ub sigma_prior SigmaNames 
%%
try
    finalValues =log([PI.par(:).finalValue]);
catch
    finalValues =log([PI.par(:).startValue]);

end

%% Simulation output
simTime = unique([PI.tspan', 1:PI.tspan(end)]);
PI=getOutput(PI,@(p)sim(p,PI.tspan(end),PI.u, simTime),exp(finalValues),...
    @(p)getPhi2(p,PI.H,length(PI.u)),...
    PI.normIndx,PI.H,'output', 'PI', 'simTime', simTime,'logTransform',false, 'errorModel','multiplicative');
PI.AIC = 2*length(PI.par)-2*likelihood_fun(finalValues)*(1);
%% Plotting output
figure('Position', [10 10 1.5e3 1e3])
ncol = ceil(sqrt(length(observables)));
nrow = ceil(length(observables)/ncol);
for i=1:length(observables)
subplot(nrow,ncol,i)
 plotSimOutput(PI,i,'all', false, 'indiv', false, 'addErrorVar', true,...
     'newFig', false, 'TimeUnit', 'days','indivData',false)
  set(gca, 'YScale','log')
% % ylim([0 .1])
end

%%
save(strjoin({cd 'CIM/PI/CIM33/SA_CIM10_PKPD_MOC1.mat'},'/'), 'PI')


