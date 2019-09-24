%% General setup for CIM model with anti-PDL1, anti-CTLA4 and combination therapy 
clear all
warning off
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))

cd('/Users/migueltenorio/Documents/MATLAB/SimBiology/CIM')

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/CIM_2.sbproj');
% Extract model
model=out.model;
variants=getvariant(model);
MOC1=variants(3);
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(cs, 'MaximumWallClock', 0.25)


%% Create function handle for simulations
% All parameters 
[name,I] = sort(get(model.Parameters, 'Name'));
value = cell2mat(get(model.Parameters, 'Value'));
value = value(I);
% Define parameters to estimate
parameters=name(value>0);
exclude_parameters = {'Avogadro' 'omega' 'vol_Tumor' 'CD25_0' 'CTLA4_CD8_0'...
    'CTLA4_Treg_0' 'KD_CD25' 'KD_antiCTLA4' 'KD_antiPDL1' 'koff_CD25' 'koff_antiCTLA4' ...
    'koff_antiPDL1'};
fixed_parameters = {'PDL1_Immune_0' 'PDL1_Tumor_0'...
    'k12' 'k23' 'k32' 'ka' 'ks_IFNg' 'ks_IL2' 'Tumor_P' 'Tumor_P_0'...
    'kdeg_CTLA4' 'kdeg_PDL1' 'kel_Debris' 'f1' 'f2' 'f3' };
uncertain_parameters = {'K_IFNg' 'K_IL2'  'kpro_Naive_max' ...
    'kel_DC' 'kel_Effector' 'kel_MDSC' 'kel_Naive' 'kel_Treg' 'kel_Tumor'};
parameters = setdiff(parameters, [exclude_parameters fixed_parameters ]);
%% Define outputs for tumor volume under different treatment conditions
observables={'TV'};
     stateVar = {'Tumor'};
 groups_subset = {'MOC1_Control' 'MOC1_Control_Mean' 'MOC1_antiPDL1' 'MOC1_antiCTLA4' 'MOC1_antiCTLA4_antiPDL1'};
 doses = {'Dose_antiPDL1' 'Dose_antiCTLA4'};
 dose = doses;
%% Define settings for control 
 observables={'TV' 'CD8_logit' 'CD107a_logit' 'MDSC_logit' 'DC_logit' 'Treg_logit' 'PDL1_Tumor_Rel' 'PDL1_Immune_Rel' };

 groups_subset = {'MOC1_Control' 'MOC1_Control_Mean'};
stateVar={'Tumor' 'CD8' 'CD107a' 'Treg' 'DC' 'GMDSC' 'Tumor_PDL1_Rel' 'Myeloid_PDL1_Rel'};

dose = {'Dose_antiPDL1'};
doses={'Control'};
%% Run pre-processing of data
run('Clavijo_Group_Pre_Processing.m')
%%
save('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/CIM/phat_TV.mat','par_hat')
