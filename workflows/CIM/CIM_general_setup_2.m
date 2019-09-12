%% General setup for CIM model with anti-PDL1, anti-CTLA4 and combination therapy 
clear all
warning off
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))

cd('/Users/migueltenorio/Documents/MATLAB/SimBiology/CIM')

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/CIM.sbproj');
% Extract model
model=out.m1;
variants=getvariant(model);
MOC1=variants(2);
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(cs, 'MaximumWallClock', 10)


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
fixed_parameters = {'K_IFNg' 'K_IL2' 'K_antiCTLA4' 'PDL1_Immune_0' 'PDL1_Tumor_0'...
    'k12' 'k23' 'k32' 'ka' 'kdep_max' 'ks_IFNg' 'ks_IL2' 'Tumor_P' 'Tumor_P_0'...
    'kdeg_CTLA4' 'kdeg_PDL1' ...
    'kdif_max' 'kel_Debris' 'kel_Effector' 'kel_Naive' 'kel_Tumor' 'kel_DC'...
    'kel_MDSC' 'kel_Treg' 'kpro_Naive_max' };
parameters = setdiff(parameters, [exclude_parameters fixed_parameters]);
% Define outputs
observables={'TV' 'CD8' 'CD107a' 'MDSC' 'DCm' 'Treg' 'PDL1_Tumor_Rel' 'PDL1_Immune_Rel' };

stateVar={'Tumor' 'CD8' 'CD107a' 'Treg' 'DC' 'GMDSC'...
    'Tumor_PDL1_Rel' 'Myeloid_PDL1_Rel'};
groups_subset = {'MOC1_Control' 'MOC1_Control_Mean' 'MOC1_antiPDL1' 'MOC1_antiCTLA4' 'MOC1_antiCTLA4_antiPDL1'};

run('Clavijo_Group_Pre_Processing.m')
doses = {'Dose_antiPDL1' 'Dose_antiCTLA4'};
