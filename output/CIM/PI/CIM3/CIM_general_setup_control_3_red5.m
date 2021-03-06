%% General setup for Kosinsky et al with k_pro as the parameter to vary for each individual
%% Search paths
warning on
clear all
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/CIM/PI/CIM3')

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/CIM_3.sbproj');
% Extract model
model=out.m1;
variants = getvariant(model);
initialStruct = struct('name', {'MOC1';'MOC2'}, 'initialValue', {5; 0.1},...
    'variant', {variants(1); variants(2)});

cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-12);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-10);
set(cs, 'MaximumWallClock', 0.25)
%% Parameter setup
parameters = {'K_CTLA4'; 'kin_CD8';'kill_CD8';...
    'kpro_Tumor';'K_MDSC'; 'kel_Effector';'KDE_MDSC'; 'K_CD8'; 'kin_Treg';...
    'kill_Treg';'K_PDL1'};
parameters = [parameters; 'T_0'];

% Define outputs% Define outputs
groups_subset = {'MOC1_Control', 'MOC1_Control_Mean', 'MOC2_Control',...
    'MOC2_Control_Mean'};
observables={'TV'  'CD8' 'CD107a' 'Treg' 'DCm'...
    'MDSC' 'PDL1_T' 'PDL1_I'};
stateVar={'Tumor'  'CD8' 'CD107a' 'Treg' 'DC'...
    'GMDSC'...
    'Tumor_PDL1_Rel' 'Myeloid_PDL1_Rel'};
doses = {'Blood.Dose_antiPDL1' 'Blood.Dose_antiCTLA4'};
%% Obtain data, simulation function and dose table

PI=getPIData2('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat',...
    stateVar,groups_subset,observables,'zeroAction', 'omit','mergePhenotypes',...
    true,'output', 'mean','maxIIV', false);
PI.variableUnits={'Volume [mL]' 'Percentage [%]' 'Percentage [%]'  'Percentage [%]' ...
     'Percentage [%]'   'Percentage [%]' ...
    'Relative units []' 'Relative units []'};
PI.normIndx = 7:8;
PI.model = 'CIM Control';
PI.observablesPlot={'Tumor volume' 'CD8+ T-cells' 'CD107a+CD8+ T-cells' 'Treg' 'DCm'...
    'MDSC' 'PDL1_T' 'PDL1_I'};

% Get initial values
[PI.x_0, PI.variants] = getInitialValues([PI.data(:).Group],...
    initialStruct);

% Get simulation function
[sim,PI.u]=initializePI(model,parameters,observables,PI,doses, 'MOC1','doseUnits', 'mole');

%% Optimization setup
% Hierarchical structure
cell_indx = [3 4];
PI.H = getHierarchicalStruct(parameters(1:end-1),PI,'n_sigma', length(observables),...
    'rand_indx', [], 'cell_indx',cell_indx, 'n_indiv', length(PI.u));
if ~isempty(PI.H.IndividualParams(1).Index)
        indivSigmaNames=arrayfun(@(x)strjoin({'omega', x.name}, '_'),PI.H.IndividualParams,'UniformOutput',false)';
else
    indivSigmaNames = [];
end
if ~isempty(PI.H.CellParams(1).Index)
    cellSigmaNames=arrayfun(@(x)strjoin({'psi', x.name}, '_'),PI.H.CellParams,'UniformOutput',false)';
else
    cellSigmaNames = [];
end
try
SigmaNames = [cellSigmaNames; indivSigmaNames];
SigmaNames(end+1:end+length(observables),1) =  cellfun(@(x) strjoin({'sigma', x}, '_'),...
    observables','UniformOutput', false);
catch
    SigmaNames=cellfun(@(x) strjoin({'sigma', x}, '_'),...
    observables','UniformOutput', false);
end
% Generating PI
alpha = [repelem(0.01, length([PI.H.CellParams(:).OmegaIndex]),1);...
    repelem(0.01, length([PI.H.IndividualParams(:).OmegaIndex]),1);...
    repelem(0.001, length(setdiff(PI.H.SigmaParams, [PI.H.CellParams(:).OmegaIndex ...
    PI.H.IndividualParams(:).OmegaIndex])),1)];
beta = [repelem(0.01, length([PI.H.CellParams(:).OmegaIndex]),1);...
    repelem(0.01, length([PI.H.IndividualParams(:).OmegaIndex]),1);...
    repelem(0.001, length(setdiff(PI.H.SigmaParams, [PI.H.CellParams(:).OmegaIndex ...
    PI.H.IndividualParams(:).OmegaIndex])),1)];
sigma_prior= [ repelem(1,length(PI.H.PopulationParams), 1);...
    repelem(1, length([PI.H.CellParams(:).Index]),1);
     repelem(1, length([PI.H.IndividualParams(:).Index]),1);...
    alpha];
PI.par = getParamStruct2(sim,PI.H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior, 'ref', 'ones');
try
    finalValues =log([PI.par(:).finalValue]);
catch
    finalValues =log([PI.par(:).startValue]);

end

% Log-ikelihood function
likelihood_fun=@(p)likelihood(exp(p),sim,PI,'censoring',false);
prior_fun=@(p)(createPriorDistribution3(exp(p),PI,PI.H,'type',{'uniform/normal/inverse gamma/inverse gamma'}));

% Residuals 
residuals_fn = @(x) getResiduals(exp(x),@(x)sim(x,PI.tspan(end),PI.u,PI.tspan),PI,...
    @(x)getPhi2(x,PI.H,length(PI.u),'initialValue',PI.x_0),exp(finalValues(end-length(observables)+1:end)),PI.normIndx);
% Parameter names for plots
% paramNames = ['kin_{CD8}' '\eta_{K_{CD8}}' '\eta_{kpro_{Tumor}}' 'kill_{max}' 'KDE_{Treg}'...
%     'KDE_{MDSC}' 'K_{pro}' 'K_{el}' 'kpro_{Tumor_{MOC1}}' 'kpro_{Tumor_{MOC2}}'...
%     {PI.par([PI.H.IndividualParams(:).Index]).name} '\lambda_{kpro_{Tumor}}'...
%     '\omega_{kin_{CD8}}' '\sigma_{TV}' '\sigma_{CD8}' '\sigma_{CD107a}' '\sigma_{Treg}'...
%     '\sigma_{DC}' '\sigma_{MDSC}' '\sigma_{PDL1_T}' '\sigma_{PDL1_I}'];
paramNames = getParamNames(PI,sim, observables);
%% Objective function

% Obj function
obj_fun=@(x)(likelihood_fun(x)*(-1)+prior_fun(x)*(-1));
tic
obj_fun((finalValues))
toc


%% Save results
save('PI_CIM_Control_4_red.mat', 'PI')
load(strjoin({cd 'PI_CIM_Control_4_red.mat'},'/'),'PI')


