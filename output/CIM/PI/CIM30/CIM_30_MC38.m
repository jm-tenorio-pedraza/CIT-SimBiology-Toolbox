%% General setup for Kosinsky et al with k_pro as the parameter to vary for each individual
%% Search paths
clear all
warning off
opengl('save', 'software')
%%
if ispc
    addpath(genpath('\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox'))
    cd('\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\output\CIM\PI\CIM30')
    out = sbioloadproject('\Users\jmten\OneDrive\Dokumente\GitHub\sbio-projects\CIM_10.sbproj');

else
    addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
    cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/CIM/PI/CIM30')
    out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/sbio-projects/CIM_10.sbproj');

end
%% Load project 
% Extract model
model=out.m1;
variants = getvariant(model);
initialStruct = struct('name', {'MOC1';'MOC2';'MC38'}, 'initialValue', {5; 0.1; 0.1},...
    'variant', {variants(1); variants(3); variants(4)});

cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(cs, 'MaximumWallClock', 0.25)
set(cs, 'Time','day')
set(cs, 'StopTime',100)

sbioaccelerate(model, cs)
%% Parameter setup
parameters = {'kin_CD8'; 'kin_Treg';'kin_DC';...
    'kin_MDSC';'kpro_Tumor'; 'kpro_Tumor_Linear'; 'kill_CD8'; 
    'K_MDSC'; 'K_CTLA4'; 'K_PDL1'; 'S_L';'ks_PDL1_Tumor';...
    'ks_PDL1_Immune'; 'kill_Treg'; 'kdep_max'; 'kin_CD45'; 'K_CD45';...
    'PDL1_Tumor_ss'; 'PDL1_Immune_ss'; 'K_IFNg'; 'kdeg_PDL1'; 'kin_Tumor';...
    'K_Tumor'; 'ks_CTLA4_CD8'; 'ks_CTLA4_Treg'};
parameters = [parameters; 'T_0'];

% Define outputs% Define outputs
groups_subset = {'MC38_Control' 'MC38_Control_Immune'};
observables={'TV'  'CD8'  'Treg' 'DCm'...
    'MDSC' 'CD107a' 'PDL1_T' 'PDL1_I'};
stateVar={'Tumor'  'CD8' 'Treg' 'DC'...
    'GMDSC' 'CD107a_Rel' 'Tumor_PDL1_Rel' 'Myeloid_PDL1_Rel'};
doses = {'Blood.Dose_antiPDL1' 'Blood.Dose_antiCTLA4'};
%% Obtain data, simulation function and dose table
if ispc
    data_ext = {'\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\PI_Clavijo_2.mat'...
        '\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\PI_Morisada_3.mat'};
else
    data_ext = {'/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo_2.mat',...
        '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Morisada_3.mat'};
end
PI=getPIData4(data_ext, stateVar,groups_subset,'output', 'mean',...
    'responseGrouping', true, 'kineticGrouping', false, 'zeroHandling', 'imput');
try
PI.data = rmfield(PI.data, {'Kinetic'});
catch
end

PI.variableUnits={'Volume [mL]' 'Percentage [%]' 'Percentage [%]'  'Percentage [%]' ...
     'Percentage [%]'   'Relative units []' ...
    'Relative units []' 'Relative units []'};
PI.observablesFields = {'TV'  'CD8'  'Treg' 'DCm'...
    'MDSC' 'CD8_E' 'PDL1_T' 'PDL1_I'};
PI.normIndx = 6:8;
PI.model = 'CIM Control';
PI.observablesPlot={'TV' 'CD8' 'Treg' 'DCm'...
    'MDSC' 'CD107a' 'PDL1_T' 'PDL1_I'};
plotData(PI, PI.observablesPlot, 'responseGrouping', true, 'kineticGrouping', false)
% Get initial values
[PI.x_0, PI.variants] = getInitialValues([PI.data(:).Group],...
    initialStruct);

% Get simulation function
[sim,PI.u]=initializePI(model,parameters,observables,PI,doses, 'PKPD_Fit_MC38','doseUnits', 'mole');
clearvars   PI1 PI2 variants
% close all
%% Optimization setup
% Hierarchical structure
PI.H = getHierarchicalStruct(parameters(1:end-1),PI,'n_sigma', length(observables),...
    'rand_indx', [] , 'cell_indx',[], 'n_indiv', length(PI.u));
SigmaNames = getVarNames(PI, stateVar);
[beta, sigma_prior] = getVarValues([.01 .01 .001], [.01 .01 .001], [1 1 1], PI);
lb=([1e-2   1e-3   1e-3   1e-4    1e-3    1e-2 1e-5    1e-3 1e-2       1e-2   1e-1    1e-4    1e-4   1e-4    0.01   1e-3     1e-3 1e-3   1e-3    0.01    1e-2    10      1e-3    1e-2    1e-2])';
ub=([1e2    1e2    1e2    1e2     1e2     2e0  1e1     1    1e2        1e2    1e3     1e3     1e4    1e2     20     2e2      10   1e3    1e3     100     1e2     1e3     1       1e2     1e3])';

PI.par = getParamStruct2(sim,PI.H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior, 'ref', 'ones','LB', lb, 'UB', ub);

PI = assignPrior(PI);
% Log-ikelihood function
likelihood_fun=@(p)likelihood(exp(p),sim,PI,'censoring',false);
prior_fun=@(p)getPriorPDFMCMC2(exp(p),PI);


paramNames = getParamNames(PI,sim, observables);
PI.paramNames  =paramNames;
clearvars beta lb ub sigma_prior SigmaNames 

%% Objective function
try
    finalValues =log([PI.par(:).finalValue]);
catch
    finalValues =log([PI.par(:).startValue]);

end
% Obj function
obj_fun=@(x)(postLogLikelihood(x, prior_fun, likelihood_fun)*(-1));
tic
obj_fun((finalValues))
toc

%% Parameter selection of inter-cell line varying params

w = arrayfun(@(x) (finalValues(x.Index)), PI.H.CellParams,'UniformOutput', false);
w = (std(cell2mat(w'),[], 2));
[w,cell_indx] =sort(w,'descend');
cell_params = [ {PI.H.CellParams(:).name}'];

z = arrayfun(@(x) (finalValues(x.Index)), PI.H.IndividualParams,'UniformOutput', false);
z = (std(cell2mat(z'),[], 2));
[z,ind_indx] =sort(z,'descend');
ind_params = [{PI.H.IndividualParams(:).name}'];

z_Cell = arrayfun(@(x) std(PI.H.CellIndx'.*finalValues(x.Index),[],2),...
    PI.H.IndividualParams,'UniformOutput', false);
z_Cell = cellfun(@mean,z_Cell');
[z_Cell,indCell_indx] =sort(z_Cell,'descend');

table([cell_params(cell_indx); ind_params(ind_indx)], [w; z; ])
table([ind_params(indCell_indx)], z_Cell)
%% Save results
save('PI_CIM30_0.mat', 'PI')
if ispc
    load(strjoin({cd 'PI_CIM30_0.mat'},'\'),'PI')

else
load(strjoin({cd 'PI_CIM5_Control_Reduced_2_0.mat'},'/'),'PI')
end
%% save MCMC results
N_i='7';
save(strjoin({cd '/PI_CIM29_3_x_' N_i '.mat'},''), strjoin({'x' N_i},''))
save(strjoin({cd '/PI_CIM29_3_p_x_' N_i '.mat'},''), strjoin({'p_x' N_i},''))
save(strjoin({cd '/PI_CIM29_3_J' N_i '.mat'},''), strjoin({'J' N_i},''))
save(strjoin({cd '/PI_CIM29_3_n_id' N_i '.mat'},''), strjoin({'n_id' N_i},''))
save(strjoin({cd '/PI_CIM29_3_stepSize' N_i '.mat'},''), strjoin({'stepSize' N_i},''))
%%
load(strjoin({cd '/PI_CIM29_3_x_' N_i '.mat'},''))
load(strjoin({cd '/PI_CIM29_3_p_x_' N_i '.mat'},''))
load(strjoin({cd '/PI_CIM29_3_J' N_i '.mat'},''))
load(strjoin({cd '/PI_CIM29_3_n_id' N_i '.mat'},''))
load(strjoin({cd '/PI_CIM29_3_stepSize' N_i '.mat'},''))

for i=1:5
    load(strjoin({cd '/PI_CIM29_3_x_' num2str(i) '.mat'},''))
    load(strjoin({cd '/PI_CIM29_3_p_x_' num2str(i) '.mat'},''))
end

