%% Search paths
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))

cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/workflows/Kosinsky2018')
warning off
%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/Kosinsky.sbproj');

% Extract model
kosinsky=out.m1;
cs=kosinsky.getconfigset;
% Extract predifined parameter contents
variants=getvariant(kosinsky);
%% load data and previous results
load('data.mat')
load('PI.mat')
data_subset=data([1:2:23 24],:);
%% Create function handle for simulations
% Define parameters to estimate
parameters={'r' 'TV_max' 'e_Td' 'k_LN' 'k_pro' 'd_0', 'K_pdl' 'S_R' 'S_L'};

% Define outputs
observables={'TV' 'CD8' 'CD107a' 'DCm' 'ISC' 'PDL1'};

% Define dose
control = sbiodose('rd');
control.TargetName = 'Dose_antiPDL1';
control.TimeUnits = 'day';
control.AmountUnits = 'micromole';
control.RateUnits = 'micromole/second';

% Create simFunction object
sim=createSimFunction(kosinsky,parameters,observables, [control],...
    'UseParallel', true);

% Create cell of doses
antiPDL1_table=table([7 12 17]', [0.0013 0.0013 0.0013]', [0 0 0]',...
    'VariableNames',{'Time' 'Amount' 'Rate'});
antiPDL1_table.Properties.VariableUnits={'day' 'micromole' 'micromole/second'};
control_table=table([7 12 17]', [0 0 0]', [0 0 0]',...
    'VariableNames',{'Time' 'Amount' 'Rate'});
control_table.Properties.VariableUnits={'day' 'micromole' 'micromole/second'};

u=[repelem({control_table},6)'; repelem({antiPDL1_table},6)'; repelem({control_table},1)'];

%% Optimization setup

% Initial value
p0=([0.2; 3.5; 1; 0.03; 1.5; 0.005; 0.001; 0.03;...
    repelem(sim.Parameters.Value(9),13)'; repelem(0.1,6)'])';
% Parameter boundaries
p_lb=[1e-2 1e-0 1e-3 1e-4 1e-3 1e-4  1e-4 1e-8 repelem(1e-8,1,13) repelem(1e-3,1,6)];
p_ub=[1e0 1e3 1e4 1e0 1e1 1e2 1e0 1e1 repelem(1e0,1,13), repelem(10,1,6)];
% Generating PI
paramNames=[sim.Parameters.Name(1:8); repelem(sim.Parameters.Name(9),13)';...
    {'b_TV'; 'b_CD8'; 'b_CD107a'; 'b_DC'; 'b_MDSC'; 'b_PDL1'}];
PI.par=struct('Name', paramNames,'minValue',num2cell((p_lb')), 'maxValue',...
    num2cell(p_ub'),'startValue', num2cell(p0'));
PI.data=data_subset;
PI.n_data=sum(cellfun(@(x)sum(sum(~isnan(x))),{PI.data.dataValue},'UniformOutput',true));

% Hierarchical structure
H.PopulationParams=1:8;
H.IndividualParams=9:21;
H.SigmaParams=22:27;

% Objective function
residuals_fun=@(p)getResiduals(p,@(x)sim(x,100,u,1:1:100),PI,...
    @(x)getPhi(x,H,size(PI.data,1)),(@(x)getCovariance(x,H)));
obj_fun=@(p)sum(residuals_fun(exp(p)));

% Optimizer options
options_fminsearch=optimset('Display','iter');
options_anneal.Verbosity=2;
options_anneal.InitTemp=2;

%% Optimizing
% Nelder-Mead
p_hat=fminsearch(obj_fun,finalValues,options_fminsearch);
% Simulated annealing
finalValues=anneal(obj_fun,p_hat,options_anneal);

% Simulation output
PI.data=getOutput(data_subset,@(p)sim(p,100,u,1:1:100),exp(finalValues),...
    @(p)getPhi(p,H,length(u)));

% Plotting tumor volume
for i=1:6
plotSimOutput(PI.data,i)
end
finalValue=num2cell(exp(finalValues'));
[PI.par(1:end).finalValue]=finalValue{:,:};


%% MCMC
% Initial walkers
w0=[finalValues;randn(300,size(finalValues,2))*0.1+repmat(finalValues,300,1)];
logL=rowfun(obj_fun,table(w0));
[~,I]=sort(logL{:,:});

w0=w0(I(1:size(w0,2)*2),:);
% Log-ikelihood function
likelihood_fun=@(p)sum(residuals_fun(exp(p))*(-1));
prior_fun=createPriorDistribution(PI,H,'type','normal');

% MCMC sampler
tic
[models,logP]=gwmcmc(w0',{prior_fun likelihood_fun},1e6, 'StepSize', 1.2,...
'ThinChain',1,'BurnIn',0);
toc
models_array=models(:,:)';
logL=sum(logP(:,:)',2);
plotTrace(logL,PI)
plotTrace(models_array,PI)
plotTrace(exp(models_array(7e5:100:end,:)),PI)
tightfig
colIndx=num2cell(repelem(1,12));
colIndx(13)={2:6};
[PI.data(1:end).colIndx]=colIndx{:,:};
[data]=getMCMCOutput(PI.data,@(p)sim(p,100,u,1:1:100),exp(models_array(9e5:1000:end,:)),@(p)getPhi(p,H,length(u)));
%% Save results

save('PI.mat', 'PI')
save('Kosinsky_MOC1_SL_models.mat','models')
save('Kosinsky_MOC1_SL_logP.mat','logP')
