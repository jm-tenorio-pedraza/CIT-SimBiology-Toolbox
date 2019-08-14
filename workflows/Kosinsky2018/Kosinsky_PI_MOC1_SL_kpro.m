%% Search paths
warning off
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))

cd('/Users/migueltenorio/Documents/MATLAB/SimBiology/Kosinsky')
%% Load previous results
load('models_SL_kpro.mat')
load('logL_SL_kpro.mat')
%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/Kosinsky.sbproj');

% Extract model
kosinsky=out.m1;
cs=kosinsky.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-14);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-13);

% Extract predifined parameter contents
variants=getvariant(kosinsky);
%% load data and previous results
load('Kosinsky_data.mat')
data_subset=Kosinsky_data([1:2:23 24],:);
PI.data=data_subset;
%% Create function handle for simulations
% Define parameters to estimate
parameters={'r' 'TV_max' 'e_Td' 'k_LN' 'K_TCD' 'K_pdl' 'S_R' 'S_L' 'k_pro'};

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
    'UseParallel', false);

% Create cell of doses
antiPDL1_table=table([7 12 17]', [0.0013 0.0013 0.0013]', [0 0 0]',...
    'VariableNames',{'Time' 'Amount' 'Rate'});
antiPDL1_table.Properties.VariableUnits={'day' 'micromole' 'micromole/second'};
control_table=table([7 12 17]', [0 0 0]', [0 0 0]',...
    'VariableNames',{'Time' 'Amount' 'Rate'});
control_table.Properties.VariableUnits={'day' 'micromole' 'micromole/second'};

u=[repelem({control_table},6)'; repelem({antiPDL1_table},6)'; repelem({control_table},1)'];

%% Optimization setup
% mu_prior
%       r      TV_max   e_Td    k_LN    K_TCD   K_pdl    S_R     S_L                 k_pro              b_0
mu=     [0.4   10       10      1e-2    0.2    .1       0.03    repelem(.0089,13)   repelem(1,13)       repelem(0.1,6)];
% sigma_prior
sigma=  [2     3e0      3e0     3e0     3e0     3       3       repelem(3,13)       repelem(3,13)       repelem(1.1,6)];
% Initial value
p0=     [0.1   10       3.5     0.01    0.2     0.02    0.03    repelem(.0089,13)   repelem(1,13)       repelem(0.1,6)];
% Parameter boundaries
p_lb=   [0.01  1e0      1e-3    1e-4    1e-4    1e-4    1e-4    repelem(1e-4,1,13)  repelem(1e-3,1,13)  repelem(1e-3,1,6)];
p_ub=   [10    1e4      1e3     1e0     1e1     1e1     1e1     repelem(1e1,1,13)   repelem(4,1,13)     repelem(10,1,6)];

% Generating PI
PI.n_data=sum(cellfun(@(x)sum(sum(~isnan(x))),{PI.data.dataValue},'UniformOutput',true));
sigmaNames={'b_TV'; 'b_CD8'; 'b_CD107a'; 'b_DC'; 'b_MDSC'; 'b_PDL1'};
PI.par=getParamStruct(sim,1:7, 8:9,13,[p_lb' p0' p_ub' mu' sigma'],sigmaNames);
% Hierarchical structure
H.PopulationParams=1:7;
H.IndividualParams=struct('name', {'S_L'; 'k_pro'}, 'Index', {8:20; 21:33});
H.SigmaParams=34:39;

% Residuals function
residuals_fun=@(p)getResiduals(p,@(x)sim(x,100,u,1:1:100),PI,...
    @(x)getPhi2(x,H,size(PI.data,1)),(@(x)getCovariance(x,H)),4:6);

% Log-ikelihood function
likelihood_fun=@(p)sum(residuals_fun(exp(p))*(-1));
prior_fun=createPriorDistribution2(PI,H,'type','normal');

% Obj function
obj_fun=@(x)(likelihood_fun(x)*(-1)+prior_fun(x)*(-1));

% Optimizer options
options_fminsearch=optimset('Display','iter','MaxFunEvals', 1e4, 'MaxIter',1e4);
options_anneal.Verbosity=2;
options_anneal.InitTemp=2;
options = optimoptions('simulannealbnd','PlotFcns',...
          {@saplotbestx,@saplotbestf,@saplotx,@saplotf},'InitialTemperature', 10, 'MaxFunctionEvaluations', 1e4);
%% Optimizing
% Nelder-Mead
p_hat=fminsearch(obj_fun,log(p0),options_fminsearch);

% Simulated annealing
finalValues=anneal(obj_fun,p_hat,options_anneal);
tic
finalValues=simulannealbnd(obj_fun,finalValues,[],[],options);
toc
% Simulation output
PI=getOutput(PI,@(p)sim(p,100,u,1:1:100),exp(finalValues),...
    @(p)getPhi2(p,H,length(u)), 4:6);

% Plotting tumor volume
for i=1:6
plotSimOutput(PI.data,i)
legend(observables(i))
end

finalValue=num2cell(exp(finalValues'));
[PI.par(1:end).finalValue]=finalValue{:,:};

%% MCMC
% Initial walkers
w0=[finalValues;randn(300,size(finalValues,2))*0.1+repmat(finalValues,300,1)];
logL=rowfun(obj_fun,table(w0));
[~,I]=sort(logL{:,:});

w0=w0(I(1:size(w0,2)*2),:);

% MCMC sampler
tic
[models,logP]=gwmcmc(w0',{prior_fun likelihood_fun},1e6, 'StepSize', 1.2,...
'ThinChain',1,'BurnIn',0);
toc

tic
[models2,logP2]=gwmcmc(models(:,:,end),{prior_fun likelihood_fun},1e6, 'StepSize', 1.1,...
'ThinChain',1,'BurnIn',0);
toc
models_array=cat(3,models,models2);
logL=cat(3,logP,logP2);

% Plotting thinned out samples
plotMCMCDiagnostics(models(:,:,:),logP(:,:,:),PI)

% Bivariate posterior thinned out samples
plotBivariateMarginals_2(models_array(6e5:100:end,:),PI)
%% Posterior predictions
simFun=@(x)getOutput(PI,@(p)sim(p,100,u,1:1:100),x,...
    @(p)getPhi2(p,H,length(u)), 4:6);
tic
PI=getPosteriorPredictions(exp(models_array(6e5:400:end,:)),PI,simFun,observables);
toc
PI=getCredibleIntervals(PI,observables, exp(models_array(6e5:400:end,:)),H);
plotPosteriorPredictions(PI,observables)
                    
%% Save results
save('PI_SL_kpro.mat', 'PI')
save('models_SL_kpro.mat','models')
save('logP_SL_kpro.mat','logP')
