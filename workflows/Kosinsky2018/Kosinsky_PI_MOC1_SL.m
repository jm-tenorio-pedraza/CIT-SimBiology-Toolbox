%% Search paths
warning off
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/MATLAB/SimBiology/Kosinsky/output/PI_SL')

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/Kosinsky.sbproj');

% Extract model
kosinsky=out.m1;
cs=kosinsky.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(cs, 'MaximumWallClock', 0.2)

%% load data and previous results
load('Kosinsky_data.mat')
data_subset=Kosinsky_data([1:2:23 24],:);
PI.data=data_subset;
PI.variableUnits={'Volume [mL]' 'Percentage [%]' 'Percentage [%]' 'Relative units []'...
    'Relative units []' 'Relative units []'};

load('models_Kosinsky_MOC1_SL.mat')
load('PI_SL.mat')

%% Create function handle for simulations
% Define parameters to estimate
parameters={'r' 'TV_max' 'e_Td' 'k_LN' 'K_TCD' 'K_pdl' 'S_R' 'k_pro' 'S_L'};

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
%       r      TV_max   e_Td    k_LN    K_TCD   K_pdl    S_R    k_pro   S_L                 b_0
mu=     [0.4   10       10      1e-2    0.2    .1       0.03    1       repelem(.0089,13)   repelem(1,6)];
% sigma_prior
sigma=  [2     3e0      3e0     3e0     3e0     3       3       3       repelem(3,13)       repelem(1.1,6)];
% Initial value
p0=     [0.1   10       3.5     0.01    0.2     0.02    0.03    1       repelem(.0089,13)   repelem(0.1,6)];
% Parameter boundaries
p_lb=   [0.01  1e0      1e-3    1e-4    1e-4    1e-4    1e-4    1e-3    repelem(1e-4,1,13)  repelem(1e-3,1,6)];
p_ub=   [10    1e4      1e3     1e0     1e1     1e1     1e1     4       repelem(1e1,1,13)   repelem(10,1,6)];

% Generating PI
PI.n_data=sum(cellfun(@(x)sum(sum(~isnan(x))),{PI.data.dataValue},'UniformOutput',true));
sigmaNames={'b_TV'; 'b_CD8'; 'b_CD107a'; 'b_DC'; 'b_MDSC'; 'b_PDL1'};
PI.par=getParamStruct(sim,1:8,9,13,[p_lb' p0' p_ub' mu' sigma'],sigmaNames);
PI.variableUnits={'Volume [mL]' 'Percentage [%]' 'Percentage [%]' 'Relative units []'...
    'Relative units []' 'Relative units []'};
% Hierarchical structure
H.PopulationParams=1:8;
H.IndividualParams=struct('name', {'S_L'}, 'Index', {9:21});
H.SigmaParams=22:27;

% Residuals function
residuals_fun=@(p)getResiduals(p,@(x)sim(x,100,u,1:1:100),PI,...
    @(x)getPhi2(x,H,size(PI.data,1)),(@(x)getCovariance(x,H)),4:6);

% Log-ikelihood function
likelihood_fun=@(p)sum(residuals_fun(exp(p))*(-1));
prior_fun=@(p)createPriorDistribution2(p,PI,H,'type','normal');

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
w0=[log(p0);randn(300,size(p0,2))*0.1+repmat(log(p0),300,1)];
logL=rowfun(obj_fun,table(w0));
[~,I]=sort(logL{:,:});

w0=w0(I(1:size(w0,2)*2),:);

% Ensemble MCMC sampler
tic
[models,logP]=gwmcmc(w0',{prior_fun likelihood_fun},1e5, 'StepSize', 2.5,...
'ThinChain',1,'BurnIn',0);
toc

% Plot entire chain
plotMCMCDiagnostics(exp(models_array),logL,PI)

% Plot thinned out samples
postSample=models_array(:,:,1.8e4:200:end);
postSample=postSample(:,:)';
plotMCMCDiagnostics(postSample,logL(:,:,1.8e4:2000:end),PI,'model', 'Kosisnky (S_L)')
figure
hold on
for i=1:13 
    histogram(exp(postSample(:,H.IndividualParams(1).Index(i))) )
end
    legend(PI.data(1:end).Group)

% Bivariate marginal plots
plotBivariateMarginals_2(postSample(:,H.PopulationParams),PI)
plotBivariateMarginals_2(postSample(:,H.IndividualParams(1).Index),PI,'names',...
    {PI.par(H.IndividualParams(1).Index).name})
plotBivariateMarginals_2(postSample(:,H.SigmaParams),PI,'names', {PI.par(H.SigmaParams).name})
%% MH MCMC sampler
tic
[params,logP, accept]=mcmc_mh(mean(postSample),likelihood_fun,prior_fun,1e4,...
    'ProposalSigma',cov(postSample),'BurnIn', 0,'Thinning', 1,'StepSize', 0.05);
toc

%% Posterior predictions
simFun=@(x)getOutput(PI,@(p)sim(p,100,u,1:1:100),x,...
    @(p)getPhi2(p,H,length(u)), 4:6);
tic
PI=getPosteriorPredictions(exp(postSample),PI,simFun,observables);
toc
PI=getCredibleIntervals(PI,observables, exp(postSample),H);
plotIndividualPosteriorPredictions(PI,observables)
plotGroupPosteriorPredictions(PI,observables(1),exp(postSample),H)

%%  Hamiltonian MCMC

logPosterior_fun=@(p)logPosterior(p,likelihood_fun,prior_fun,PI,H);
smp = hmcSampler(logPosterior_fun,mean(postSample),'NumSteps',50);
[MAPpars,fitInfo] = estimateMAP(smp,'VerbosityLevel',1);
% Tune sampler
tic
[smp,tuneinfo] = tuneSampler(smp,'Start',MAPpars,'Verbosity', 5);
toc
figure;
plot(tuneinfo.StepSizeTuningInfo.StepSizeProfile);
xlabel('Iteration');
ylabel('Step size');
set(gca,'YScale','log')
accratio = tuneinfo.StepSizeTuningInfo.AcceptanceRatio
% Draw samples using 4 chains
NumChains = 1;
chains = cell(NumChains,1);
Burnin = 500;
NumSamples = 100;
for c = 1:NumChains
    if (c == 1)
        level = 1;
    else
        level = 0;
    end
    chains{c} = drawSamples(smp,'Start',MAPpars + randn(size(MAPpars)), ...
        'Burnin',Burnin,'NumSamples',NumSamples, ...
        'VerbosityLevel',5,'NumPrint',300);
end
%% Save results
save('PI_SL.mat', 'PI')
save('models_Kosinsky_MOC1_SL.mat','models_array')
save('logP_Kosinsky_MOC1_SL.mat','logL')
