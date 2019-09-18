%% Optimization of Kosinsky et al's model to find initial parameter estimate

finalValues = log([PI.par(:).startValue]);
% Optimizer options
options_fminsearch=optimset('Display','iter','MaxFunEvals', 1e4, 'MaxIter',1e4, 'TolFun', 1e-4);
options_anneal.Verbosity=2;
options_anneal.InitTemp=100;
options = optimoptions('simulannealbnd','PlotFcns',...
          {@saplotbestx,@saplotbestf,@saplotx,@saplotf},'InitialTemperature', 100, 'MaxFunctionEvaluations', 1e4);
delta = 1;
while delta >1e-4
% Nelder-Mead
[p_hat, fval_fminsearch]=fminsearch(obj_fun,finalValues,options_fminsearch);

% Simulated annealing
[finalValues, fval_anneal]=anneal(obj_fun,params,options_anneal);
delta = abs(fval_fminsearch - fval_anneal);
end
tic
finalValues=simulannealbnd(obj_fun,p_hat,[],[],options);
toc
% Stochastic EM
[params, logL] = saem(finalValues, likelihood_fun, prior_fun,H,'m', 1e3,...
    'StepSize',0.00005,'MinFunc', 'fminsearch','OutputFn', @(x)getOutput(PI,@(p)sim(p,100,u,1:100),exp(x),...
    @(p)getPhi2(p,H,length(u)),7:8,1:100));

%% Simulation output
PI=getOutput(PI,@(p)sim(p,100,u,1:1:100),exp(finalValues),...
    @(p)getPhi2(p,H,length(u)), length(observables)-1:length(observables),1:100);
PI = getOutput(PI,@(p)sim(p,100,u,1:100),exp(finalValues),...
    @(p)getPhi2(p,H,length(u)),7:8,1:100);
 
E_likelihood = (@(x)sum(getErrors(PI,exp(x)))*(-1));
E_prior = @(x) prior_fun([finalValues([H.PopulationParams]) finalValues(H.IndividualParams.Index) x]);

[sigma, logP_sigma, ~] = mcmc_mh(finalValues(setdiff(H.SigmaParams, H.IndividualParams.OmegaIndex)),E_likelihood, E_prior, 1000);
      
% Plotting tumor volume
for i=1:length(observables)
plotSimOutput(PI.data,i)
legend(observables(i))
end

finalValue=num2cell(exp(params'));
[PI.par(1:end).finalValue]=finalValue{:,:};
%% Save results
save('PI.mat', 'PI')

load('PI_CIM.mat','PI')
%% Create variant object
MOC1_optimized = createVariant(PI,H,'MOC1_optimized');
variant = addvariant(model, 'MOC1_optimized');

addcontent(variant, MOC1_optimized.Content);
out.m1 = model;
sbiosaveproject '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/CIM_2.sbproj' 'model'