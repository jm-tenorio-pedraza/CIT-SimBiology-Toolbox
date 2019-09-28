
%% Objective function
finalValues = log([PI.par(:).startValue]);

% Obj function
obj_fun=@(x)(likelihood_fun(x)*(-1)+prior_fun(x')*(-1));
tic
obj_fun(finalValues)
toc

%% Optimization

% Optimizer options
options_fminsearch=optimset('Display','iter','MaxFunEvals', 1e4, 'MaxIter',1e4, 'TolFun', 1e-4);
options_anneal.Verbosity=2;
options_anneal.InitTemp=100;
options = optimoptions('simulannealbnd','PlotFcns',...
          {@saplotbestx,@saplotbestf,@saplotx,@saplotf},'InitialTemperature', 100, 'MaxFunctionEvaluations', 1e4);
delta = 1;

%% Local optimisation

ub = log([PI.par([H.PopulationParams H.IndividualParams.Index]).maxValue]);
lb = log([PI.par([H.PopulationParams H.IndividualParams.Index]).minValue]);

[p_hat, ~] = lsqnonlin(residuals_fn,finalValues([H.PopulationParams H.IndividualParams.Index]), lb,ub, options_fminsearch);
finalValues([H.PopulationParams H.IndividualParams.Index])=p_hat;
%% Global optimisation
while delta >1e-2
% Nelder-Mead
[p_hat, fval_fminsearch]=fminsearch(obj_fun,finalValues,options_fminsearch);
% Simulated annealing
[finalValues, fval_anneal]=anneal(obj_fun,p_hat,options_anneal);
delta = abs(fval_fminsearch - fval_anneal);
end
tic
finalValues=simulannealbnd(obj_fun,finalValues,[],[],options);
toc
% Stochastcic EM
[params, logL] = saem(p_hat', likelihood_fun, prior_fun,H,'m', 2e3,...
    'StepSize',0.05,'MinFunc', 'fminsearch','OutputFn', @(x)getOutput(PI,@(p)sim(p,100,u,1:100),exp(x),...
    @(p)getPhi2(p,H,length(u)),7:8,1:100));

%% Simulation output
PI=getOutput(PI,@(p)sim(p,100,u,1:1:100),exp(p_hat),...
    @(p)getPhi2(p,H,length(u)), [],1:100);
 
      
% Plotting tumor volume
for i=1:length(observables)
plotSimOutput(PI.data,i)
legend(observables(i))

end

finalValue=num2cell(exp(p_hat'));
[PI.par(1:end).finalValue]=finalValue{:,:};

%% Create variant object
MOC1_optimized = createVariant(PI,H,'MOC1_optimized');
variant = addvariant(model, 'MOC1_optimized');

addcontent(variant, MOC1_optimized.Content);
out.m1 = model;
sbiosaveproject '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/CIM_2.sbproj' 'model'