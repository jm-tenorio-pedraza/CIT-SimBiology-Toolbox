
%% Optimization

% Optimizer options
options_fminsearch=optimset('Display','iter','MaxFunEvals', 1e4, 'MaxIter',1e4, 'TolFun', 1e-4);
options_anneal.Verbosity=2;
options_anneal.InitTemp=100;
options = optimoptions('simulannealbnd','PlotFcns',...
          {@saplotbestx,@saplotbestf,@saplotx,@saplotf},'InitialTemperature', 100, 'MaxFunctionEvaluations', 1e4);
delta = 1;


%% Global optimisation

while delta >1e-4

% Population parameter optimization
obj_fun_pop = @(x)obj_fun([x finalValues([PI.H.CellParams.Index PI.H.IndividualParams.Index...
    PI.H.SigmaParams])]);
% Simulated annealing
[p_hat_pop, ~]=anneal(obj_fun_pop,finalValues([PI.H.PopulationParams]),options_anneal);
% Nelder-Mead
[p_hat_pop, ~]=fminsearch(obj_fun_pop,p_hat_pop,options_fminsearch);

finalValues([PI.H.PopulationParams]) = p_hat_pop;

% Individual Parameter optimization
obj_fun_indiv = @(x)obj_fun([finalValues([PI.H.PopulationParams]) x finalValues([PI.H.SigmaParams])]);
% Simulated annealing
[p_hat_indiv, fval_anneal]=anneal(obj_fun_indiv,finalValues([PI.H.CellParams.Index...
    PI.H.IndividualParams.Index]),options_anneal);
% Nelder-Mead
[p_hat_indiv, fval_fminsearch]=fminsearch(obj_fun_indiv,p_hat_indiv,options_fminsearch);

finalValues([PI.H.CellParams.Index PI.H.IndividualParams.Index]) = p_hat_indiv;

[finalValues, fval_anneal]=anneal(obj_fun, finalValues, options_anneal);
delta = abs(fval_anneal - fval_fminsearch);
end

 %% SAEM                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          %% SAEM
residuals_fx = @(x,sigma) getResiduals(exp(x),@(x)sim(x,PI.tspan(end),PI.u,PI.tspan),PI,...
    @(x)getPhi2(x,PI.H,length(PI.u),'initialValue',PI.x_0),...
    sigma,PI.normIndx);
tic
[finalValues, logL] = saem(finalValues, residuals_fx,likelihood_fun, prior_fun, PI.H, PI,...
    'm', 2e3,'StepSize',2.38^2,'MinFunc', 'fminsearch','OutputFn',...
    @(x)getOutput(PI,@(p)sim(p,PI.tspan(end),PI.u,PI.tspan),exp(x),...
    @(p)getPhi2(p,PI.H,length(PI.u),'initialValue',PI.x_0), PI.normIndx,PI.H));
toc
PI.AIC = 2*length(PI.par)-2*obj_fun(finalValues)*(-1);

%% Simulation output
PI=getOutput(PI,@(p)sim(p,PI.tspan(end),PI.u,PI.tspan),exp(finalValues),...
    @(p)getPhi2(p,PI.H,length(PI.u),'initialValue',PI.x_0), PI.normIndx,PI.H);
 
     
%% Plotting output
for i=1:length(observables)
plotSimOutput(PI,i)
end
%%
finalValue=num2cell(exp(finalValues'));
[PI.par(1:end).finalValue]=finalValue{:,:};
plotFit(PI,'sigma', exp(finalValues(setdiff(PI.H.SigmaParams,[PI.H.IndividualParams.OmegaIndex PI.H.CellParams.OmegaIndex]))))
%% Plotting errors
for i=1:length(observables)
plotError(exp(finalValues(setdiff(PI.H.SigmaParams,...
    [PI.H.IndividualParams.OmegaIndex, PI.H.CellParams.OmegaIndex]))),PI,i)
legend(observables(i),'Location', 'best')
end

%% Create variant object
MOC1_optimized = createVariant(PI,H,'MOC1_optimized');
variant = addvariant(model, 'MOC1_optimized');
addcontent(variant, MOC1_optimized.Content);
out.m1 = model;
