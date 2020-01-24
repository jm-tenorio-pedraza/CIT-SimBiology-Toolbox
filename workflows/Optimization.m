
%% Optimization

% Optimizer options
options_fminsearch=optimset('Display','iter','MaxFunEvals', 1e4, 'MaxIter',1e4, 'TolFun', 1e-4);
options_anneal.Verbosity=2;
options_anneal.InitTemp=100;
options = optimoptions('simulannealbnd','PlotFcns',...
          {@saplotbestx,@saplotbestf,@saplotx,@saplotf},'InitialTemperature', 100, 'MaxFunctionEvaluations', 1e4);
delta = 1;

%% Local optimisation

ub = log([PI.par([PI.H.PopulationParams PI.H.CellParams.Index PI.H.IndividualParams.Index]).maxValue]);
lb = log([PI.par([PI.H.PopulationParams PI.H.CellParams.Index PI.H.IndividualParams.Index]).minValue]);
[phat, ~] = lsqnonlin(residuals_fn,(finalValues([PI.H.PopulationParams...
    PI.H.CellParams.Index PI.H.IndividualParams.Index])), lb,ub, options_fminsearch);

finalValues([PI.H.PopulationParams PI.H.CellParams.Index PI.H.IndividualParams.Index])=phat;
 finalValues([PI.H.IndividualParams.OmegaIndex]) = arrayfun(@(x)log(std(finalValues(x.Index))), PI.H.IndividualParams);
% finalValues([PI.H.CellParams.OmegaIndex]) = arrayfun(@(x)log(std(finalValues(x.Index))), PI.H.CellParams);

%% Global optimisation
while delta >1e-4
% Nelder-Mead
[finalValues, fval_fminsearch]=fminsearch(obj_fun,p_hat,options_fminsearch);

% Simulated annealing
[p_hat, fval_anneal]=anneal(obj_fun,finalValues,options_anneal);

% Local optimization
% Residuals 
residuals_fn = @(x) getResiduals(exp(x),@(x)sim(x,PI.tspan(end),PI.u,PI.tspan),PI,...
    @(x)getPhi2(x,PI.H,length(PI.u),'initialvalue',PI.x_0),exp(p_hat(setdiff([PI.H.SigmaParams],...
    [PI.H.CellParams(:).OmegaIndex PI.H.IndividualParams(:).OmegaIndex]))),PI.normIndx);

[phat, ~] = lsqnonlin(residuals_fn,(p_hat([PI.H.PopulationParams...
    PI.H.CellParams.Index PI.H.IndividualParams.Index])), lb,ub, options_fminsearch);
 finalValues=p_hat;

finalValues([PI.H.IndividualParams.OmegaIndex]) = arrayfun(@(x)std(finalValues(x.Index)), PI.H.IndividualParams);
finalValues([PI.H.CellParams.OmegaIndex]) = arrayfun(@(x)std(finalValues(x.Index)), PI.H.CellParams);

% % [finalValues, fval_fminunc] = fminunc(obj_fun,finalValues,options_fminsearch);
% fval_fminunc = obj_fun(finalValues);
delta = abs(fval_anneal - fval_fminsearch);
end

 %% SAEM                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          %% SAEM
residuals_fx = @(x,sigma) getResiduals(exp(x),@(x)sim(x,PI.tspan(end),PI.u,PI.tspan),PI,...
    @(x)getPhi2(x,PI.H,length(PI.u),'initialValue',PI.x_0),...
    sigma,PI.normIndx);
tic
[finalValues, logL] = saem(finalValues, residuals_fx,likelihood_fun, prior_fun, PI.H, PI,...
    'm', 1e4,'StepSize',2.38^2,'MinFunc', 'fminsearch','OutputFn',...
    @(x)getOutput(PI,@(p)sim(p,PI.tspan(end),PI.u,PI.tspan),exp(x),...
    @(p)getPhi2(p,PI.H,length(PI.u),'initialValue',PI.x_0), PI.normIndx,PI.H));
toc1
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
