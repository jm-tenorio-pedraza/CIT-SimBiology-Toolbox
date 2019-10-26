
%% Objective function

% Obj function
obj_fun=@(x)(likelihood_fun(x)*(-1)+prior_fun(x)*(-1));
tic
obj_fun((finalValues))
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

ub = log([PI.par([H.PopulationParams H.CellParams.Index H.IndividualParams.Index]).maxValue]);
lb = log([PI.par([H.PopulationParams H.CellParams.Index H.IndividualParams.Index]).minValue]);
[phat, ~] = lsqnonlin(residuals_fn,(finalValues([H.PopulationParams...
    H.CellParams.Index H.IndividualParams.Index])), lb,ub, options_fminsearch);

finalValues([H.PopulationParams H.CellParams.Index H.IndividualParams.Index])=phat;

%% Global optimisation
while delta >1e-4
% Nelder-Mead
  [finalValues, fval_fminsearch]=fminsearch(obj_fun,finalValues,options_fminsearch);

% Simulated annealing
[p_hat, fval_anneal]=anneal(obj_fun,finalValues,options_anneal);

% Local optimization
% Residuals 
residuals_fn = @(x) getResiduals(exp(x),@(x)sim(x,144,u,1:.1:144),PI,...
    @(x)getPhi2(x,H,length(u),'initialvalue',x_0),exp(p_hat(setdiff([H.SigmaParams],...
    [H.CellParams(:).OmegaIndex H.IndividualParams(:).OmegaIndex]))),normIndx);
[phat, ~] = lsqnonlin(residuals_fn,(p_hat([H.PopulationParams...
    H.CellParams.Index H.IndividualParams.Index])), lb,ub, options_fminsearch);
finalValues=p_hat;

finalValues([H.PopulationParams H.CellParams.Index H.IndividualParams.Index])=phat;
finalValues([H.CellParams.OmegaIndex]) = arrayfun(@(x)std(finalValues(x.Index)), H.CellParams);
finalValues([H.IndividualParams.OmegaIndex]) = arrayfun(@(x)std(finalValues(x.Index)), H.IndividualParams);

% [finalValues, fval_fminunc,~,~,grad,hessian] = fminunc(obj_fun,finalValues,options_fminsearch);
fval_fminunc = obj_fun(finalValues);
delta = abs(fval_anneal - fval_fminunc);
end

%% SAEM
[params, logL] = saem(finalValues, residuals_func, prior_fun,H,PI,...
    'm', 5e3,'StepSize',0.1,'MinFunc', 'lsqnonlin','OutputFn',...
    @(x)getOutput(PI,@(p)sim(p,100,u,1:100),exp(x),...
    @(p)getPhi2(p,H,length(u),x_0),normIndx,1:100));

%% Simulation output
PI=getOutput(PI,@(p)sim(p,100,u,1:1:PI.tspan(end)),exp(finalValues),...
    @(p)getPhi2(p,H,length(u),'initialValue',x_0), normIndx,H);
 
     
%% Plotting output
for i=1:length(observables)
plotSimOutput(PI,i)
legend(observables(i),'Location', 'best')
end
%%
finalValue=num2cell(exp(finalValues'));
[PI.par(1:end).finalValue]=finalValue{:,:};

%% Plotting errors
for i=1:length(observables)
plotError(exp(finalValues(setdiff(H.SigmaParams,...
    [H.IndividualParams.OmegaIndex, H.CellParams.OmegaIndex]))),PI,i)
legend(observables(i),'Location', 'best')
end

%% Create variant object
MOC1_optimized = createVariant(PI,H,'MOC1_optimized');
variant = addvariant(model, 'MOC1_optimized');
addcontent(variant, MOC1_optimized.Content);
out.m1 = model;
