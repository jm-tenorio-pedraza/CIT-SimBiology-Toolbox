
%% Objective function

% Obj function
obj_fun=@(x)(likelihood_fun(x)*(-1)+prior_fun(x')*(-1));
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

ub = log([PI.par([H.PopulationParams H.IndividualParams.Index]).maxValue]);
lb = log([PI.par([H.PopulationParams H.IndividualParams.Index]).minValue]);
[phat, ~] = lsqnonlin(residuals_fn,(finalValues([H.PopulationParams...
    H.IndividualParams.Index])), lb,ub, options_fminsearch);

finalValues([H.PopulationParams H.IndividualParams.Index])=phat;

%% Global optimisation
while delta >1e-2
% Nelder-Mead
[finalValues, fval_fminsearch]=fminsearch(obj_fun,finalValues,options_fminsearch);

% Simulated annealing
[p_hat, fval_anneal]=anneal(obj_fun,finalValues,options_anneal);

% Local optimization
% Residuals 
residuals_fn = @(x) getResiduals(exp(x),@(x)sim(x,144,u,1:.1:144),PI,...
    @(x)getPhi2(x,H,length(u),'initialvalue',x_0),exp(p_hat(setdiff([H.SigmaParams],...
    [H.IndividualParams(:).OmegaIndex]))),normIndx);
[phat, ~] = lsqnonlin(residuals_fn,(p_hat([H.PopulationParams...
    H.IndividualParams.Index])), lb,ub, options_fminsearch);
finalValues=p_hat;

finalValues([H.PopulationParams H.IndividualParams.Index])=phat;

[finalValues, fval_fminunc,~,~,grad,hessian] = fminunc(obj_fun,finalValues,options_fminsearch);
delta = abs(fval_fminsearch - fval_fminunc);
end
tic
finalValues=simulannealbnd(obj_fun,finalValues,[],[],options);
toc
residuals_func = @(x,sigma, addSigma) getResiduals(exp(x),@(x)sim(x,100,u,1:1:100),PI,...
    @(x)getPhi2(x,H,length(u),x_0),sigma,normIndx, 'addSigma', addSigma);

% Stochastcic EM
[params, logL] = saem(finalValues, residuals_func, prior_fun,H,PI,...
    'm', 5e3,'StepSize',0.1,'MinFunc', 'lsqnonlin','OutputFn',...
    @(x)getOutput(PI,@(p)sim(p,100,u,1:100),exp(x),...
    @(p)getPhi2(p,H,length(u),x_0),normIndx,1:100));

%% Simulation output
PI=getOutput(PI,@(p)sim(p,150,u,1:.1:150),exp(finalValues),...
    @(p)getPhi2(p,H,length(u),'initialValue',x_0), normIndx,H);
 
      
% Plotting tumor volume
for i=1:length(observables)
plotSimOutput(PI,i)
legend(observables(i))
end

finalValue=num2cell(exp(finalValues'));
[PI.par(1:end).finalValue]=finalValue{:,:};

%% Create variant object
MOC1_optimized = createVariant(PI,H,'MOC1_optimized');
variant = addvariant(model, 'MOC1_optimized');
addcontent(variant, MOC1_optimized.Content);
out.m1 = model;
