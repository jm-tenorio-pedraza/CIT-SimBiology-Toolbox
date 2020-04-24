% Local Optimization
options_fminsearch=optimset('Display','iter','MaxFunEvals', 5e4, 'MaxIter',5e4, 'TolFun', 1e-4);
ub = log([PI.par([PI.H.PopulationParams PI.H.CellParams.Index PI.H.IndividualParams.Index]).maxValue]);
lb = log([PI.par([PI.H.PopulationParams PI.H.CellParams.Index PI.H.IndividualParams.Index]).minValue]);

% Residuals 
residuals_fn = @(x) getResiduals(exp(x),@(x)sim(x,PI.tspan(end),PI.u,PI.tspan),PI,...
    @(x)getPhi2(x,PI.H,length(PI.u),'initialValue',PI.x_0),...
    exp(finalValues(end-length(observables)+1:end)),exp(finalValues([PI.H.CellParams.OmegaIndex])),...
    exp(finalValues([PI.H.IndividualParams.OmegaIndex])),PI.normIndx);


% Esitmate Population params
residuals_pop = @(x) residuals_fn([x finalValues([PI.H.CellParams.Index PI.H.IndividualParams.Index])]);

[p_hat_pop, ~] = lsqnonlin(residuals_pop,finalValues([PI.H.PopulationParams]),...
    log([PI.par([PI.H.PopulationParams]).minValue]),...
    log([PI.par([PI.H.PopulationParams]).maxValue]), options_fminsearch);
finalValues([PI.H.PopulationParams]) = p_hat_pop;
% Estimate Individual params
residuals_indiv = @(x) residuals_fn([finalValues([PI.H.PopulationParams]) x]);

[p_hat_indiv, ~] = lsqnonlin(residuals_indiv,(finalValues([PI.H.CellParams.Index...
    PI.H.IndividualParams.Index])), log([PI.par([PI.H.CellParams.Index...
    PI.H.IndividualParams.Index]).minValue]),log([PI.par([PI.H.CellParams.Index...
    PI.H.IndividualParams.Index]).maxValue]), options_fminsearch);
finalValues([PI.H.CellParams.Index...
    PI.H.IndividualParams.Index]) = p_hat_indiv;


% Estimate both
[p_hat, ~] = lsqnonlin(residuals_fn,(finalValues([PI.H.PopulationParams...
    PI.H.CellParams.Index PI.H.IndividualParams.Index])), lb,ub, options_fminsearch);

finalValues([PI.H.PopulationParams PI.H.CellParams.Index...
    PI.H.IndividualParams.Index]) = p_hat;
%% Estimate all parameters
[finalValues, fval_fminunc] = fminunc(obj_fun, finalValues,options_fminsearch);