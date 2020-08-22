% Local Optimization
popParamsIndx = [PI.H.PopulationParams PI.H.CellParams.Index];
indivParamsIndx = [PI.H.IndividualParams.Index];
options_fminsearch=optimset('Display','iter','MaxFunEvals', 5e4, 'MaxIter',5e4, 'TolFun', 1e-4);
ub = log([PI.par([popParamsIndx indivParamsIndx]).maxValue]);
lb = log([PI.par([popParamsIndx indivParamsIndx]).minValue]);

% Residuals 
residuals_fn = @(x) getResiduals(exp(x),@(x)sim(x,PI.tspan(end),PI.u,PI.tspan),PI,...
    @(x)getPhi2(x,PI.H,length(PI.u),'initialValue',PI.x_0),...
    exp(finalValues(end-length(observables)+1:end)),exp(finalValues([PI.H.CellParams.OmegaIndex])),...
    exp(finalValues([PI.H.IndividualParams.OmegaIndex])),PI.normIndx);


%% Estimate Individual params
residuals_indiv = @(x) residuals_fn([finalValues(popParamsIndx) x]);

[p_hat_indiv, ~] = lsqnonlin(residuals_indiv,finalValues(indivParamsIndx),...
    lb(indivParamsIndx),...
    ub(indivParamsIndx), options_fminsearch);
finalValues(indivParamsIndx) = p_hat_indiv;

%% Esitmate Population params
residuals_pop = @(x) residuals_fn([x finalValues(indivParamsIndx)]);

[p_hat_pop, ~] = lsqnonlin(residuals_pop,finalValues(popParamsIndx),...
    lb(popParamsIndx),...
    ub(popParamsIndx), options_fminsearch);
finalValues(popParamsIndx) = p_hat_pop;
%% Estimate both
[p_hat, ~] = lsqnonlin(residuals_fn,(finalValues([popParamsIndx indivParamsIndx])),...
    lb,ub, options_fminsearch);

finalValues([ popParamsIndx indivParamsIndx ]) = p_hat;
%% Estimate all parameters
[finalValues, fval_fminunc] = fminunc(obj_fun, finalValues,options_fminsearch);
                                                                                           %% 
clearvars residuals_fn residuals_indiv 