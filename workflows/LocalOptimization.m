delta = 1;
popParamsIndx = [PI.H.PopulationParams];
cellParamsIndx = [PI.H.CellParams.Index];
indivParamsIndx = [PI.H.IndividualParams.Index];
respParamsIndx = [PI.H.RespParams.Index];

options_lsqnonlin=optimoptions('lsqnonlin','Display','iter','MaxFunEvals',...
    5e4, 'MaxIter',5e4, 'TolFun', 1e-4, 'StepTolerance', 1e-6);
ub = log([PI.par([popParamsIndx cellParamsIndx indivParamsIndx respParamsIndx]).maxValue]);
lb = log([PI.par([popParamsIndx cellParamsIndx indivParamsIndx respParamsIndx]).minValue]);

while delta>1e-3
phat_0 = finalValues;
% Local Optimization

% Residuals 
residuals_fn = @(x) getResiduals(exp(x),@(x)sim(x,PI.tspan(end),PI.u,PI.tspan),PI,...
    @(x)getPhi2(x,PI.H,length(PI.u),'initialValue',PI.x_0),...
    exp(phat_0(end-length(observables)+1:end)),...
    exp(phat_0([PI.H.CellParams.OmegaIndex])),...
    exp(phat_0([PI.H.IndividualParams.OmegaIndex])),...
    exp(phat_0([PI.H.RespParams.OmegaIndex])),PI.normIndx,...
    'log', false, 'errorModel', 'multiplicative');


%% Estimate Individual params
try
residuals_indiv = @(x) residuals_fn([finalValues(popParamsIndx) finalValues(cellParamsIndx) x  finalValues(respParamsIndx)]);

[p_hat_indiv, ~] = lsqnonlin(residuals_indiv,finalValues(indivParamsIndx),...
    lb(indivParamsIndx),...
    ub(indivParamsIndx), options_lsqnonlin);
finalValues(indivParamsIndx) = p_hat_indiv;
catch
end
%% Estimate cell params
try
residuals_cell = @(x) residuals_fn([finalValues(popParamsIndx) x finalValues(indivParamsIndx)  finalValues(respParamsIndx)]);

[p_hat_cell, ~] = lsqnonlin(residuals_cell,finalValues(cellParamsIndx),...
    lb(cellParamsIndx),...
    ub(cellParamsIndx), options_lsqnonlin);
finalValues(cellParamsIndx) = p_hat_cell;
catch
end
%% Esitmate Population params

residuals_pop = @(x) residuals_fn([x finalValues(cellParamsIndx) finalValues(indivParamsIndx)  finalValues(respParamsIndx)]);

[p_hat_pop, ~] = lsqnonlin(residuals_pop,finalValues(popParamsIndx),...
    lb(popParamsIndx),...
    ub(popParamsIndx), options_lsqnonlin);
finalValues(popParamsIndx) = p_hat_pop;
%% Estimate pop, cell, indiv and resp parameters
[p_hat, ~] = lsqnonlin(residuals_fn,(finalValues([popParamsIndx cellParamsIndx indivParamsIndx respParamsIndx])),...
    lb,ub, options_lsqnonlin);

finalValues([ popParamsIndx cellParamsIndx indivParamsIndx respParamsIndx]) = p_hat;
%% Estimate all parameters
options_fminunc = optimoptions('fminunc','Display','iter','MaxFunEvals',...
    5e4, 'MaxIter',5e4, 'TolFun', 1e-4, 'StepTolerance', 1e-8);
[finalValues, fval_fminunc] = fminunc(obj_fun, finalValues,options_fminunc);
delta = obj_fun(phat_0) - fval_fminunc;
end                                                                             %% 
clearvars residuals_fn residuals_indiv 