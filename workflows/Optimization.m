%% Optimization
options_fminsearch=optimset('Display','iter','MaxFunEvals', 5e4, 'MaxIter',...
    5e4, 'TolFun', 1e-4);
options_anneal.Verbosity=2;
options_anneal.InitTemp=100;
%% Partitioned optimisation
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
delta = abs(fval_anneal - fval_fminsearch);
end
%% Joint optimization
 [finalValues, fval_anneal]=anneal(obj_fun,finalValues,options_anneal);
 [finalValues, fval_fminsearch]=fminsearch(obj_fun,finalValues,options_fminsearch);
%% Simulation output
simTime = unique([PI.tspan', 1:PI.tspan(end)]);
PI=getOutput(PI,@(p)sim(p,PI.tspan(end),PI.u, simTime),exp(finalValues),...
    @(p)getPhi2(p,PI.H,length(PI.u),'initialValue',PI.x_0),...
    PI.normIndx,PI.H,'output', 'PI', 'simTime', simTime);
PI.AIC = 2*length(PI.par)-2*likelihood_fun(finalValues)*(1);
%% Plotting output
figure('Position', [10 10 1.5e3 1e3])
ncol = ceil(sqrt(length(observables)));
nrow = ceil(length(observables)/ncol);
for i=1:length(observables)
 %subplot(nrow,ncol,i)
 plotSimOutput(PI,i,'all', false, 'indiv', true, 'addErrorVar', false,...
     'newFig', true, 'TimeUnit', 'days')
%  set(gca, 'YScale','log')
end
%%
finalValue=num2cell(exp(finalValues'));
[PI.par(1:end).finalValue]=finalValue{:,:};
plotFit(PI)
%%                                                                                                                                                                 %% Plotting errors
figure('Position', [10 10 1.5e3 1e3])
ncol = ceil(sqrt(length(observables)));
nrow = ceil(length(observables)/ncol);
for i=1:length(observables)
 subplot(nrow,ncol,i)
 plotError2(PI,i,'all', false, 'indiv', false, 'addErrorVar', false, 'newFig', ...
     false, 'group', 'Cell', 'TimeUnit', 'hours')
end
