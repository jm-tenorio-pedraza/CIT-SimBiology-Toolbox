%% Optimization
options_fminsearch=optimset('Display','iter','MaxFunEvals', 2e4, 'MaxIter',...
    5e4, 'TolFun', 1e-4);
options_anneal.Verbosity=2;
options_anneal.InitTemp=100;
%% Joint optimization
[finalValues,fval_anneal]=anneal(obj_fun,finalValues,options_anneal);
[finalValues,fval_fminsearch]=fminsearch(obj_fun,finalValues,options_fminsearch);
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
 subplot(nrow,ncol,i)
 plotSimOutput(PI,i,'all', false, 'indiv', false, 'addErrorVar', false,...
     'newFig', false, 'TimeUnit', 'days')
%  set(gca, 'YScale','log')
end

%%
finalValue=num2cell(exp(finalValues'));
[PI.par(1:end).finalValue]=finalValue{:,:};
% plotFit(PI)

%% Plotting errors
figure('Position', [10 10 1.5e3 1e3])
ncol = ceil(sqrt(length(observables)));
nrow = ceil(length(observables)/ncol);
for i=1:length(observables)
 subplot(nrow,ncol,i)
 plotError2(PI,i,'all', false, 'indiv', false, 'addErrorVar', false, 'newFig', ...
     false, 'group', 'Cell', 'TimeUnit', 'hours')
end
