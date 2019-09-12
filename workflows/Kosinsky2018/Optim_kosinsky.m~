%% Optimization of Kosinsky et al's model to find initial parameter estimate

p0 = log([PI.par(:).startValue]);
% Optimizer options
options_fminsearch=optimset('Display','iter','MaxFunEvals', 1e4, 'MaxIter',1e4, 'TolFun', 1e-6);
options_anneal.Verbosity=2;
options_anneal.InitTemp=100;
options = optimoptions('simulannealbnd','PlotFcns',...
          {@saplotbestx,@saplotbestf,@saplotx,@saplotf},'InitialTemperature', 100, 'MaxFunctionEvaluations', 1e4);

% Nelder-Mead
p_hat=fminsearch(obj_fun,finalValues,options_fminsearch);

% Simulated annealing
finalValues=anneal(obj_fun,p_hat,options_anneal);
tic
finalValues=simulannealbnd(obj_fun,finalValues,[],[],options);
toc

% Stochastic EM
[params, logL] = saem(p_hat, likelihood_fun, prior_fun,H,'m', 1e3);

% Simulation output
PI=getOutput(PI,@(p)sim(p,100,u,1:1:100),exp(finalValues),...
    @(p)getPhi2(p,H,length(u)), length(observables)-1:length(observables),1:100);

% Plotting tumor volume
for i=1:length(observables)
plotSimOutput(PI.data,i)
legend(observables(i))
end

finalValue=num2cell(exp(finalValues'));
[PI.par(1:end).finalValue]=finalValue{:,:};
%% Save results
save('PI_Kosinsky.mat', 'PI')
