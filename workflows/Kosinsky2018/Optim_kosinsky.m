%% Optimization of Kosinsky et al's model to find initial parameter estimate


% Optimizer options
options_fminsearch=optimset('Display','iter','MaxFunEvals', 1e4, 'MaxIter',1e4);
options_anneal.Verbosity=2;
options_anneal.InitTemp=2;
options = optimoptions('simulannealbnd','PlotFcns',...
          {@saplotbestx,@saplotbestf,@saplotx,@saplotf},'InitialTemperature', 10, 'MaxFunctionEvaluations', 1e4);

% Nelder-Mead
p_hat=fminsearch(obj_fun,log(p0),options_fminsearch);

% Simulated annealing
finalValues=anneal(obj_fun,p_hat,options_anneal);
tic
finalValues=simulannealbnd(obj_fun,finalValues,[],[],options);
toc
% Simulation output
PI=getOutput(PI,@(p)sim(p,100,u,1:1:100),exp(finalValues),...
    @(p)getPhi2(p,H,length(u)), 4:6);

% Plotting tumor volume
for i=1:6
plotSimOutput(PI.data,i)
legend(observables(i))
end

finalValue=num2cell(exp(finalValues'));
[PI.par(1:end).finalValue]=finalValue{:,:};
%% Save results
save('PI_kpro.mat', 'PI')
