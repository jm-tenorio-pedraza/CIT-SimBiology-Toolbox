
%% Objective function
finalValues = log([PI.par(:).startValue]);

% Obj function
obj_fun=@(x)(likelihood_fun(x)*(-1)+prior_fun(x')*(-1));
tic
obj_fun(finalValues)
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
sigma_prior = @(x) prior_fun([finalValues([H.PopulationParams]) finalValues([H.IndividualParams.Index]) x]);
sigma_indx = size(H.IndividualParams.OmegaIndex,1)+1:length(H.SigmaParams);
sigma_likelihood = (@(x)sum(getErrors(PI,exp(x(sigma_indx))))*(-1));
sigma_obj = @(x)((sigma_prior(x)+sigma_likelihood(x))*(-1));

ub = log([PI.par([H.PopulationParams H.IndividualParams.Index]).maxValue]);
lb = log([PI.par([H.PopulationParams H.IndividualParams.Index]).minValue]);

[p_hat, ~] = lsqnonlin(residuals_fn,finalValues([H.PopulationParams H.IndividualParams.Index]), lb,ub, options_fminsearch);

lhs = lhsdesign(100,length(ub));
p0 = unifinv(lhs,repelem(lb,100,1),repelem(ub,100,1));
p_hat = nan(size(p0)); fval_lsqnonlin = nan(size(p0,1),1);
for i = 1:size(p0,1)
[p_hat(i,:), fval_lsqnonlin(i)] = lsqnonlin(loglikelihood_fn,p0(i,:), lb,ub, options_fminsearch);
end
[fval_sorted,I]=sort(fval_lsqnonlin);

p_hat_sorted = p_hat(I,:);
finalValues=[p_hat_sorted(1,:) finalValues([H.SigmaParams])];
set(gca,'YScale','log')
%% Global optimisation
while delta >1e-4
% Nelder-Mead
[p_hat, fval_fminsearch]=fminsearch(obj_fun,p_hat,options_fminsearch);
% Simulated annealing
[finalValues, fval_anneal]=anneal(obj_fun,p_hat,options_anneal);
delta = abs(fval_fminsearch - fval_anneal);
end
tic
finalValues=simulannealbnd(obj_fun,finalValues,[],[],options);
toc
% Stochastcic EM
[params, logL] = saem(p_hat', likelihood_fun, prior_fun,H,'m', 2e3,...
    'StepSize',0.05,'MinFunc', 'fminsearch','OutputFn', @(x)getOutput(PI,@(p)sim(p,100,u,1:100),exp(x),...
    @(p)getPhi2(p,H,length(u)),7:8,1:100));

%% Simulation output
PI=getOutput(PI,@(p)sim(p,100,u,1:1:100),exp(p_hat),...
    @(p)getPhi2(p,H,length(u)), [],1:100);
 
      
% Plotting tumor volume
for i=1:length(observables)
plotSimOutput(PI.data,i)
legend(observables(i))

end

finalValue=num2cell(exp(p_hat'));
[PI.par(1:end).finalValue]=finalValue{:,:};

%% Create variant object
MOC1_optimized = createVariant(PI,H,'MOC1_optimized');
variant = addvariant(model, 'MOC1_optimized');

addcontent(variant, MOC1_optimized.Content);
out.m1 = model;
sbiosaveproject '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/CIM_2.sbproj' 'model'