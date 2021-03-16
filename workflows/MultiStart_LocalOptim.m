%% Multi-start local optimization
% Set seed
rng('default')
s = rng;

%% Set residuals function
popParamsIndx = [PI.H.PopulationParams];
cellParamsIndx = [PI.H.CellParams.Index];
indivParamsIndx = [PI.H.IndividualParams.Index];
respParamsIndx = [PI.H.RespParams.Index];

options_lsqnonlin=optimoptions('lsqnonlin','Display','iter','MaxFunEvals',...
    5e4, 'MaxIter',5e4, 'TolFun', 1e-4, 'StepTolerance', 1e-6);
ub = log([PI.par([popParamsIndx cellParamsIndx indivParamsIndx respParamsIndx]).maxValue]);
lb = log([PI.par([popParamsIndx cellParamsIndx indivParamsIndx respParamsIndx]).minValue]);


phat_0 = finalValues;
% Local Optimization

% Residuals 
residuals_fn = @(x) getResiduals(exp(x),@(x)sim(x,PI.tspan(end),PI.u,PI.tspan),PI,...
    @(x)getPhi2(x,PI.H,length(PI.u),'initialValue',PI.x_0),...
    exp(phat_0(end-length(observables)+1:end)),...
    exp(phat_0([PI.H.CellParams.OmegaIndex])),...
    exp(phat_0([PI.H.IndividualParams.OmegaIndex])),...
    exp(phat_0([PI.H.RespParams.OmegaIndex])),PI.normIndx);
%% 
n_samples = 200;
lhs = lhsdesign(n_samples,length(ub));
p0 = unifinv(lhs,repelem(exp(lb),n_samples,1),repelem(exp(ub),n_samples,1));
p_hat = nan(size(p0)); fval_lsqnonlin = nan(size(p0,1),1);
for i = 1:size(p0,1)
    try
    [p_hat(i,:), fval_lsqnonlin(i)] = lsqnonlin(residuals_fn,log(p0(i,:)), lb,ub, options_lsqnonlin);
    catch
        p_hat(i,:)= p0(i,:);
        fval_lsqnonlin(i) = 1e9;
    end
end
[fval_sorted,I]=sort(fval_lsqnonlin);

% Sort parameter samples wrt rss
p_hat_sorted = p_hat(I,:);
pstart_sorted = p0(I,:);
%% Add results to PI
p0_cell = mat2cell(pstart_sorted, ones(size(pstart_sorted,1),1));
p_hat_cell =  mat2cell(p_hat_sorted, ones(size(p_hat_sorted,1),1));
fval_cell = num2cell(fval_sorted);
PI.MultiStartOptim = struct('p_start', p0_cell, 'p_end', p_hat_cell, 'fval_lsqnonlin', fval_cell);
%% Plot results
figure
hold on
ax =gca;
cutoffValue = find([PI.MultiStartOptim(:).fval_lsqnonlin]>[PI.MultiStartOptim(1).fval_lsqnonlin]*2,1);
set(ax, 'YScale', 'log')
plot([PI.MultiStartOptim(:).fval_lsqnonlin])
plot(repelem(PI.MultiStartOptim(cutoffValue).fval_lsqnonlin, length(PI.MultiStartOptim),1));
%% Select acceptable samples
X0 = cell2mat({PI.MultiStartOptim(1:cutoffValue).p_end}');
cv = std(exp(X0))./abs(mean(exp(X0)));
plotHistogram(X0, paramNames([PI.H.PopulationParams PI.H.CellParams.Index PI.H.IndividualParams.Index]))

%% Add best results to PI.par
x0_cell = num2cell(exp(X0(1,:)));
mean_cell = num2cell(mean(exp((X0))));
cv_cell = num2cell(cv);
[PI.par([PI.H.PopulationParams PI.H.CellParams.Index PI.H.IndividualParams.Index]).phat_MultiStartOptim] = (x0_cell{:,:});
[PI.par([PI.H.PopulationParams PI.H.CellParams.Index PI.H.IndividualParams.Index]).mean_MultiStartOptim] = (mean_cell{:,:});

[PI.par([PI.H.PopulationParams]).cv_MultiStartOptim] = (cv_cell{:,:});

%% Use selected samples as initial walkers for MCMC
sigma_0 = repmat(finalValues([H.SigmaParams]),sum(cutoffIndx),1);
X0=[log(finalValues);X0 log(sigma_0)];

