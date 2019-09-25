%% Multi-start local optimization

ub = log([PI.par([H.PopulationParams H.IndividualParams.Index]).maxValue]);
lb = log([PI.par([H.PopulationParams H.IndividualParams.Index]).minValue]);

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
