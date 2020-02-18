function logPrior = evaluateIndivPrior(p,PI)

psi = num2cell(repelem(exp(p([PI.H.CellParams.OmegaIndex])),1, length(PI.H.CellParams(1).Index)));
omega = num2cell(repelem(exp(p([PI.H.IndividualParams.OmegaIndex])), 1, length(PI.H.IndividualParams(1).Index)));
sigma = [psi omega];
p = num2cell(p);

[PI.par(1:end).p] = p{:,:};
[PI.par([PI.H.CellParams.Index PI.H.IndividualParams.Index]).sigma] = sigma{:,:};
logPrior = arrayfun(@(x) x.priorHandle(x.p,x.sigma),...
    PI.par([PI.H.CellParams.Index PI.H.IndividualParams.Index]));
logPrior = sum(log(logPrior));
return
