function logPrior = evaluatePrior(p,PI)

p = num2cell(p);
[PI.par(1:end).p] = p{:,:};
logPrior = arrayfun(@(x) x.priorHandle(x.p),PI.par([PI.H.PopulationParams PI.H.SigmaParams]));
logPrior = sum(log(logPrior));
return
