function PI = assignPrior(PI)
prior = repelem({'nan'}, length(PI.par));

prior(PI.H.PopulationParams) = {'U'};
prior(PI.H.SigmaParams) = {'IG'};
try
    prior([PI.H.IndividualParams.Index]) = {'N_z'};
catch
end
try
    prior([PI.H.CellParams.Index]) = {'N_z'};
catch
end

try
    prior([PI.H.RespParams.Index]) = {'N_z'};
catch
end

[PI.par(1:end).prior] = prior{:,:};
return
