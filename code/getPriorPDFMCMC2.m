function priorPDF = getPriorPDFMCMC2(p,PI)

prior = {PI.par(:).prior};
priorPDF = nan(length(prior),1);
U_indx = ismember(prior, 'U');
N_indx = ismember(prior, 'N');

% Prior functions
N = @(x,mu,sigma) exp(-(x-mu).^2./(2*sigma.^2))./(2*pi*sigma);
N_z = @(x,sigma) exp(-(x.^2)./(2*sigma.^2))./(2*pi*sigma);
U = @(x,a,b) (x>=a && x<=b);
IG = @(x,a,b) x.^(-a-1).*exp(-b./x);

% Expand sigma parameters of IEV and IIV parameters
psi = reshape((repelem(p([PI.H.CellParams.OmegaIndex])', length(PI.H.CellParams(1).Index),1)),[],1);
omega = (repelem(p([PI.H.IndividualParams.OmegaIndex])', length(PI.H.IndividualParams(1).Index),1));
nu =  reshape((repelem(p([PI.H.RespParams.OmegaIndex])', length(PI.H.RespParams(1).Index),1)),[],1);
% Transform p vector into a cell array and assing to PI.par structure
curr_p = num2cell(p);
[PI.par(1:end).curr_p] = curr_p{:,:};

% Assing the current value of the sigma for IEV and IIV parameters as the
% sigma_prior
psi = num2cell(psi);
omega = num2cell(omega);
nu = num2cell(nu);
[PI.par([PI.H.CellParams.Index]).sigma_prior] = psi{:,:};
[PI.par([PI.H.IndividualParams.Index]).sigma_prior] = omega{:,:};
[PI.par([PI.H.RespParams.Index]).sigma_prior] = nu{:,:};

% Evaluate function
priorPDF(U_indx,1) = arrayfun(@(y) U(y.curr_p, y.minValue, y.maxValue), PI.par(U_indx));
priorPDF(N_indx,1) = arrayfun(@(y) N(log(y.curr_p), log(y.mu_prior), y.sigma_prior), PI.par(N_indx));
priorPDF([PI.H.CellParams.Index]) = arrayfun(@(y) N_z(log(y.curr_p), y.sigma_prior^2), PI.par([PI.H.CellParams.Index]));
priorPDF([PI.H.IndividualParams.Index]) = arrayfun(@(y) N_z(log(y.curr_p), y.sigma_prior^2), PI.par([PI.H.IndividualParams.Index]));
priorPDF([PI.H.RespParams.Index]) = arrayfun(@(y) N_z(log(y.curr_p), y.sigma_prior^2), PI.par([PI.H.RespParams.Index]));

priorPDF(PI.H.SigmaParams,1) = arrayfun(@(y) IG(y.curr_p.^2, y.mu_prior, y.sigma_prior), PI.par(PI.H.SigmaParams));

priorPDF = sum(log(priorPDF));
return



