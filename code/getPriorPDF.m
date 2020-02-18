function priorPDF = getPriorPDF(p,PI, prior)

priorSigma = repelem({'IG'}, 1,length(PI.H.SigmaParams));
priorZ = repelem({'N_z'}, 1, length([PI.H.IndividualParams.Index]));
priorW = repelem({'N_z'}, 1, length([PI.H.CellParams.Index]));

prior = [prior priorW priorZ priorSigma];

N = @(x,mu,sigma) exp(-(x-mu).^2./(2*sigma.^2))./(2*pi*sigma);
N_z = @(x,sigma) exp(-(x.^2)./(2*sigma.^2))./(2*pi*sigma);

U = @(x,a,b) prod(and(x>=a, x<=b));
IG = @(x,a,b) b.^a./gamma(a).*x.^(-a-1).*exp(-b./x);

mu = [PI.par(:).mu_prior];
sigma = [PI.par(:).sigma_prior];
lb = [PI.par(:).minValue];
ub = [PI.par(:).maxValue];
priorPDF = nan(length(prior),1);

U_indx = ismember(prior, 'U');
N_indx = ismember(prior, 'N');

psi = reshape(exp(repelem(p([PI.H.CellParams.OmegaIndex]), length(PI.H.CellParams(1).Index),1)),[],1);
omega = exp(repelem(p([PI.H.IndividualParams.OmegaIndex]), length(PI.H.IndividualParams(1).Index),1));

priorPDF(U_indx,1) = U(p(U_indx), log(lb(U_indx)),log(ub(U_indx)));
priorPDF(N_indx,1) = N(p(N_indx), log(mu(N_indx)), sigma(N_indx));
priorPDF(PI.H.SigmaParams,1) = IG(exp(p(PI.H.SigmaParams)).^2, mu(PI.H.SigmaParams), sigma(PI.H.SigmaParams));
priorPDF([PI.H.CellParams.Index]) = N_z(p([PI.H.CellParams.Index])', psi);
priorPDF([PI.H.IndividualParams.Index]) = N_z(p([PI.H.IndividualParams.Index])', omega);

priorPDF = sum(log(priorPDF));
return



