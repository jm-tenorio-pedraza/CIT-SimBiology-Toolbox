function theta_hat = samplePrior(p,PI,varargin)
% Inputs:
%   p: is a vector of probabilities, e.g., the output from LHSdesign
%   PI: is the parameter identification structure with PI.par and H
%   structures
% Outputs:
%   theta_hat : vector of samples from prior distribution

par=inputParser;
par.addParameter('LB', [])
par.parse(varargin{:})
par = par.Results;

prior = {PI.par(:).prior};
theta_hat = nan(1,length(prior));
U_indx = ismember(prior, 'U');
N_indx = ismember(prior, {'N'});
IG_indx = ismember(prior, 'IG');

% Transform p vector into a cell array and assing to PI.par structure
curr_p = num2cell(p);
[PI.par(1:end).curr_p] = curr_p{:,:};


theta_hat(U_indx) = arrayfun(@(y) unifinv(y.curr_p, y.minValue, y.maxValue), PI.par(U_indx));
theta_hat(N_indx) = arrayfun(@(y) norminv((y.curr_p), log(y.mu_prior), y.sigma_prior), PI.par(N_indx));
theta_hat(IG_indx) = arrayfun(@(y) sqrt(1./gaminv(y.curr_p,1./y.mu_prior, y.sigma_prior)),PI.par(IG_indx));

% Expand sigma parameters of IEV and IIV parameters
psi = reshape((repelem(theta_hat([PI.H.CellParams.OmegaIndex])', length(PI.H.CellParams(1).Index),1)),[],1);
omega = (repelem(theta_hat([PI.H.IndividualParams.OmegaIndex])', length(PI.H.IndividualParams(1).Index),1));
nu =  reshape((repelem(theta_hat([PI.H.RespParams.OmegaIndex])', length(PI.H.RespParams(1).Index),1)),[],1);
% Assing the current value of the sigma for IEV and IIV parameters as the
% sigma_prior
psi = num2cell(psi);
omega = num2cell(omega);
nu = num2cell(nu);
[PI.par([PI.H.CellParams.Index]).sigma_prior] = psi{:,:};
[PI.par([PI.H.IndividualParams.Index]).sigma_prior] = omega{:,:};
[PI.par([PI.H.RespParams.Index]).sigma_prior] = nu{:,:};

theta_hat([PI.H.CellParams.Index]) = arrayfun(@(y) exp(norminv((y.curr_p), (y.sigma_prior))), PI.par([PI.H.CellParams.Index]));
theta_hat([PI.H.IndividualParams.Index]) = arrayfun(@(y) exp(norminv((y.curr_p), (y.sigma_prior))), PI.par([PI.H.IndividualParams.Index]));
theta_hat([PI.H.RespParams.Index]) = arrayfun(@(y) exp(norminv((y.curr_p), (y.sigma_prior))), PI.par([PI.H.RespParams.Index]));

end





