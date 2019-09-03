function [logpdf, gradlogpdf]= logPosterior(param,loglikelihood, logprior,PI, H)
%% Log-posterior
logpdf=loglikelihood(param)+logprior(param);



%% Gradient

% Extract indexes from normally-distributed parameters
indiv_indx=arrayfun(@(x)(x.Index),H.IndividualParams,'UniformOutput',false);
indiv_indx=reshape(cell2mat(indiv_indx),1,[]);
norm_indx=[H.PopulationParams indiv_indx];
norm_gradlogpdf=-(param(norm_indx)-log([PI.par(norm_indx).mu_prior]))./[PI.par(norm_indx).sigma_prior].^2;

% Lognormal prior for the standar deviations
sigma_gradlogpdf = -(param(H.SigmaParams)-log([PI.par(H.SigmaParams).mu_prior]))./[PI.par(H.SigmaParams).sigma_prior].^2;

gradlogpdf=NaN(size(param));
gradlogpdf(norm_indx)=norm_gradlogpdf;
gradlogpdf(H.SigmaParams)=sigma_gradlogpdf;
return