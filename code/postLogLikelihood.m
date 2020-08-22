function postLogLikelihood = postLogLikelihood(p, prior, likelihood)
prior_pdf = prior(p);
if isinf(prior_pdf) || isnan(prior_pdf) || ~isreal(prior_pdf)
    postLogLikelihood = -1e9;
    return
else
    postLogLikelihood = prior_pdf + likelihood(p);
end
return