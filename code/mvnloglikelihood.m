function loglikelihood = mvnloglikelihood(p, H,Y,X,Z)
fixed_index = H.PopulationParams;
sigma_index = setdiff(H.SigmaParams, [H.IndividualParams.OmegaIndex]);
omega_index = H.IndividualParams.OmegaIndex;
eta_index = [H.IndividualParams(:).Index];
if size(p,1)>size(p,2)
else
    p=p';
end
mu = repmat(p(fixed_index),1,size(Y,2));
eta = repmat(p(eta_index),1,size(Y,2));
sigma = exp(p(sigma_index))';
omega = exp(p(omega_index));
loglikelihood =  size(X,1)*sum(log(sigma)+log(sqrt(2*pi)))+sum(sum((X*mu+Z*eta-Y).^2./(2*sigma.^2)))+...
    sum(log(omega)+log(sqrt(pi*2)) + ((p(eta_index) - p(H.PopulationParams(end))).^2./(2*omega^2)));
return


