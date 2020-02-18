function PI = getPriorDistribution(PI, prior)

priorSigma = repelem({'IG'}, 1,length(PI.H.SigmaParams));
priorZ = repelem({'N_z'}, 1, length([PI.H.IndividualParams.Index]));
priorW = repelem({'N_z'}, 1, length([PI.H.CellParams.Index]));

prior = [prior priorW priorZ priorSigma];

N = @(x,mu,sigma) exp(-(x-mu)^2/(2*sigma^2))/(2*pi*sigma);
N_z = @(x,sigma) exp(-(x)^2/(2*sigma^2))/(2*pi*sigma);

U = @(x,a,b) prod(x>a && x<=b);
IG = @(x,a,b) b^a/gamma(a)*x^(-a-1)*exp(-b/x);

mu = [PI.par(:).mu_prior];
sigma = [PI.par(:).sigma_prior];
lb = [PI.par(:).minValue];
ub = [PI.par(:).maxValue];
for i=1:length(prior)
    prior_i = prior{i};
    switch prior_i
        case 'N'
            PI.par(i).priorHandle = @(x)N(x,log(mu(i)), sigma(i));
        case 'N_z'
            PI.par(i).priorHandle = N_z;
        case 'U'
            PI.par(i).priorHandle = @(x)U(x,log(lb(i)), log(ub(i)));
        case 'IG'
            PI.par(i).priorHandle = @(x)IG(exp(x)^2, mu(i), sigma(i));
    end
end
return