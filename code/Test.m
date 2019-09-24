%% Code for DREAM MCMC testing using multivariat normal with random effects

% Model = Y = XB + ZN + E
% Y is an mx7 matrix for 7 different variables
% X is an mx20 matrix
% B is a 20X7 matrix of coefficients
% Z is an mX1 vector of random effects where we assume the same random
% effect affects all variables similarly
% E is an mX7 matrix of measurement errors 
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
rng(1)
%%
H.PopulationParams = 1:20;
H.IndividualParams.Index = 21:30;
H.IndividualParams.EtaIndex = 20;
H.IndividualParams.OmegaIndex = 31;
H.SigmaParams = 31:39;

m = 200;
k = 8; % number of variables measured
b = length(H.PopulationParams); % Number of fixed effects
z = 20; % Number of observations per individual
n = m/z; % Number of distinct observations (random effects) number of individuals
B_sigma = exp(-rand(b,1));
B = repmat(exp(rand(b,1).*B_sigma),1,k); % weight of each dimension i on each variable j
O = exp(-rand(1))*2;
N = repmat(exp(log(B(end,1))+randn(n,1).*O),1,k); % random effect of each individual
Z = kron(eye(n),ones(z,1));
S = exp(-rand(1,k));
E = randn(m,k).*S;  % Error in measurements one sigma for each k variable
X = randn(m,b); % Random uniform weights of each variable on each observation
Y = exp(X*log(B) + Z*log(N) + E);
p = [B(:,1); N(:,1); O; S'];
XB=X*log(B);
ZN = Z*log(N);
names = [repelem({'b'}, b, 1); repelem({'eta'}, n,1); repelem({'sigma'}, k+1,1)];

%%
ncol = ceil(sqrt(k));
nrow = ceil(k/ncol);
figure;
for i=1:k
    for j=1:n
        subplot(nrow,ncol, i)
        histogram(log(Y((j-1)*n+1:(j-1)*n+z,i)),'Normalization', 'probability')
        hold on
    end
end

%%
likelihood = @(p) (mvnloglikelihood((p),H,log(Y),X,Z)*(-1));
mu = p([H.PopulationParams H.SigmaParams])+randn(length([H.PopulationParams H.SigmaParams]),1);
prior = @(x) (sum(((x([H.PopulationParams H.SigmaParams])-mu').^2/10)...
    +log(sqrt(pi*10)))*(-1)); % sigma for all parameters' priors is sqrt(10/2)
%likelihood= @(p) (prior(p) + likelihood(p));

p0 = (randn(b+n+k+1,b+n+k+1));

[x, p_x,accept] = dreamH(p0,likelihood,prior,size(p0,1),4e3,size(p0,2),'BurnIn',...
    1e4,'H', H);
[x, p_x,accept] = dreamHParallel(p0,likelihood,prior,size(p0,1),4e3,size(p0,2),'BurnIn',...
    1e5,'H', H);
[x2, p_x2,accept2] = dream(p0,likelihood,prior,size(p0,1),4e3,size(p0,2),'BurnIn',...
    1e5);

plot(-p_x2)
set(gca,'YScale', 'log')

%%
plotMCMCDiagnostics(x,p_x,'name', names)
x_mat = x(:,:)';
plotBivariateMarginals_2(x_mat(1e5:1e3:end,H.PopulationParams(1:10)),'names', names(H.PopulationParams(1:10)))
plotBivariateMarginals_2(x_mat(1e5:1e3:end,H.PopulationParams(11:20)),'names', names(H.PopulationParams(11:20)))
plotBivariateMarginals_2(x_mat(1e5:1e3:end,H.IndividualParams.Index),'names', names(H.IndividualParams.Index))
plotBivariateMarginals_2(x_mat(1e5:1e3:end,H.SigmaParams),'names', names(H.SigmaParams))

plotCorrMat(x_mat(1e5:1e3:end,:),names)
%% Comparing estimates with true values

phat = mean(x_mat(1e5:1e3:end,:))';
p_hat=table(p,exp(phat));
loglog(p,exp(phat),'+')
text(p,exp(phat), names)
hold on
loglog(min(p):max(p), min(p):max(p),'-k')
norm(p-exp(phat))