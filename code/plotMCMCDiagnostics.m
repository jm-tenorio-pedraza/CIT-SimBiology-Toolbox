function plotMCMCDiagnostics(params, logP, varargin)
if nargin<3
    error('AMCMCH:toofewinputs','AMCMC requires atleast 3 inputs.')
end
p=inputParser;
p.addParameter('model','',@ischar);
p.addParameter('name',{});

if sum(imag(logP)>0)>0
    imag_indx=sum(imag(logP)>0);
    fprintf('Warning: There were %.f imaginary numbers in the log-likelihood',imag_indx)
    params=params(~imag(logP)>0,:);
    logP=logP(~imag(logP)>0);
else
   

end
p.addParameter('steps',1:1:size(params,1));
p.addParameter('BurnIn',1e5);
p.addParameter('Thinning',40);
p.addParameter('AdaptSteps',100);
p.addParameter('initialSigma', repelem(0.2/size(params,2),1, size(params,2)));

p.parse(varargin{:});
p=p.Results;

if size(params,3)>1
    phat=params(:,:)';
else
    phat=params;
end
if size(logP,3)>1
    logP=sum(logP(:,:))';
end
% LogL traceplot
figure
plot(-logP)
xlabel('MCMC step')
ylabel('Negative Log-Likelihood')
title(strjoin({'Log-posterior trace plot of',p.model},''), 'interpreter', 'none')
set(gca, 'YScale', 'log')

% Autocorrelation
figure
[C,lags,ESS]=eacorr(params);
plot(lags,C)
ylim([-1 1])
title('Autocorrelation of posterior samples of')

% Param traceplot
figure('Renderer', 'painters', 'Position', [10 10 1500 600])
plotTrace(phat,'names', p.name,'ESS',ESS)

% Correlation matrix
figure
plotCorrMat(phat, p.name, 'model', p.model)
end