function plotMCMCDiagnostics(params, logP, varargin)
if nargin<3
    error('AMCMCH:toofewinputs','AMCMC requires atleast 3 inputs.')
end
p=inputParser;
p.addParameter('model','',@ischar);
p.addParameter('name',{});
p.addParameter('plots',{'trace' 'autocorr' 'corr'});
p.addParameter('interpreter','none');

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
[C,lags,ESS]=eacorr(params);
% LogL traceplot
if ismember('trace', p.plots)
figure
h=plot(-logP);
xlabel('MCMC step')
ylabel('Negative Log-Likelihood')
title(strjoin({'Log-posterior density trace plot of',p.model},' '), 'interpreter', p.interpreter)
set(gca, 'YScale', 'log')
colors = linspecer(size(logP,2));
for i=1:size(logP,2)
    h(i).Color = colors(i,:);
end
% Param traceplot
figure('Renderer', 'painters', 'Position', [10 10 1500 600])
plotTrace(params,'names', p.name,'ESS',ESS,'interpreter', p.interpreter)
end
% Autocorrelation
if ismember('autocorr', p.plots)
figure('Renderer', 'painters', 'Position', [10 10 1500 600])


ncol = ceil(sqrt(length(ESS)));
nrow = ceil(length(ESS)/ncol);
colors = linspecer(2);
for i=1:length(ESS)
    subplot(nrow,ncol,i)
    plot(lags,C(:,i))
    hold on
    plot([lags(1) lags(end)],zeros(1,2),'-k', 'Color', colors(1,:))
    ylim([-1 1])
    tau = cumsum(C(:,i));
    indx = find(tau*5<(1:length(tau))');
    title(p.name{i},'Interpreter', p.interpreter)
    try
    legend (strjoin({'\tau =' num2str(1+2*sum(C(1:indx(1),i)))},' '))
    catch
       legend (strjoin({'\tau =' num2str(1+2*sum(C(1:end,i)))},' '))

    end
end
end

% Correlation matrix
if ismember('corr', p.plots)
figure('Renderer', 'painters', 'Position', [10 10 1500 600])
plotCorrMat(phat, p.name, 'model', p.model,'interpreter', p.interpreter)
end
end