function plotFit(PI,varargin)
input= inputParser;
input.addParameter('sigma',repelem(1,1,size(PI.data(1).dataValue,2)))
input.parse(varargin{:})
input = input.Results;

figure
hold on
f=arrayfun(@(x)plot(x.dataValue, x.y_hat,'*'), PI.data,'UniformOutput', false);
plot(1e-3:.1:10,1e-3:.1:10,'k')
plot(1e-3:.1:10,(1e-3:.1:10)*10,'--k')
plot(1e-3:.1:10,(1e-3:.1:10)*.1,'--k')

set(gca,'XScale', 'log', 'YScale', 'log')
xlabel('Data')
ylabel('Fitted value')
title(strjoin({'Model fit to data', PI.model},' '))

figure
hold on
h=arrayfun(@(x)histogram((log(x.dataValue)-log(x.y_hat))./repmat(input.sigma,size(x.dataValue,1),1),'Normalization', 'probability'),PI.data,'UniformOutput',false)
xlabel('Errors')
ylabel('Probability')
title(strjoin({'Residual error distribution (Log-transformed and normalized)', PI.model},' '))

figure
hold on
h=arrayfun(@(x)histogram((x.dataValue)-(x.y_hat),'Normalization', 'probability'),PI.data,'UniformOutput',false)
xlabel('Errors')
ylabel('Probability')
title(strjoin({'Residual error distribution (Absolute errors)', PI.model},' '))


return
