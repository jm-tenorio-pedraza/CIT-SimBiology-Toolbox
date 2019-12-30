function plotFit(PI)

figure
hold on
f=arrayfun(@(x)plot(x.dataValue, x.y_hat,'*'), PI.data,'UniformOutput', false);
plot(1e-3:.1:10,1e-3:.1:10,'k')
set(gca,'XScale', 'log', 'YScale', 'log')
xlabel('Data')
ylabel('Fitted value')
title(strjoin({'Model fit to data', PI.model},' '))

figure
hold on
h=arrayfun(@(x)histogram(x.dataValue-x.y_hat,'Normalization', 'probability'),PI.data,'UniformOutput',false)
xlabel('Errors')
ylabel('Probability')
title(strjoin({'Residual error distribution', PI.model},' '))

return
