function plotFit(PI,varargin)
input= inputParser;
input.addParameter('sigma',repelem(1,1,size(PI.data(1).dataValue,2)))
input.parse(varargin{:})
input = input.Results;

figure
hold on
try
    f=arrayfun(@(x)errorbar(x.dataValue, x.y_hat,x.SD, 'horizontal','s'), PI.data,'UniformOutput', false);

catch
    f=arrayfun(@(x)plot(x.dataValue, x.y_hat,'*'), PI.data,'UniformOutput', false);
end
ax = gca;
x_lim = ax.XLim;
y_lim = ax.YLim;
xl = min([x_lim(1) y_lim(1)]);
xu = max([x_lim(2) y_lim(2)]);
if xl==0
    xl = 1e-3;
end
rx = ((xl):((xu-xl)/1e3):(xu));
plot(rx,rx,'k')
plot(rx,rx*10,'--k')
plot(rx,rx*.1,'--k')
% plot(rx,rx*1.2,'--k')
% plot(rx,rx*0.8, '--k')
% plot(rx,rx*1.5,'--g')
% plot(rx,rx*0.5, '--g')
set(gca,'XScale', 'log', 'YScale', 'log')
xlabel('Data')
ylabel('Fitted value')
title(strjoin({'Model fit to data (', PI.model, ')'},''))
ylim([xl xu])
xlim([xl xu])
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
