function plotFit(PI,varargin)
input= inputParser;
input.addParameter('sigma',repelem(1,1,size(PI.data(1).dataValue,2)))
input.addParameter('newFig', true)
input.parse(varargin{:})
input = input.Results;

% Obtain % of errors within error bands
% errors_10 = sum(arrayfun(@(x) sum(sum(abs(log10(x.dataValue./x.y_hat))<log10(1/0.9))),PI.data, 'UniformOutput', true))/PI.n_data*100;
% errors_20 = sum(arrayfun(@(x) sum(sum(abs(log10(x.dataValue./x.y_hat))<log10(1/0.8))),PI.data, 'UniformOutput', true))/PI.n_data*100;
% errors_50 = sum(arrayfun(@(x) sum(sum(abs(log10(x.dataValue./x.y_hat))<log10(1/0.5))),PI.data, 'UniformOutput', true))/PI.n_data*100;

errors_10 = sum(arrayfun(@(x) sum(sum(abs((x.dataValue-x.y_hat)./x.y_hat)<0.1)),PI.data, 'UniformOutput', true))/PI.n_data*100;
errors_20 = sum(arrayfun(@(x) sum(sum(abs((x.dataValue-x.y_hat)./x.y_hat)<0.2)),PI.data, 'UniformOutput', true))/PI.n_data*100;
errors_50 = sum(arrayfun(@(x) sum(sum(abs((x.dataValue-x.y_hat)./x.y_hat)<0.5)),PI.data, 'UniformOutput', true))/PI.n_data*100;

%
% Plot data points to fits
if input.newFig
figure
else
end
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
colors = linspecer(3);
% Add error reference lines
h = [];
h(1) = plot(rx,rx*1.11,'-', 'Color', colors(3,:),'LineWidth', 2);
plot(rx,rx*0.9, '-', 'Color', colors(3,:),'LineWidth', 2)

h(2) = plot(rx,rx*1.25,'--', 'Color', colors(1,:),'LineWidth', 2);
plot(rx,rx*.8,'--', 'Color', colors(1,:),'LineWidth', 2)

h(3) = plot(rx,rx*2,'-', 'Color', colors(2,:),'LineWidth', 2);
plot(rx,rx*0.5, '-', 'Color', colors(2,:),'LineWidth', 2)

set(gca,'XScale', 'log', 'YScale', 'log')
xlabel('Data')
ylabel('Fitted value')
title(strjoin({'Model fit to data (', PI.model, ')'},''))
ylim([xl xu])
xlim([xl xu])

errors = [{strjoin({'\beta = 10%: ', num2str(errors_10), '%'},'')} ...
    {strjoin({'\beta = 20%: ',num2str(errors_20), '%'},'')}...
    {strjoin({'\beta = 50%: ',num2str(errors_50), '%'},'')}];
lgnd = legend(h,errors,'Location', 'best');
title(lgnd, '% of model fits within \beta % of data value')


% figure
% hold on
% h=arrayfun(@(x)histogram((log(x.dataValue)-log(x.y_hat))./repmat(input.sigma,size(x.dataValue,1),1),'Normalization', 'probability'),PI.data,'UniformOutput',false)
% xlabel('Errors')
% ylabel('Probability')
% title(strjoin({'Residual error distribution (Log-transformed and normalized)', PI.model},' '))
% 
% figure
% hold on
% h=arrayfun(@(x)histogram((x.dataValue)-(x.y_hat),'Normalization', 'probability'),PI.data,'UniformOutput',false)
% xlabel('Errors')
% ylabel('Probability')
% title(strjoin({'Residual error distribution (Absolute errors)', PI.model},' '))


return
