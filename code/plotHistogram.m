function plotHistogram(x, paramNames)
ncol = ceil(sqrt(size(x,2)));
nrow = ceil(size(x,2)/ncol);
norm_pdf = @(x, mu, sigma) exp(-(x-mu).^2./(2*sigma.^2))./(2*pi*sigma);
for i=1:size(x,2)
subplot(nrow, ncol, i)
hold on
h=histogram(x(:,i), 'Normalization', 'pdf');
y = norm_pdf(x(:,i), mean(x(:,i)), std(x(:,i)));
[xsorted, sort_indx] = sort(x(:,i));


p = plot(xsorted,y(sort_indx));
title(paramNames{i})
end

legend({'Posterior samples', 'Normal approximation'})
% Q-Q plots
figure('Position', [10 10 1000 800])
for i=1:size(x,2)
subplot(nrow, ncol, i)
qqplot((x(:,i)-mean(x(:,i)))./std(x(:,i)));
ylabel('Posterior quantiles')
xlabel('Z quantiles')
title(strjoin({'Q-Q plot (' paramNames{i} ')'}, ''))
end
% tightfig()
end


