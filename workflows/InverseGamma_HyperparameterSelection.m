%% Inverse gamma prior
invgamma_prior= @(x,a,b)((x.^(-a-1).*exp(-b./x)));
x = 0.01:0.01:5;
y = invgamma_prior(x.^2,.01,.01);
plot(x,-log(y))
hold on
log(y(1))
log(y(end))
% set(gca, 'YScale', 'log')
%%
a = 0.01:0.1:1;
b = 0.01:0.1:1;
[X, Y] = meshgrid(a,b);
colors = linspecer(size(Y,2));
ncol = ceil(sqrt(size(X,2)));
nrow = ceil(size(X,2)/ncol);
indx = 1;
for j = 1:size(X,2)
    for i = 1:size(X,1)
        y = invgamma_prior(x.^2, X(i,j), Y(i,j));
        subplot(nrow,ncol, indx)
        hold on
        plot(x,y/max(y), 'Color', colors(i,:));
    end
    child = get(gca, 'Children');
    legend(cellfun(@(x)strjoin({'b = ' x},''), cellstr(num2str(Y(:,j))), 'UniformOutput', false))
    hold off
    title(strjoin({'a = ', num2str(X(1,j))}, ''))
    indx = indx + 1;
end