function R_hat = getGelmanRubinStatistic(x,parameters,varargin)
T = size(x,3);
N = size(x,2);
d = size(x,1);
p = inputParser;
p.addParameter('threshold', 1.2);
p.parse(varargin{:});

p=p.Results;

x_bar = 2*sum(x(:,:,ceil(T/2):end), 3)/(T-2);
x_doublebar = mean(x_bar,2);
%W = sum(sum(x(:,:,ceil(T/2):end).^2,3) - sum(x(:,:,ceil(T/2):end).*x_bar,3),2) + ceil(T/2)*sum(x_bar.^2,2);
W = 2/(N*(T-2))*sum(sum((x(:,:,ceil(T/2):end)-x_bar).^2,3),2);
B = 1/(2*(N-1))*sum((x_bar-x_doublebar).^2,2)/T;

sigma_plus = (T-2)/T*W+2*B;
R_hat = sqrt((N+1)/N*sigma_plus./W - (T-2)/(N*T));
colors = linspecer(length(R_hat));
figure;
hbar = bar(R_hat','FaceColor','flat');
hbar(1).BaseValue = min(min(R_hat));
haxes = hbar(1).Parent;
haxes.XTick = 1:length(R_hat);
haxes.XTickLabel = parameters;
haxes.XTickLabelRotation = 30;
% for k = 1:length(R_hat)
%     hbar(k).CData = colors(k,:);
% end
hold on
plot(1:length(parameters), repelem(p.threshold,1,length(parameters)), '.-k')
set(gca,'YScale','log')
return


