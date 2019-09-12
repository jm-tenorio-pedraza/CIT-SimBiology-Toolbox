function pars = plotPSS(pcs,pc_indx,parameters,varargin)

p=inputParser;
p.addParameter('threshold', -1);
p.parse(varargin{:})
p=p.Results;

figure;
colors=linspecer(pc_indx);
[~,I] = sort(log10(abs(pcs(:,1:pc_indx))),'descend');
pcs_sorted = log10(abs(pcs(I(:,1),1:pc_indx)));
hbar = bar(pcs_sorted,'FaceColor','flat');
hbar(1).BaseValue = -7;
haxes = hbar(1).Parent;
haxes.XTick = 1:length(parameters);
haxes.XTickLabel = parameters(I);
haxes.XTickLabelRotation = 30;
for k = 1:pc_indx
    hbar(k).CData = colors(k,:);
end
ylabel('Relative PC weight log10(|W_{ij}|/max(W_{ij}))');
title('Parameter sensitivity spectrum (PSS)')
phi = parameters(I(:,1));
pars = [];

hold on
plot(1:length(parameters), repelem(p.threshold,1,length(parameters)), '.-k')
hold off
legend_i={};
for i=1:pc_indx
pars(i).p_hat = phi(hbar(i).YData>p.threshold);
legend_i(i)={strjoin({'pc' num2str(i)},'')};
end
legend_i(end+1) = {'Threshold'};
legend(legend_i,'Location','NorthEastOutside');

return