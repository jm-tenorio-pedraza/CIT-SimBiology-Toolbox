function plotCI(PI, model, varargin)
if nargin<1
    error('plotCI:toofewinputs','plotCI requires atleast 2 inputs.')
end
p=inputParser;
p.addParameter('method','HPD');
p.addParameter('confidence','95%');
p.addParameter('algorithm','MCMC');
p.addParameter('name', {PI.par(:).name});
p.addParameter('interpreter', 'none')
p.parse(varargin{:});
p=p.Results;

names=p.name;
n_par=length(names);
colors=linspecer(n_par);
indx=n_par:-1:1;

figure
for i=1:n_par
    if strcmp(p.algorithm,'MCMC')
        central=PI.par(indx(i)).posterior_mean;
    else
        central=PI.par(indx(i)).mean_boot;
    end
    h=plot(([PI.par(indx(i)).LB central PI.par(indx(i)).UB ]), repelem(indx(i),1,3), 'd');
    h.LineStyle='-';
    h.Color =colors(indx(i),:);
    h.MarkerFaceColor=colors(indx(i),:);
    h.MarkerEdgeColor=colors(indx(i),:);
    legend_i=strjoin({num2str(central,3), ' ', '[',num2str((PI.par(indx(i)).LB),3), ',',num2str((PI.par(indx(i)).UB),3),']'},'');
    text(PI.par(indx(i)).LB, i+0.5, legend_i)
    hold on
    
end
set(gca,'XScale', 'log','YAxisLocation', 'right','YTick', 1:n_par,'YTickLabel', names(indx), 'TickLabelInterpreter', p.interpreter)
% xlim([1e-16 1e2])
title(strjoin({p.confidence 'Credible intervals for', model, '(' p.method, ')'},' '))

