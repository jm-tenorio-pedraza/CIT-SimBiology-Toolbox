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
indx=1:n_par;

figure
for i=1:n_par
    if strcmp(p.algorithm,'MCMC')
        central=PI.par(i).posterior_mean;
    else
        central=PI.par(i).mean_boot;
    end
    h=plot(([PI.par(i).LB central PI.par(i).UB ]), repelem(indx(i),1,3), 'd');
    h.LineStyle='-';
    h.Color =colors(i,:);
    h.MarkerFaceColor=colors(i,:);
    h.MarkerEdgeColor=colors(i,:);
    legend_i=strjoin({num2str(central,3), ' ', '[',num2str((PI.par(i).LB),3), ',',num2str((PI.par(i).UB),3),']'},'');
    text(PI.par(i).LB, i+0.5, legend_i)
    hold on
    
end
set(gca,'XScale', 'log','YAxisLocation', 'right','YTick', indx,'YTickLabel', names, 'TickLabelInterpreter', p.interpreter)
% xlim([1e-16 1e2])
title(strjoin({p.confidence 'Credible intervals for', model, '(' p.method, ')'},' '))

