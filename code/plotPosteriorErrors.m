function plotPosteriorErrors(PI,colIndx,varargin)
if nargin<2
    error('GWMCMC:toofewinputs','AMCMC requires atleast 2 inputs.')
end
p=inputParser;
p.addParameter('simTime',PI.tspan);
p.addParameter('outputs','indiv');
p.addParameter('all', false);
p.addParameter('central', 'median')
p.addParameter('newFig', true)
p.addParameter('TimeUnit', 'days')
p.addParameter('indiv', true);
p.addParameter('interpreter', 'tex');
p.addParameter('color', 'group');
p.parse(varargin{:});
p=p.Results;

% Plot all simulations or just those with data available?
if p.all
    n_sim=size(PI.data,1);
    sim_indx = ones(n_sim,1);
    n_col=ceil(sqrt(n_sim));
    n_row=ceil(n_sim/n_col);
else
    sim_indx=(arrayfun(@(x) ~all(isnan(x.dataValue(:,colIndx))), PI.data));
    n_sim = sum(sim_indx);
    n_col = ceil(sqrt(n_sim));
    n_row = ceil(n_sim/n_col);
end


treatments=([PI.data(sim_indx).Group]);
if length(treatments)~=n_sim
    treatments = unique({PI.data(:).Group});
else
    treatments = unique(treatments);
end
if strcmp(p.color, 'dataset')
    treatments = {PI.data(:).Name};
    treatment_colors=linspecer(length(treatments));

elseif strcmp(p.color, 'cell')
    cellIndx = 1:length(PI.CellTypes);
    treatments = PI.H.CellIndx*cellIndx';
    treatment_colors= PI.H.CellIndx*linspecer(length(cellIndx));
else
    treatment_colors=linspecer(length(treatments));

end

simIndx = find(sim_indx);

if p.newFig
    figure('Position', [10 10 900 1000])
end
% Identify which output to plot
output_i = strjoin({PI.observablesFields{colIndx} '_predError'},'');

try
    minX = min(arrayfun(@(x) min(min(x.(output_i))), PI.outputPred, 'UniformOutput', true));
    maxX = max(arrayfun(@(x) max(max(x.(output_i))), PI.outputPred, 'UniformOutput', true));
    ylimit = max(abs([minX maxX]));
catch
end
% For each valid simulation do:
for j=1:length(simIndx)
    if strcmp(p.outputs,'indiv')                                            % Extract prediction intervals only if individual outputs are desired
        subplot(n_row,n_col,j)
    else     
    end
        col_i=treatment_colors(simIndx(j),:);

    boxplot(PI.outputPred(simIndx(j)).(output_i),'Colors', col_i);
   hold on
     %% Axes specs
    ylabel(PI.variableUnits(colIndx),'interpreter', 'none')
    xlabel(strjoin({'Time [' p.TimeUnit ']'},''))
    if strcmp(p.outputs,'indiv')
        title(PI.data(simIndx(j)).Name,'interpreter','none')
    end
    ax = gca;                                                                   % Add title and legends

    ax.XTickLabel = num2str(PI.data(simIndx(j)).dataTime);
%     ylim([-ylimit, ylimit]);

end


legends = [PI.data(sim_indx).Name];
if length(legends)~= length(PI.data)
    legends = {PI.data(sim_indx).Name};
end
if strcmp(p.outputs,'indiv')
    %legend(ax.Children(end:-4:4),legends,'Interpreter', 'none','location','best')
else
    title(PI.observablesPlot(colIndx),'interpreter', p.interpreter)
    legend(legends,'Interpreter', 'none', 'location','best')
end

return