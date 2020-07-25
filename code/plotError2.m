function plotError2(PI,colIndx, varargin)
par=inputParser;
par.addParameter('all', true)
par.addParameter('indiv', true)
par.addParameter('addErrorVar', true)
par.addParameter('newFig', true)
par.addParameter('interpreter', 'tex')
par.addParameter('group', 'byDataset')
par.addParameter('TimeUnit', 'days')

par.parse(varargin{:})
par = par.Results;

% Plot all simulations or just those with data available?
if par.all
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
if strcmp(par.group, 'byDataset')
    treatments=([PI.data(sim_indx).Group]);
    treatment_colors=linspecer(length(treatments));
    if length(treatments)~=n_sim
        treatments = unique({PI.data(:).Group});
    else
        treatments = unique(treatments);
    end
elseif strcmp(par.group, 'Cell')
    treatments = PI.H.CellIndx(:,:);
    treatment_colors = PI.H.CellIndx(:,:)*linspecer(size(PI.H.CellIndx,2));
  
end

if par.indiv
else
    treatments = {PI.data(:).Name};
end

    simIndx = find(sim_indx);

try
    minX = min(arrayfun(@(x) min(x.dataValue(:,colIndx)-x.y_hat(:,colIndx)), PI.data, 'UniformOutput', true));
    maxX = max(arrayfun(@(x) max(x.dataValue(:,colIndx)-x.y_hat(:,colIndx)), PI.data, 'UniformOutput', true));
catch
    
end
if par.newFig
figure('Position', [10 10 900 1000])
end
for i=1:n_sim
    if par.indiv
        subplot(n_row,n_col,i)
    else
    end
    hold on

        dat=plot(PI.data(simIndx(i)).dataTime, PI.data(simIndx(i)).dataValue(:,colIndx) - PI.data(simIndx(i)).y_hat(:,colIndx));
    if par.indiv
            title(PI.data(simIndx(i)).Name,'interpreter', par.interpreter)

        col_i=treatment_colors(ismember(treatments,PI.data(simIndx(i)).Group),:);
    else
        title(PI.observablesPlot(colIndx), 'interpreter', par.interpreter)
        col_i = treatment_colors(simIndx(i),:);
    end
     dat.LineStyle='none';
    dat.Color = col_i;
    dat.MarkerFaceColor=col_i;
    dat.MarkerEdgeColor=col_i;
    dat.Marker='d';
   
    ylabel(PI.variableUnits(colIndx))
    if strcmp(par.TimeUnit, 'days')
    xlabel('Time [days]')
    elseif strcmp(par.TimeUnit, 'hours')
        xlabel('Time [hours]')
    end
    ax = gca;
    if par.indiv
    try
        legend(ax.Children(3:1),{PI.observablesPlot{colIndx} 'Data' 'SD'}, 'location', 'best', 'interpreter','none')
    catch
        legend(ax.Children(2:1),{PI.observablesPlot{colIndx} 'SD'}, 'location', 'best', 'interpreter', 'none')
    end
    end
    %
    % Change axis if there is large variability in output
    if std(log10(PI.data(i).simOutput(:,colIndx)))>2
        %         set(gca,'YScale','log')
    end
    % set(gca,'YScale','log')
    %        ylim([0.1, 100])
    try
        %ylim(10.^([floor(log10(min(PI.data(simIndx(i)).dataValue(:,colIndx)))) ceil(log2(max(PI.data(simIndx(i)).dataValue(:,colIndx))))]))
    catch
    end
    ylimit = max(abs([minX maxX]));
    
    ylim([floor(-ylimit), ceil(ylimit)])
    
  
end
        abline = plot(PI.data(simIndx(i)).dataTime, zeros(1, length(PI.data(simIndx(i)).dataTime)));

    abline.Color = 'k';

ax = gca;
if ~par.indiv
    ncells = 1:size(PI.H.CellIndx,2);
    dataLegend = PI.H.CellTypes(PI.H.CellIndx(simIndx,:)*ncells');
    [dataLegend_Unique, uniqueIndx] = unique(dataLegend,'stable');
    legendIndx = uniqueIndx+1;
    legend(ax.Children(legendIndx),dataLegend_Unique, 'location', 'best', 'interpreter','none')        
end
ax.FontSize = 14;
return
