function plotSimOutput(PI,colIndx, varargin)
par=inputParser;
par.addParameter('all', true)
par.addParameter('indiv', true)
par.addParameter('addErrorVar', true)
par.addParameter('newFig', true)
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

treatments=([PI.data(sim_indx).Group]);
if length(treatments)~=n_sim
    treatments = unique({PI.data(:).Group});
else
    treatments = unique(treatments);
end
if par.indiv
else
    treatments = {PI.data(:).Name};
end

treatment_colors=linspecer(length(treatments));
simIndx = find(sim_indx);

try
    minX = min(arrayfun(@(x) min(x.dataValue(:,colIndx)), PI.data, 'UniformOutput', true));
    maxX = max(arrayfun(@(x) max(x.dataValue(:,colIndx)), PI.data, 'UniformOutput', true));
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
    sim=plot(PI.data(simIndx(i)).simTime, PI.data(simIndx(i)).simOutput(:,colIndx));
    try
        dat = errorbar(PI.data(simIndx(i)).dataTime, PI.data(simIndx(i)).dataValue(:,colIndx), PI.data(simIndx(i)).SD(:,colIndx));
    catch
        dat=plot(PI.data(simIndx(i)).dataTime, PI.data(simIndx(i)).dataValue(:,colIndx));
        
    end
    if par.addErrorVar
        try
            X = [PI.data(simIndx(i)).simTime; PI.data(simIndx(i)).simTime(end:-1:1)];
            Y = [PI.data(simIndx(i)).ub(:,colIndx); PI.data(simIndx(i)).lb(end:-1:1,colIndx)];
            nanindx = isnan(Y);
            error = patch('XData', X(~nanindx),...
                'YData',Y(~nanindx));
            
        catch
            error = [];
        end
    else
        error = [];
    end
    title(PI.data(simIndx(i)).Name,'interpreter', 'none')
    if par.indiv
        col_i=treatment_colors(ismember(treatments,PI.data(simIndx(i)).Group),:);
    else
        col_i = treatment_colors(ismember(treatments,PI.data(simIndx(i)).Name),:);
    end
    sim.Color=col_i;
    sim.LineWidth= 2;
    dat.LineStyle='none';
    dat.Color = col_i;
    dat.MarkerFaceColor=col_i;
    dat.MarkerEdgeColor=col_i;
    
    dat.Marker='d';
    error.LineStyle = 'none';
    error.FaceColor = col_i;
    error.FaceAlpha = 0.2;
    
    ylabel(PI.variableUnits(colIndx))
    xlabel('Time [days]')
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
    
    ylim([floor(minX), ceil(maxX)])
    
    if strcmp(PI.variableUnits{colIndx}, 'Volume [ml]')
        ylim([1e-2, 3])
        
    elseif strcmp(PI.variableUnits{colIndx}, 'Percentage [%]')
        %ylim([1e-1, 30])
    else
    end
    
end
ax = gca;
if ~par.indiv
    if par.addErrorVar
        legend(ax.Children(end-2:-3:1),{PI.data(simIndx).Name}, 'location', 'best', 'interpreter','none')
    else
        legend(ax.Children(end-1:-2:1),{PI.data(simIndx).Name}, 'location', 'best', 'interpreter','none')

    end
        
end
return
