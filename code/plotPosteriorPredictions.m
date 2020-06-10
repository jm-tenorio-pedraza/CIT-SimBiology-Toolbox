function plotPosteriorPredictions(PI,colIndx,varargin)
if nargin<2
    error('GWMCMC:toofewinputs','AMCMC requires atleast 2 inputs.')
end
p=inputParser;
p.addParameter('simTime',PI.tspan);
p.addParameter('outputs','indiv');
p.addParameter('all', false);
p.addParameter('indx', 1:length(PI.data));
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
    cellIndx = 1:length(PI.H.CellTypes);
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
outputs = PI.observablesFields;
output_i=char(outputs(colIndx));

try
    minX = min(arrayfun(@(x) min(x.dataValue(:,colIndx)), PI.data, 'UniformOutput', true));
    maxX = max(arrayfun(@(x) max(x.dataValue(:,colIndx)), PI.data, 'UniformOutput', true));
    
catch
end
% For each valid simulation do:
for j=1:length(simIndx)
    simTime=[p.simTime p.simTime(end:-1:1)];                                % Extract the simulation time
    ci_data=[PI.CI(simIndx(j)).(output_i){'UB',:},...                       % Extract the CIs
        PI.CI(simIndx(j)).(output_i){'LB',:}(end:-1:1)];
    if strcmp(p.central, 'median')                                          % Extract the central tendency value desired
        m_data =PI.CI(simIndx(j)).(output_i){'Median',:};
    elseif strcmp(p.central, 'mean')
        m_data =PI.CI(simIndx(j)).(output_i){'Mean',:};
    end
    ci_nan=or(isnan(ci_data),ci_data==0);                                   % Determine whether nan's are observed to eliminate them 
    m_nan=or(isnan(m_data), m_data==0);
    
    if strcmp(p.outputs,'indiv')                                            % Extract prediction intervals only if individual outputs are desired
        subplot(n_row,n_col,j)
        pi_data=[PI.CI(simIndx(j)).(output_i){'Pred_UB',:},...
            PI.CI(simIndx(j)).(output_i){'Pred_LB',:}(end:-1:1)];
        pi_nan=or(isnan(pi_data), pi_data==0);
        pi=patch('XData', simTime(~pi_nan),'YData',pi_data(~pi_nan));       % Plot PIs
    else
        if all(isnan( PI.data(simIndx(j)).dataValue(:,colIndx)))            
            continue
        else
            pi_data=[PI.CI(simIndx(j)).(output_i){'Pred_UB',:},...
                PI.CI(simIndx(j)).(output_i){'Pred_LB',:}(end:-1:1)];
            pi_nan=or(isnan(pi_data), pi_data==0);
            
            %Plotting prediction interval
            pi=patch('XData', simTime(~pi_nan),'YData',pi_data(~pi_nan));      
        end      
    end
    hold on
    ci=patch('XData', simTime(~ci_nan),'YData',ci_data(~ci_nan));           % Credible interval
    sim=plot(simTime(~m_nan),m_data(~m_nan));                               % Mean or median of simulations
    try
        dat= errorbar(PI.data(simIndx(j)).dataTime,...                      % Plot data point + SD error bar
            PI.data(simIndx(j)).dataValue(:,colIndx),...
            PI.data(simIndx(j)).SD(:,colIndx));
    catch
        dat=plot(PI.data(simIndx(j)).dataTime, ...                          % Plot data point only
            PI.data(simIndx(j)).dataValue(:,colIndx));
    end
    %% Line specs
    col_i=treatment_colors(simIndx(j),:);
    pi.FaceColor=col_i;
    pi.EdgeColor=col_i;
    pi.FaceAlpha=0.2;
    pi.LineStyle='none';
    
    ci.FaceColor=col_i;
    ci.EdgeColor=col_i;
    ci.FaceAlpha=0.4;
    ci.LineStyle='--';
    ci.LineWidth=1;
    
    sim.Color=col_i;
    
    dat.LineStyle='none';
    dat.Color = col_i;
    dat.MarkerFaceColor=col_i;
    dat.MarkerEdgeColor=col_i;
    dat.Marker='d';
    %% Axes specs
    ylabel(PI.variableUnits(colIndx),'interpreter', 'none')
    xlabel(strjoin({'Time [' p.TimeUnit ']'},''))
    grid on
    hold on

    ax=gca;
    if strcmp(p.outputs,'indiv')
        title(PI.data(simIndx(j)).Name,'interpreter','none')
%         try
%             legend(ax.Children, {'Data' PI.observablesPlot{colIndx}...
%                 '95% Credible Interval' '95% Prediction Interval'},...
%             'interpreter', p.interpreter,'Location', 'best')
%         catch
%             legend(ax.Children, {PI.observablesPlot{colIndx}...
%                 '95% Credible Interval' '95% Prediction Interval'},...
%                 'interpreter', p.interpreter,'Location', 'best')
%         end
    end
    if maxX<1
        ylim([(minX), (maxX)])
        %set(gca,'YScale','log')
    else
        ylim([floor(minX), ceil(maxX)])
    end
    %set(ax, 'YScale', 'log')
    if max(m_data)/min(m_data)>1e1
       % set(gca,'YScale','log')
    end
end


ax = gca;                                                                   % Add title and legends

legends = [PI.data(sim_indx).Name];
if length(legends)~= length(PI.data)
    legends = {PI.data(sim_indx).Name};
end

if strcmp(p.outputs,'indiv')
    %legend(ax.Children(end:-4:4),legends,'Interpreter', 'none','location','best')
else
    title(PI.observablesPlot(colIndx),'interpreter', p.interpreter)
    legend(ax.Children(end:-3:3),legends,'Interpreter', 'none', 'location','best')
end

return