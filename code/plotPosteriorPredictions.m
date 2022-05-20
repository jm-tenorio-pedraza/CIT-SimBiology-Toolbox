function plotPosteriorPredictions(PI,colIndx,varargin)
if nargin<2
    error('GWMCMC:toofewinputs','AMCMC requires atleast 2 inputs.')
end
p=inputParser;
p.addParameter('simTime',PI.tspan);
p.addParameter('outputs','indiv');
p.addParameter('plot', 'data');
p.addParameter('indx', 1:length(PI.data));
p.addParameter('central', 'median')
p.addParameter('newFig', true)
p.addParameter('TimeUnit', 'days')
p.addParameter('interpreter', 'tex');
p.addParameter('color', 'group');
p.addParameter('YScale', 'linear');
p.addParameter('indivData', false) 

p.parse(varargin{:});
p=p.Results;

% Plot all simulations or just those with data?
if strcmp(p.plot, 'all')
    n_sim=size(PI.data,1);
    sim_indx = ones(n_sim,1);
elseif strcmp(p.plot,'data')
    sim_indx=(arrayfun(@(x) ~all(isnan(x.dataValue(:,colIndx))), PI.data));
    n_sim = sum(sim_indx);
else
    sim_indx=(arrayfun(@(x) ~all(isnan(x.dataValue(:,colIndx))), PI.data));
    treatment_unique = unique([PI.data(:).Group]);
    n_sim = length(treatment_unique);
end

% Determine number of rows and columns
n_col=ceil(sqrt(n_sim));
n_row=ceil(n_sim/n_col);


% Determine colors
if strcmp(p.color, 'dataset')
    treatments = {PI.data(:).Name};
    treatment_colors=linspecer(length(treatments));
    
elseif strcmp(p.color, 'cell')
    cellIndx = 1:length(PI.H.CellTypes);
    treatment_colors= PI.H.CellIndx*linspecer(length(cellIndx));
else
    if strcmp(p.color, 'outcome')
        outcomes = {PI.data(:).Response};
        outcome_unique = unique(outcomes);
        
        ResponseM = zeros(length(PI.data),length(outcome_unique));
        for i=1:length(PI.data)
            ResponseM(i,:) = ismember(outcome_unique,outcomes(i));
        end
        treatment_colors = ResponseM*linspecer(length(outcome_unique));
    end
end

% Determine valid indexes
simIndx = find(sim_indx);

if p.newFig
    figure('Position', [10 10 900 1000])
end

% Identify which output to plot
outputs = PI.observablesFields;
output_i=char(outputs(colIndx));

% Determine consistent data ranges for all simulations
try
    minX = min(arrayfun(@(x) min(x.dataValue(:,colIndx)), PI.data, 'UniformOutput', true));
    maxX = max(arrayfun(@(x) max(x.dataValue(:,colIndx)), PI.data, 'UniformOutput', true));
    
catch
end

if strcmp(p.outputs, 'outcome')
    groups = unique([PI.data(simIndx).Group]);
    groupIndx = zeros(length(simIndx), length(groups));
    for i=1:length(simIndx)
        groupIndx(i,:) = ismember(groups, PI.data(simIndx(i)).Group);
    end
    
    for j=1:length(groups)
        dataIndx = find(groupIndx(:,j));
        n_dataIndx = length(dataIndx);
        simTime=[p.simTime p.simTime(end:-1:1)];                                % Extract the simulation time
        subplot(n_row,n_col,j)
        hold on
        for k=1:n_dataIndx
            k_indx = simIndx(dataIndx(k));
            ci_data=[PI.CI(k_indx).(output_i){'UB',:},...                       % Extract the CIs
                PI.CI(k_indx).(output_i){'LB',:}(end:-1:1)];
            if strcmp(p.central, 'median')                                          % Extract the central tendency value desired
                m_data =PI.CI(k_indx).(output_i){'Median',:};
            elseif strcmp(p.central, 'mean')
                m_data =PI.CI(k_indx).(output_i){'Mean',:};
            end
            ci_nan=or(isnan(ci_data),ci_data==0);                                   % Determine whether nan's are observed to eliminate them
            m_nan=or(isnan(m_data), m_data==0);
            pi_data=[PI.CI(k_indx).(output_i){'Pred_UB',:},...
                PI.CI(k_indx).(output_i){'Pred_LB',:}(end:-1:1)];
            pi_nan=or(isnan(pi_data), pi_data==0);
            pi=patch('XData', simTime(~pi_nan),'YData',pi_data(~pi_nan));       % Plot PIs
            
            ci=patch('XData', simTime(~ci_nan),'YData',ci_data(~ci_nan));           % Credible interval
            sim=plot(simTime(~m_nan),m_data(~m_nan));                               % Mean or median of simulations
            if p.indivData
                dataIndx = ismember({PI.IndivData(1:end).Name},PI.data(k_indx).Name);
                dat=arrayfun(@(x)plot(x.dataTime,x.dataValue(:,colIndx)),PI.IndivData(dataIndx),'UniformOutput',false);
            else
                try
                    dat= errorbar(PI.data(k_indx).dataTime,...                      % Plot data point + SD error bar
                        PI.data(k_indx).dataValue(:,colIndx),...
                        PI.data(k_indx).SD(:,colIndx));
                catch
                    dat=plot(PI.data(k_indx).dataTime, ...                          % Plot data point only
                        PI.data(k_indx).dataValue(:,colIndx));
                end
            end
            % Line specs
            col_i=treatment_colors(k_indx,:);
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
            if p.indivData
                for i=1:length(dat)
                    dat{i}.Color =col_i;
                    dat{i}.Marker='d';
                    dat{i}.MarkerFaceColor=col_i;
                    dat{i}.MarkerEdgeColor=col_i;
                end
            else
                dat.LineStyle='none';
                dat.Color = col_i;
                dat.MarkerFaceColor=col_i;
                dat.MarkerEdgeColor=col_i;
                dat.Marker='d';
            end
           
        end
        % Axes specs
            ylabel(PI.variableUnits(colIndx),'interpreter', 'none')
            xlabel(strjoin({'Time [' p.TimeUnit ']'},''))
            grid on
            ax=gca;
            title(groups(j),'interpreter',p.interpreter)
%         if maxX<1
%             ylim([(minX), (maxX)])
%         else
%             ylim([floor(minX), ceil(maxX)])
%         end
        set(ax, 'YScale', p.YScale)  
    end

else
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
        
        if ismember(p.outputs,{'indiv'})                                        % Extract prediction intervals only if individual outputs are desired
            subplot(n_row,n_col,j)
            pi_data=[PI.CI(simIndx(j)).(output_i){'Pred_UB',:},...
                PI.CI(simIndx(j)).(output_i){'Pred_LB',:}(end:-1:1)];
            pi_nan=or(isnan(pi_data), pi_data==0);
            pi=patch('XData', simTime(~pi_nan),'YData',pi_data(~pi_nan));       % Plot PIs
        elseif strcmp(p.outputs, 'group')
            if all(isnan( PI.data(simIndx(j)).dataValue(:,colIndx)))
                continue
            else
                pi_data=[PI.CI(simIndx(j)).(output_i){'Pred_UB',:},...
                    PI.CI(simIndx(j)).(output_i){'Pred_LB',:}(end:-1:1)];
                pi_nan=or(isnan(pi_data), pi_data==0);
                
                %Plotting prediction interval
%                 pi=patch('XData', simTime(~pi_nan),'YData',pi_data(~pi_nan));
            end
        else
            
        end
        hold on
        ci=patch('XData', simTime(~ci_nan),'YData',ci_data(~ci_nan));           % Credible interval
        sim=plot(simTime(~m_nan),m_data(~m_nan));                               % Mean or median of simulations
        if p.indivData
            dataIndx = ismember({PI.IndivData(1:end).Name},PI.data(simIndx(j)).Name);
            dat=arrayfun(@(x)plot(x.dataTime,x.dataValue(:,colIndx)),PI.IndivData(dataIndx),'UniformOutput',false);
        else
            try
                dat= errorbar(PI.data(simIndx(j)).dataTime,...                      % Plot data point + SD error bar
                    PI.data(simIndx(j)).dataValue(:,colIndx),...
                    PI.data(simIndx(j)).SD(:,colIndx));
            catch
                dat=plot(PI.data(simIndx(j)).dataTime, ...                          % Plot data point only
                    PI.data(simIndx(j)).dataValue(:,colIndx));
            end
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
        
        if p.indivData
            for i=1:length(dat)
                dat{i}.Color =col_i;
                dat{i}.Marker='d';
                dat{i}.MarkerFaceColor=col_i;
                dat{i}.MarkerEdgeColor=col_i;
            end
        else
            dat.LineStyle='none';
            dat.Color = col_i;
            dat.MarkerFaceColor=col_i;
            dat.MarkerEdgeColor=col_i;
            dat.Marker='d';
        end
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
%         if maxX<1
%             ylim([(minX), (maxX)])
%         else
%             ylim([floor(minX), ceil(maxX)])
%         end
        set(ax, 'YScale', p.YScale)
        
    end
    
    
    ax = gca;                                                                   % Add title and legends
    
    legends = {PI.data(sim_indx).Name};
    if iscell(legends{1})
        legends = [PI.data(sim_indx).Name];
    end
    
    
    
    
    if strcmp(p.outputs,'indiv')
        %legend(ax.Children(end:-4:4),legends,'Interpreter', 'none','location','best')
    else
        title(PI.observablesPlot(colIndx),'interpreter', p.interpreter)
        legend(ax.Children(end:-3:3),legends,'Interpreter', 'none', 'location','best')
    end
    
    ax.FontSize = 14;
end

return