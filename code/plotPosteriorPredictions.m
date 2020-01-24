function plotPosteriorPredictions(PI,outputs,varargin)
if nargin<2
    error('GWMCMC:toofewinputs','AMCMC requires atleast 2 inputs.')
end
p=inputParser;
p.addParameter('simTime',PI.tspan);
p.addParameter('outputs','indiv');

p.parse(varargin{:});
p=p.Results;
n_sim=size(PI.data,1);
treatments=([PI.CI(:).Group]);

if length(treatments)~=n_sim
    treatments = unique({PI.data(:).Group});
else
    treatments = unique(treatments);
end
treatment_colors=linspecer(length(treatments));
for i=1:length(outputs)
    output_i=char(outputs(i));
    figure('Renderer', 'painters', 'Position', [10 10 1000 800])
    ncol=ceil(sqrt(length(PI.CI)));
    nrow=ceil(length(PI.CI)/ncol);

    for j=1:size(PI.CI,2)
        simTime=[p.simTime p.simTime(end:-1:1)];
        ci_data=[PI.CI(j).(output_i){'UB',:}, PI.CI(j).(output_i){'LB',:}(end:-1:1)];
        m_data =PI.CI(j).(output_i){'Median',:};
        ci_nan=or(isnan(ci_data),ci_data==0);
        m_nan=or(isnan(m_data), m_data==0);
       
        if strcmp(p.outputs,'indiv')
            subplot(nrow,ncol,j)
            pi_data=[PI.CI(j).(output_i){'Pred_UB',:}, PI.CI(j).(output_i){'Pred_LB',:}(end:-1:1)];
            pi_nan=or(isnan(pi_data), pi_data==0);
                        
             % Plotting prediction interval
            pi=patch('XData', simTime(~pi_nan),'YData',pi_data(~pi_nan));
        else
            if all(isnan( PI.data(j).dataValue(:,i)))
                continue
            end
        end
        hold on
        
        pi_data=[PI.CI(j).(output_i){'Pred_UB',:}, PI.CI(j).(output_i){'Pred_LB',:}(end:-1:1)];
            pi_nan=or(isnan(pi_data), pi_data==0);
                        
             % Plotting prediction interval
%             pi=patch('XData', simTime(~pi_nan),'YData',pi_data(~pi_nan));
        % Plotting credible interval
        ci=patch('XData', simTime(~ci_nan),'YData',ci_data(~ci_nan));

        % Plotting simulation median
        sim=plot(simTime(~m_nan),m_data(~m_nan));
        
        % Plotting data
        try
            dat= errorbar(PI.data(j).dataTime, PI.data(j).dataValue(:,i), PI.data(j).SD(:,i));
        catch
            dat=plot(PI.data(j).dataTime, PI.data(j).dataValue(:,i));

        end
       
        
        col_i=treatment_colors(ismember(treatments,PI.CI(j).Group),:);
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
        
        ylabel(PI.variableUnits(i),'interpreter', 'none')
        xlabel('Time [days]')
        grid on
        hold on
        if max(m_data)/min(m_data)>1e1
%             set(gca,'YScale','log')
        else
        end
        ax=gca;
            if strcmp(p.outputs,'indiv')
                title(PI.data(j).Name,'interpreter','none')
                try
                    legend(ax.Children, {'Data' output_i '95% Credible Interval' '95% Prediction Interval'},'interpreter', 'none')
                catch
                    legend(ax.Children, {output_i '95% Credible Interval' '95% Prediction Interval'},'interpreter', 'none')

                end
            end


    end
        ax = gca;
    if (strcmp('%',PI.variableUnits{i}))
        ylim([0 100])
    end
%     try
%         ylim(2.^([floor(log2(min(PI.data(j).dataValue(:,i)))) ceil(log2(max(PI.data(j).dataValue(:,i))))]))
%     catch
%     end

    
    if strcmp(p.outputs,'indiv')
    else
        title(output_i,'interpreter', 'none')
        legend(ax.Children(end:-4:3),[PI.data(:).Group],'Interpreter', 'none')
    end

end

end
