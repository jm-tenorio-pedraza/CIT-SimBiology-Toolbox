function plotIndividualPosteriorPredictions(PI,outputs,varargin)
if nargin<2
    error('GWMCMC:toofewinputs','AMCMC requires atleast 2 inputs.')
end
p=inputParser;
p.addParameter('simTime',0:1:100);
p.parse(varargin{:});
p=p.Results;

treatments=unique({PI.CI.Group});
treatment_colors=linspecer(length(treatments));

for i=1:length(outputs)
    output_i=char(outputs(i));
    figure
    ncol=ceil(sqrt(length(PI.CI)));
    nrow=ceil(length(PI.CI)/ncol);
    
    for j=1:length(PI.CI)
        subplot(nrow,ncol,j)
        hold on
        ci_data=[PI.CI(j).(output_i){'UB',:}, PI.CI(j).(output_i){'LB',:}(end:-1:1)];
        pi_data=[PI.CI(j).(output_i){'Pred_UB',:}, PI.CI(j).(output_i){'Pred_LB',:}(end:-1:1)];
        m_data =PI.CI(j).(output_i){'Median',:};
        pi_nan=or(isnan(pi_data), pi_data==0);
        ci_nan=or(isnan(ci_data),ci_data==0);
        m_nan=or(isnan(m_data), m_data==0);
        simTime=[p.simTime p.simTime(end:-1:1)];
        % Plotting prediction interval
        pi=patch('XData', simTime(~pi_nan),'YData',pi_data(~pi_nan));

        % Plotting credible interval
        ci=patch('XData', simTime(~ci_nan),'YData',ci_data(~ci_nan));

        % Plotting simulation median
        sim=plot(simTime(~m_nan),m_data(~m_nan));
        
        % Plotting data
        dat=plot(PI.data(j).dataTime, PI.data(j).dataValue(:,i));
        
       
        
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
        dat.MarkerFaceColor=col_i;
        dat.MarkerEdgeColor=col_i;
        dat.Marker='d';
        
        title(PI.CI(j).Name,'interpreter','none')
        ylabel(PI.variableUnits(i))
        xlabel('Time [days]')
        grid on
        hold on
        if max(m_data)/min(m_data)>100
            set(gca,'YScale','log')
        else
        end
    end
end
