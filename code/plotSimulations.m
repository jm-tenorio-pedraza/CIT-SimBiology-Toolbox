function plotSimulations(PI,varargin)
inputs=inputParser;
inputs.addParameter('YScale','log')
inputs.addParameter('YLabel',PI.variableUnits)

inputs.parse(varargin{:})
inputs=inputs.Results;
 figure
    nrow=ceil(sqrt(length(PI.observablesPlot)));
    ncol=ceil(length(PI.observablesPlot)/nrow);
    colors=linspecer(length(PI.data));
    for i=1:length(PI.observablesPlot)
        output_i=PI.observablesPlot{i};
        subplot(nrow,ncol,i)
        hold on
        for j=1:length(PI.data)
            %             plot(PI.tspan,mean(PI.output(j).(output_i),1,'omitnan'),'Color',colors(j,:));
            sim_output_ij= PI.output(j).(output_i);
            sim_CI_LB = quantile(sim_output_ij, .025,1);
            sim_CI_UB = quantile(sim_output_ij, .975,1);
            sim_median = median(sim_output_ij,1,'omitnan');
            h=plot(PI.tspan,sim_median,'Color',colors(j,:),'LineWidth',2);
            ci=patch('XData',[PI.tspan PI.tspan(end:-1:1)], ...
                'YData',[sim_CI_UB sim_CI_LB(end:-1:1)]);
            ci.FaceColor = colors(j,:);
            ci.EdgeColor = colors(j,:);
            ci.FaceAlpha = .1;
            ci.LineStyle = 'none';
            
        end
        title(output_i)
        set(gca,'YScale',inputs.YScale)
        ylabel(inputs.YLabel{i})
    end
    ax =gca;
    legend(ax.Children(end:-2:1),[PI.data(:).Group],'interpreter','none')
end