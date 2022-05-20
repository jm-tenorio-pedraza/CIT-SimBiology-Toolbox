function plotSimValue(PI,varargin)
inputs=inputParser;
inputs.addParameter('scale','log')
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
            plot(PI.tspan,(PI.data(j).simValue(:,i)),'Color',colors(j,:));
        end
        title(output_i)
        set(gca,'YScale',inputs.scale)
    end
    ax =gca;
    legend(ax.Children,[PI.data(:).Group])
end