function plotORR(PI, treatments, varargin)

inputs = inputParser;
inputs.addParameter('output', 'SLD')
inputs.addParameter('lines', 'mean')
inputs.addParameter('credibility', .9)
inputs.addParameter('prediction', .95)
inputs.parse(varargin{:});
inputs = inputs.Results;
outputs = inputs.output;

lb_cred=(1-inputs.credibility)/2;
ub_cred=inputs.credibility + (1-inputs.credibility)/2;

lb_pred=(1-inputs.prediction)/2;
ub_pred=inputs.prediction + (1-inputs.prediction)/2;
for k = 1:length(outputs)
    output_k = outputs{k};
    output_sigma_k = strjoin({output_k, 'sigma'}, '_');
    figure('Renderer', 'painters', 'Position', [10 10 800 700])
    colors = linspecer(4);
    lineStyle = {'-' '-.' '--' ':'};
    ncol = ceil(sqrt(length(treatments)));
    nrow = ceil(length(treatments)/ncol);
    
    minX = min(arrayfun(@(x)min(min(x.(output_k))), PI.output, 'UniformOutput',true));
    maxX = max(arrayfun(@(x)max(max(x.(output_k))), PI.output, 'UniformOutput',true));

    for i =1:length(treatments)
        subplot(nrow,ncol,i)
        hold on
        if strcmp(inputs.lines, 'mean')
        for j=1:4
            mean_ij=mean(PI.output(i).(output_k)(PI.output(i).Response(:,j),:),1,'omitnan');
            ci_ub = quantile(PI.output(i).(output_k)(PI.output(i).Response(:,j),:),ub_cred);
            ci_lb = quantile(PI.output(i).(output_k)(PI.output(i).Response(:,j),:),lb_cred);
            try
            pi_ub = quantile(PI.output(i).(output_sigma_k)(PI.output(i).Response(:,j),:),ub_pred);
            pi_lb = quantile(PI.output(i).(output_sigma_k)(PI.output(i).Response(:,j),:),lb_pred);
                        pi_plot = patch('XData', [PI.tspan/7 PI.tspan(end:-1:1)/7], 'YData', [pi_ub pi_lb(end:-1:1)]);

            catch
                pi_plot =[];
            end
            
            h=plot(PI.tspan/7, mean_ij);
            ci_plot = patch('XData', [PI.tspan/7 PI.tspan(end:-1:1)/7], 'YData', [ci_ub ci_lb(end:-1:1)]);
            
            h.Color = colors(j,:);
            h.LineWidth = 2;
            ci_plot.FaceColor = colors(j,:);
            ci_plot.LineStyle = lineStyle{j};
            ci_plot.EdgeColor = colors(j,:);
            ci_plot.FaceAlpha = 0.2;
            
            pi_plot.LineStyle = 'none';
            pi_plot.FaceColor = colors(j,:);
            pi_plot.FaceAlpha = 0.2;
        end
        ax = gca;
        elseif strcmp(inputs.lines, 'individual')
            indx = zeros(4,1);
            h=plot(PI.tspan/7, PI.output(i).(output_k));
            colorsMat = PI.output(i).Response*colors;
            for j=1:length(h)
                h(j).Color = colorsMat(j,:);
            end
            for j=1:4
                try
                indx(j) = find(PI.output(i).Response(:,j),1);
                catch
                end
                indx = indx(indx~=0);
            end
            legend(h(indx),{'PD'; 'SD'; 'PR'; 'CR'})
        end
        if strcmp(output_k, 'SLD')
            set(gca, 'YSCale', 'linear','YLim', [-100 100])
        else
            set(gca, 'YScale', 'linear', 'YLim', [minX maxX])
        end
        plot(PI.tspan/7, repelem(20,1,length(PI.tspan)), '--k')
        plot(PI.tspan/7, repelem(-30,1,length(PI.tspan)), '--k')
        legend(ax.Children(end:-3:3),{'PD'; 'SD'; 'PR'; 'CR'})

        xlabel('Time [weeks]')
        ylabel('Change in tumor diameter wrt baseline')
        title (strjoin({'Mean tumor response in virtual patients (' treatments{i} ')'},''),'interpreter', 'none')
    end
end
return


