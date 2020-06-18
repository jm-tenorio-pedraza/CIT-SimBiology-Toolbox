function plotSLD(PI, treatments)
% Plot mean SLD
figure('Renderer', 'painters', 'Position', [10 10 800 700])
colors = linspecer(4);
lineStyle = {'-' '-.' '--' ':'};
ncol = ceil(sqrt(length(treatments)));
nrow = ceil(length(treatments)/ncol);
for i =1:length(treatments)
   subplot(nrow,ncol,i)
   hold on
   for j=1:4
       mean_ij=mean(PI.output(i).SLD(PI.output(i).Response(:,j),:),'omitnan');
       ci_ub = quantile(PI.output(i).SLD(PI.output(i).Response(:,j),:),.975);
       ci_lb = quantile(PI.output(i).SLD(PI.output(i).Response(:,j),:),.025);
       
       pi_ub = quantile(PI.output(i).SLD_sigma(PI.output(i).Response(:,j),:),.975);
       pi_lb = quantile(PI.output(i).SLD_sigma(PI.output(i).Response(:,j),:),.025);
       
    h=plot(PI.tspan/7, mean_ij);
    ci_plot = patch('XData', [PI.tspan/7 PI.tspan(end:-1:1)/7], 'YData', [ci_ub ci_lb(end:-1:1)]);
    pi_plot = patch('XData', [PI.tspan/7 PI.tspan(end:-1:1)/7], 'YData', [pi_ub pi_lb(end:-1:1)]);
    h.Color = colors(j,:);
    h.LineWidth = 2;
    ci_plot.FaceColor = colors(j,:);
    ci_plot.LineStyle = lineStyle{j};
    ci_plot.EdgeColor = colors(j,:);
    ci_plot.FaceAlpha = 0.2;
    
    pi_plot.LineStyle = 'none';
    pi_plot.FaceColor = colors(j,:);
    pi_plot.FaceAlpha = 0.01;
   end
   ax = gca;
   legend(ax.Children(end:-3:1),{'PD'; 'SD'; 'PR'; 'CR'})

 set(gca, 'YSCale', 'linear','YLim', [-100 100])
 xlabel('Time [weeks]')
ylabel('Change in tumor diameter wrt baseline')
title (strjoin({'Mean tumor response in virtual patients (' treatments{i} ')'},''),'interpreter', 'none')

end

