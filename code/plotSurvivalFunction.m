function plotSurvivalFunction(PI, cutoff_point, groups)
colors=linspecer(length(groups));
cutoff_indx = find(PI.tspan == cutoff_point);
hold on
for i=1:length(groups)
    ci_ub = PI.output(i).S_t+sqrt(PI.output(i).V_t)*1.96;
    ci_lb = PI.output(i).S_t-sqrt(PI.output(i).V_t)*1.96;
    
    h = plot(PI.tspan(1:cutoff_indx),PI.output(i).S_t);
    h.Color = colors(i,:);
    h.LineWidth=2;
    ci_plot = patch('XData', [PI.tspan(1:cutoff_indx) PI.tspan(cutoff_indx:-1:1)],...
        'YData', [ci_ub ci_lb(end:-1:1)]);
    ci_plot.FaceColor = colors(i,:);
    ci_plot.EdgeColor = colors(i,:);
    ci_plot.FaceAlpha = 0.5;
    
end
ax = gca;
legend(ax.Children(end:-2:1),groups, 'interpreter', 'none', 'FontSize', 12)
title('Progression-free survival (6 months)')
ylabel('S(t)')
xlabel('Time [days]')
ylim([0 1])
