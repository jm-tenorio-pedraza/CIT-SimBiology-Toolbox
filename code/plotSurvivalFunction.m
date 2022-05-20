function plotSurvivalFunction(PI, cutoff_point, groups,varargin)
inputs=inputParser;
inputs.addParameter('survivalType','PFS')
inputs.parse(varargin{:});
inputs=inputs.Results;
colors=linspecer(length(groups));
cutoff_indx = find(PI.tspan >= cutoff_point,1);
figure
hold on
for i=1:length(groups)
    if strcmp(inputs.survivalType,'PFS')
        title(strjoin({'Progression-free survival (' num2str(cutoff_point/30) ' months)'}, ''))
        S_t=PI.output(i).S_t(1:cutoff_indx);
        V_t=PI.output(i).V_t(1:cutoff_indx);

        nanIndx = isnan(S_t);
        S_t = S_t(~nanIndx);
        ci_ub = S_t(~nanIndx)+sqrt(V_t(~nanIndx))*1.96;
        ci_lb = S_t(~nanIndx)-sqrt(V_t(~nanIndx))*1.96;
        tspan= PI.tspan(1:cutoff_indx);
        tspan=tspan(~nanIndx);
        
    else
                title(strjoin({'Overall survival (' num2str(cutoff_point/30) ' months)'}, ''))

        S_t=PI.output(i).OS_S_t(1:cutoff_indx);
        V_t=PI.output(i).OS_V_t(1:cutoff_indx);
        nanIndx = isnan(S_t);
        S_t = S_t(~nanIndx);

        ci_ub = S_t(~nanIndx)+sqrt(V_t(~nanIndx))*1.96;
        ci_lb = S_t(~nanIndx)-sqrt(V_t(~nanIndx))*1.96;
         tspan= PI.tspan(1:cutoff_indx);
        tspan=tspan(~nanIndx);
        
    end
    
    h = plot(tspan,S_t);
    ci_plot = patch('XData', [tspan tspan(end:-1:1)],...
        'YData', [ci_ub ci_lb(end:-1:1)]);
    
    h.Color = colors(i,:);
    h.LineWidth=2;
   
    ci_plot.FaceColor = colors(i,:);
    ci_plot.EdgeColor = colors(i,:);
    ci_plot.FaceAlpha = 0.5;
    
end
ax = gca;
legend(ax.Children(end:-2:1),groups, 'interpreter', 'none', 'FontSize', 12)

ylabel('S(t)')
xlabel('Time [days]')
ylim([0 1])
