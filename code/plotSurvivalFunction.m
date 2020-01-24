function plotSurvivalFunction(PI, cutoff_point, groups)
colors=linspecer(length(groups));
cutoff_indx = find(PI.tspan == cutoff_point);
hold on
for i=1:length(groups)
    h = errorbar(PI.tspan(1:cutoff_indx),PI.output(i).S_t, sqrt(PI.output(i).V_t)*1.96);
    h.Color = colors(i,:);
end

legend(groups, 'interpreter', 'none')
title('Progression-free survival (6 months)')
ylabel('S(t)')
xlabel('Time [days]')
ylim([0 1])
