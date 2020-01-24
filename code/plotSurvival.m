function plotSurvival(T, censor, theta, groups)
N = size(theta,1);
figure()
colors=linspecer(length(groups));
hold on
for i=1:length(groups)
    indx = (i-1)*N+1;
    [f,x] = ecdf(T(indx:indx+N-1,1),'Censoring',censor(indx:indx+N-1),'function','survivor');
    h=stairs(x,f,'--','LineWidth', 2);
    h.Color = colors(i,:);
end

legend(groups, 'interpreter', 'none')
title('Progression-free survival (6 months)')
ylabel('S(t)')
xlabel('Time [days]')
ylim([0 1])
