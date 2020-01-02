function plotIIVParams(postSamples, PI,varargin)

for j=1:length(PI.H.IndividualParams)
    figure
hold on

for i=1:length(PI.data)
    histogram(exp(postSamples(:,PI.H.IndividualParams(j).Index(i))),...
        'FaceColor',PI.data(i).colors,'FaceAlpha', 0.5, 'Normalization', 'probability')
    
end
legend({PI.data(:).Name},'interpreter', 'none')
ylabel('prob')
xlabel('Deviations wrt mean parameter')
title(strjoin({'Inter-individual variation in', PI.H.IndividualParams(j).name},' '))
end
% histogram(exp(postSamples(:,PI.H.IndividualParams(1).EtaIndex)),...
%     'FaceColor',PI.data(i).colors,'FaceAlpha', 0.5, 'Normalization', 'probability')

