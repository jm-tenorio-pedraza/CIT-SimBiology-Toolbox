function plotIIVParams(postSamples, PI,varargin)
try
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
catch
    close gcf
    

    
end
try
     for j=1:length(PI.H.CellParams)
        figure
        hold on
        
        for i=1:length(PI.H.CellParams(j).Index)
            histogram((postSamples(:,PI.H.CellParams(j).Index(i))),...
                'FaceColor',PI.data(i).colors,'FaceAlpha', 0.5, 'Normalization',...
                'probability')
            
        end
        legend({PI.data(:).Name},'interpreter', 'none')
        ylabel('prob')
        xlabel('Deviations wrt mean parameter')
        title(strjoin({'Inter-model variation in', PI.H.CellParams(j).name},' '))
       % set(gca, 'XScale', 'log')
     end
catch
        close gcf

end
return
% histogram(exp(postSamples(:,PI.H.IndividualParams(1).EtaIndex)),...
%     'FaceColor',PI.data(i).colors,'FaceAlpha', 0.5, 'Normalization', 'probability')

