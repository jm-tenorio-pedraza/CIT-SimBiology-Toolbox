function [S] = shmPlot2(F,group,time, observables,varargin)
p = inputParser;
p.addParameter('tau', 0.01)
p.parse(varargin{:});
p=p.Results;
if size(time,1)>size(time,2)
else
    time=time';
end
if size(observables,1)>size(observables,2)
else
    observables=observables';
end

if size(group,1)>size(group,2)
else
    group=group';
end

observablesIndx = repmat(observables, length(group),1);
groupIndx = repelem(group, length(observables),1);
tau =p.tau*max(max(F));

n_u= sum(sum((F>tau)>0)>0);
S = [];
for k=1:n_u
    u_i = reshape(F(:,k), length(time),[]);
    plotindx = max(u_i)>tau;                            % Identify variables with F values above tau
    nrows = sum(plotindx);
    F_t_k = u_i(:,plotindx);
    
    a_m_k = max(F_t_k) - min(F_t_k);
    
    [a_m_k,I] = sort(a_m_k,'descend');
    F_t_k = F_t_k(:,I);
    group_k = groupIndx(plotindx);
    group_k = group_k(I);
    observables_k = observablesIndx(plotindx);
    observables_k = observables_k(I);
        
    figure('Renderer', 'painters', 'Position', [10 10 1500 600])
    for i=1:sum(plotindx)
        subplot(nrows,1,i)
        imagesc(F_t_k(:,i)'/a_m_k(i))
        ax=gca;
        ax.XTick = linspace(1,length(time),10);
        ax.XTickLabel = cellfun(@(x)num2str(x),num2cell(linspace(min(time),max(time),10)),'UniformOutput',false);
        ax.YTickLabel ={};
        colormap(jet)
        ylab=ylabel(strjoin({group_k{i,1}, observables_k{i,1}},'_'),'Interpreter', 'none');
        ylab.Rotation = 0;
        ylab.HorizontalAlignment='right';
%         text(max(time), max(F_t_k(:,i)/a_m_k(i)), num2str(a_m_k(i)))
    end
    S(k).F_t = F_t_k;
    S(k).a_m = a_m_k;
    
    
end
return
