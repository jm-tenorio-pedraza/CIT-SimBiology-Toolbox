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

observablesIndx = repmat(observables, length(group),1);             % Cell array with p rows, one for each group
groupIndx = repelem(group, length(observables),1);                  % Cell array with m rows, one for each observable
tau =p.tau*max(max(F));

n_u= sum(sum((F>tau)>0)>0);                                         % Determine how many elements per column of F(u,t) are greater than a certain threshold
S = [];
for k=1:n_u                         
    u_i = reshape(F(:,k)', length(time),[]);                         % Reshape the k column of F to a kxm, with one column for each observable and one row for each time point
    plotindx = max(u_i)>tau;                                        % Identify variables with F values above tau
    nrows = sum(plotindx);
    F_t_k = u_i(:,plotindx);                                        % Select columns of F to plot
    
    a_m_k = max(F_t_k) - min(F_t_k);                                % determine the a_m_ks for each column of F
    
    [a_m_k,I] = sort(a_m_k,'descend');                              % Sort the a_m_k to plot the largest one first
    F_t_k = F_t_k(:,I);                                             % Re-order F according to a
    group_k = groupIndx(plotindx);                                  % Select the appropriate indexes of the groups
    group_k = group_k(I);                                           % Reorder the groups according to a
    observables_k = observablesIndx(plotindx);                      % Select the appropriate indexes of the observables
    observables_k = observables_k(I);                               % Reorder the observables according to a
        
    figure('Renderer', 'painters', 'Position', [10 10 1500 600])
    if sum(plotindx)>20
        ncol=2;
        nrows=ceil(nrows/ncol);
        indx = [1:2:sum(plotindx) 2:2:sum(plotindx)];
    else
        ncol = 1;
        indx = 1:sum(plotindx);
    end
    
    for i=1:sum(plotindx)
        
        subplot(nrows,ncol,indx(i))
        imagesc(F_t_k(:,i)'/a_m_k(i))                               % Heat map of F' rescaled by a_m_k
        ax=gca;
        ax.XTickLabel = {};
        ax.YTickLabel ={};
        colormap(jet)
        ylab=ylabel(strjoin({group_k{i,1}, observables_k{i,1}},'_'),'Interpreter', 'none');
        ylab.Rotation = 0;
        ylab.HorizontalAlignment='right';
        ylab.FontSize = 11-0.1*sum(plotindx)/ncol;
%         text(max(time), max(F_t_k(:,i)/a_m_k(i)), num2str(a_m_k(i)))
        if i==1||i==nrows+1
            if k==1
            title(strjoin({'Sensitivity Heat Map of the ' num2str(k) '^{st} time series' },''))
            elseif k==2
                title(strjoin({'Sensitivity Heat Map of the ' num2str(k) '^{nd} time series' },''))
            else
                if k==3
                    title(strjoin({'Sensitivity Heat Map of the ' num2str(k) '^{rd} time series' },''))
                else
                    title(strjoin({'Sensitivity Heat Map of the ' num2str(k) '^{th} time series' },''))
                                                    
                end

            end
        end
        if mod(i,nrows)==0 || mod(i,sum(plotindx))==0
            ax.XTick = linspace(1,length(time),10);
            ax.XTickLabel = cellfun(@(x)num2str(x),num2cell(linspace(min(time),max(time),10)),'UniformOutput',false);
        end
    end
    S(k).F_t = F_t_k;
    S(k).a_m = a_m_k;
    
end
return
