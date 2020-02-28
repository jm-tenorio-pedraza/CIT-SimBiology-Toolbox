function h = plotTrace(p_hat,varargin)
if nargin<2
    error('GWMCMC:toofewinputs','AMCMC requires atleast 2 inputs.')
end
p=inputParser;
p.addParameter('type','individual');
p.addParameter('steps',1:1:size(p_hat,1));
p.addParameter('names',{});
p.addParameter('ESS',[]);
p.addParameter('interpreter','none');

p.parse(varargin{:});
p=p.Results;

dim_phat=size(p_hat);

n_dim=length(dim_phat);
switch n_dim
    case 2
        if dim_phat(1)>dim_phat(2)
        else
            p_hat=p_hat';
        end
        n_p=dim_phat(2);

    case 3
        n_p=dim_phat(1);
        n_w=dim_phat(2);
        p.steps=1:dim_phat(3);
end
n_col=ceil(sqrt(n_p));
n_row=ceil(n_p/n_col);

% Limits of plots
max_y=max(max(p_hat(:,:)));
min_y=min(min(p_hat(:,:)));
% Plotting trace plots
switch n_dim
    case 2
if strcmp(p.type,'joint')
    for i=1:n_p
         plot(p.steps,p_hat(:,i));
         xlabel('MCMC step')
         text(p.steps(end), p_hat(end,i),p.names(i), 'interpreter', p.interpreter, 'FontSize', 8)
    end
elseif strcmp(p.type, 'individual')
    for i=1:n_p
        subplot(n_row, n_col,i)

         plot(p.steps,p_hat(:,i));
         try
         title(p.names(i), 'interpreter',p.interpreter)
         catch
         end
         xlabel('MCMC step')

         if ~isempty(p.ESS)
             legend(strjoin({'ESS=' num2str(p.ESS(i))},''),'Location','best')
         end

              ylim([min_y max_y])
%              
    end

end
    case 3
        for i=1:n_p
            subplot(n_row, n_col,i)
            hold on
%             if dim_phat(2)<=dim_phat(3)
            p_ij=reshape(p_hat(i,:,:),n_w,dim_phat(3),1);                  % NxT matrix with the numer of walkers/ensemble in the rows and the columns representing the steps 
                colors=linspecer(n_w);

%             else
%                 p_ij=reshape(p_hat(:,i,:),dim_phat(3),dim_phat(1),1);
%                 colors=linspecer(dim_phat(2));
% 
%             end
            h=plot(p.steps,p_ij);
            for j=1:length(h)
                set(h(j),'color',colors(j,:))
            end
            
             if ~isempty(p.ESS)
                legend(strjoin({'ESS=' num2str(p.ESS(i))},''),'Location','best')
            end
%             ylim([min_y max_y])

            try
            title(p.names(i), 'interpreter', p.interpreter)
            catch
            end
            xlabel('MCMC step')
            %ylim([min_y max_y])
        end
            
end

return
    