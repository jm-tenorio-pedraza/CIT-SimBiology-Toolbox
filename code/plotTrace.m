function plotTrace(p_hat,PI,varargin)
if nargin<2
    error('GWMCMC:toofewinputs','AMCMC requires atleast 2 inputs.')
end
p=inputParser;
p.addParameter('type','individual');
p.addParameter('steps',1:1:size(p_hat,1));
p.addParameter('names',{PI.par(:).name});
p.addParameter('ESS',[]);
p.parse(varargin{:});
p=p.Results;
if size(p_hat,1) >size(p_hat,2)
else
    p_hat=p_hat';
end


n_p=size(p_hat,2);
n_col=ceil(sqrt(n_p));
n_row=ceil(n_p/n_col);

% Limits of plots
% max_y=max(max(p_hat));
% min_y=min(min(p_hat));
% Plotting trace plots
if strcmp(p.type,'joint')


    for i=1:n_p

         plot(p.steps,p_hat(:,i));
         text(p.steps(end), p_hat(end,i),p.names(i), 'interpreter', 'none', 'FontSize', 8)
    end
else
    for i=1:n_p
        subplot(n_row, n_col,i)

         plot(p.steps,p_hat(:,i));
         title(p.names(i), 'interpreter', 'none')
         if ~isempty(p.ESS)
             legend(strjoin({'ESS=' num2str(p.ESS(i))},''),'Location','best')
         end
%          if any(p_hat(:,i)<0)
%              ylim([min_y max_y])
%          else
%              ylim([min_y max_y])
%              
%          end
    end
end
xlabel('MCMC step')
end
    