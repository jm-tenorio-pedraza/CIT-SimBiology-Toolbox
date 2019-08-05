function plotTrace(p_hat,PI,varargin)
if nargin<2
    error('GWMCMC:toofewinputs','AMCMC requires atleast 2 inputs.')
end
p=inputParser;
p.addParameter('type','individual');
p.addParameter('steps',1e5:40:size(p_hat,1)*40+1e5-1);
p.addParameter('names',{PI.par(:).name});
p.parse(varargin{:});
p=p.Results;

n_p=size(p_hat,2);
n_r=size(p_hat,1);
n_col=ceil(sqrt(n_p));
n_row=ceil(n_p/n_col);
if strcmp(p.type,'joint')
% [ax_handles,~] = getNormFigure(1, 1,[],'figureFormat','landscape');
% set(gcf, 'CurrentAxes', ax_handles(1))

    for i=1:n_p

         plot(p.steps,p_hat(:,i));
         text(p.steps(end), p_hat(end,i),p.names(i), 'interpreter', 'none', 'FontSize', 8)
    end
else
%     [ax_handles,~] = getNormFigure(n_row, n_col,[],'figureFormat','landscape');
    for i=1:n_p
        subplot(n_row, n_col,i)

%          set(gcf, 'CurrentAxes', ax_handles(i))
         plot(p.steps,p_hat(:,i));
         title(p.names(i), 'interpreter', 'none')
    end
end
xlabel('MCMC step after burn-in')
end
    