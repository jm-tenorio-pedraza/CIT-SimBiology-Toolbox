function plotCorrMat(p_hat, Names, varargin)
if nargin<2
    error('GWMCMC:toofewinputs','AMCMC requires atleast 2 inputs.')
end
p=inputParser;
p.addParameter('model','MOC1');
p.addParameter('type','spearman');
p.parse(varargin{:});
p=p.Results;

n_col=size(p_hat,2);
n_row=size(p_hat,1);
if n_col==n_row
    p_hat_corr=p_hat;
else
p_hat_corr = corr(p_hat, 'type', p.type);
end
n_var = size(p_hat_corr,1);
imagesc(p_hat_corr, [-1 1])
set(gca, 'XTick', 1:n_var,'YTick', 1:n_var); % center x-axis and y-axis ticks on bins
set(gca, 'XTickLabel',Names, 'TickLabelInterpreter', 'none','XTickLabelRotation', 45); % set x-axis labels
set(gca, 'FontSize', 18); 

set(gca, 'YTickLabel', Names, 'TickLabelInterpreter', 'none'); % set y-axis labels
title(strjoin({'Correlation matrix for', p.model},' ' ), 'FontSize', 20); % set title

% Defining two-color gradient color map
blue1 = [0, 0, 1];
white=[1, 1, 1];
red = [1, 0, 0];
length=50;
colors_p1 = [linspace(blue1(1),white(1),length)', linspace(blue1(2),white(2),length)', linspace(blue1(3),white(3),length)'];
colors_p2 = [linspace(white(1),red(1),length)', linspace(white(2),red(2),length)', linspace(white(3),red(3),length)'];

colormap([colors_p1; colors_p2]);
colorbar;
fontSize=200/n_col;
for i=1:n_var
    for j=1:n_var
        if j>i
            if abs(p_hat_corr(i,j))>0.5
            text(i-0.25,j,num2str(round(p_hat_corr(i,j),2)), 'FontSize', fontSize*(abs(p_hat_corr(i,j))))
            else
            end
        else
        end
    end
end

