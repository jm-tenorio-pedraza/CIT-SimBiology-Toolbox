function plotBivariateMarginals_2(posterior, PI,varargin)
if nargin<2
    error('plotBivariateMarginals_2:toofewinputs','plotBivariateMarginals_2 requires atleast 2 inputs.')
end
param=inputParser;
param.addParameter('colors', linspecer(1))
param.addParameter('names', {PI.par(:).name})
param.parse(varargin{:})
param=param.Results;
p = size(posterior,2);
for i=1:p
    for j=1:p
        if i>j
            subplot(p,p,(i-1)*p+j)
           
            Z=gkde2(posterior(:,[i j]));
            contour(Z.x,Z.y,Z.pdf)

           
             xlabel(param.names(i),'Fontsize',5,...
                'Fontweight','bold')
             ylabel(param.names(j),'Fontsize',5,...
                'Fontweight','bold')
        elseif i==j
            subplot(p,p, (i-1)*p+j)
            histogram(posterior(:,i), 'FaceColor', param.colors, 'Normalization', 'probability')
            xlabel(param.names(i),'Fontsize',5,...
                'Fontweight','bold')
            ylabel('pdf','Fontsize',5,'Fontweight','bold')
        else
        end
        
    end
end
end
            