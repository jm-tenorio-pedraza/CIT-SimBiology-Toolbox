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
figure
for i=1:p
    for j=1:p
        if i>j
            subplot(p,p,(i-1)*p+j)
           if license('test','statistics_toolbox')
               [z,xi]=ksdensity(posterior(:,[i j]));
               z = reshape(z,ceil(sqrt(900)),[]);
               x = reshape(xi(:,2),ceil(sqrt(900)),[]);
               y =reshape(xi(:,1),ceil(sqrt(900)),[]);
               contour(x,y,z)
           else
            Z=gkde2(posterior(:,[i j]));
            contour(Z.x,Z.y,Z.pdf)

           end

%            coor_ij=corr(posterior(:,[i j]));
%            text(
           if j==1
             ylab=ylabel(param.names(i),'Fontsize',10,...
                'Fontweight','normal','Interpreter','none');
            ylab.Rotation=90;
           elseif i==p
             xlab=xlabel(param.names(j),'Fontsize',10,...
                'Fontweight','normal','Interpreter','none');
            xlab.Rotation=0;
           end
             if j==1&&i==p
             xlab=xlabel(param.names(j),'Fontsize',10,...
                 'Fontweight','normal','Interpreter','none');
             xlab.Rotation=0;
             
             end
        elseif i==j
           
            subplot(p,p, (i-1)*p+j)
            h=histogram(posterior(:,i), 'FaceColor', param.colors, 'Normalization', 'probability');
            p_max=max(h.Values);
            m=mean(posterior(:,i));
            sigma=std(posterior(:,i));
            if p<11
%             text(m,p_max,{strjoin({'\mu = ' num2str(m)},'') strjoin({'\sigma = ' num2str(sigma)},'') })
%             legend({strjoin({'\mu = ' num2str(m)},'') strjoin({'\sigma = ' num2str(sigma)},'') },'Location','best')
            end
             if j==1||i==p
             xlab=xlabel(param.names(j),'Fontsize',10,...
                 'Fontweight','normal','Interpreter','none');
             xlab.Rotation=0;
             else
             end
             ylabel('pdf','Fontsize',5,'Fontweight','normal')

        else
        end
        
    end
end
tightfig;
end
            