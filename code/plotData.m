function plotData(PI, stateVar, varargin)
p = inputParser;
p.addParameter('responseGrouping', false)
p.addParameter('kineticGrouping', false)
p.parse(varargin{:})
p=p.Results;
if  ~p.responseGrouping && ~(p.kineticGrouping)
    case_i = 1;
elseif p.responseGrouping && ~(p.kineticGrouping)
    case_i = 2;
else
    if ~p.responseGrouping && (p.kineticGrouping)
        case_i = 3;
    else
        case_i = 4;
    end
end

%% Plot
nvar = length(stateVar);
colors = table2cell(table(linspecer(length(PI.data))));
[PI.data(1:end).colors] = colors{:,:};
 figure('Position', [10 10 1000 900])
switch case_i
    case 1
            
        ncol = ceil(sqrt(nvar));
        nrow = ceil(sqrt(nvar/ncol));
        for i=1:length(stateVar)
            subplot(nrow, ncol, i)
            try
                 arrayfun(@(x)errorbar(x.dataTime, x.dataValue(:,i),x.SD(:,i),'Color',...
                    x.colors,'Marker','*'),PI.data)
            catch
                hold on
                arrayfun(@(x)plot(x.dataTime, x.dataValue(:,i),'Color',...
            x.colors,'Marker','*'),PI.data)
            end
            if i==length(stateVar)
        legend({PI.data(:).Name},'Interpreter', 'none')
            end
        end
        
    case {2 3 4}
        ncol = ceil(sqrt(nvar));
        nrow = ceil((nvar/ncol));
        for i=1:length(stateVar)
            subplot(nrow, ncol, i)
            hold on
             try
                 arrayfun(@(x)errorbar(x.dataTime, x.dataValue(:,i),x.SD(:,i),'Color',...
                    x.colors,'Marker','*'),PI.data)
            catch
                arrayfun(@(x)plot(x.dataTime, x.dataValue(:,i),'Color',...
            x.colors,'Marker','*'),PI.data)
            end
            if i==length(stateVar)
                legend({PI.data(:).Name},'Interpreter', 'none')
            end
        title(stateVar(i))
        end
        
end
       


return
