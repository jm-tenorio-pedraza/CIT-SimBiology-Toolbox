function T=detectPD(x, tspan,varargin)
inputs=inputParser;
inputs.addParameter('type','PFS')
inputs.addParameter('minTV',0)
inputs.addParameter('CC',[])
inputs.parse(varargin{:})
inputs=inputs.Results;
if strcmp(inputs.type,'PFS')
    indx = find(and(x>20, x-inputs.minTV>.5),1);
else
    indx = find(x>=inputs.CC,1);
end
    
if isempty(indx)
    T = tspan(end);
else
    T = tspan(indx);
end
