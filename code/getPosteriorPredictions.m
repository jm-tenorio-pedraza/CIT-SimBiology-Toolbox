function PI=getPosteriorPredictions(params,PI, output_handle,outputs,varargin)
if nargin<2
    error('GWMCMC:toofewinputs','AMCMC requires atleast 2 inputs.')
end
p=inputParser;
p.addParameter('simTime',0:1:100);
p.parse(varargin{:});
p=p.Results;

if size(params,3)>1
    params=params(:,:)';
end
% Set up outputs
PI.output=struct('Name',{PI.data(1:end).Name}');
out=repelem({NaN(size(params,1),length(p.simTime))},length(PI.output),1);
for i=1:length(outputs)
    output_i=char(outputs(i));
    [PI.output(1:end).(output_i)]=out{:,:};
end
% Obtaining simulations
for i=1:size(params,1)
    PI_i=output_handle(params(i,:));
    
    for j=1:length(outputs)
         output_j=char(outputs(j));
         out_j=arrayfun(@(x)x.simOutput(:,j)',PI_i.data,'UniformOutput',false);
         for k=1:length(PI.output)
             PI.output(k).(output_j)(i,:)=out_j{k,:};
         end
         
    end
end