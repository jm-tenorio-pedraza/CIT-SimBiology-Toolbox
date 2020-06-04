function PI=getPosteriorPredictions(params,PI, output_handle,outputs,varargin)
if nargin<2
    error('GWMCMC:toofewinputs','AMCMC requires atleast 2 inputs.')
end
p=inputParser;
p.addParameter('simTime',PI.tspan);
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
OutputStruct = [];
dataOutput = repelem({'nan'}, size(params,1),1);
[OutputStruct(1:size(params,1)).dataOutput]=dataOutput{:,:} ;
for i=1:size(params,1)
    OutputStruct(i).dataOutput=output_handle(params(i,:));
end


for j=1:length(outputs)
    output_j=char(outputs(j));
    out_j=arrayfun(@(x)cellfun(@(y) y(:,j)',x.dataOutput, 'UniformOutput',false),OutputStruct,'UniformOutput',false);
    
    
    for k=1:length(PI.output)
        out_jk = cell2mat(cellfun(@(x)x{k,:},out_j, 'UniformOutput',false)');
        PI.output(k).(output_j) = out_jk;
    end
    
end
end