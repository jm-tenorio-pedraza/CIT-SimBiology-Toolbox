function phi=getPhi2(p,H,n_sim,varargin)
% Input:    1) parameter vector of dimension (px1)
%           2) H, structure array with k fields, one for each variant +1 ... 
%               for joint sigma parameters:
%               SigmaParams: Indexes of parameters that quantify variable
%                               measurement error standard deviation
%               Each k field contains:
%               PopulationParams: Indexes of population-wise parameters in
%                                   p
%               IndividualParams: Structure with fields that contain the...
%               indexes of parameters that vary at theindividual level
%               
% Output:   Matrix (n_simxp) of parameters
par=inputParser;
par.addParameter('initialValue',[])
par.parse(varargin{:})
par=par.Results;
if nargin<3
    error('getPhi requires 3 inputs')
end
if size(p,1)>size(p,2)
    p=p';
end
p_indiv=arrayfun(@(x)p(x.Index),H.IndividualParams,'UniformOutput',false);
p_indiv=cell2mat([p_indiv'])';
if size(p_indiv,1)>=size(p_indiv,2)
else
    p_indiv=p_indiv';
end
etaIndx = [H.IndividualParams.EtaIndex];

phi = repmat(p(H.PopulationParams),n_sim,1);
try
    phi(:,etaIndx) = p_indiv;
catch 
end
if ~isempty(par.initialValue)
    phi(:,end+1) = par.initialValue;
else
end