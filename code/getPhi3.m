function phi=getPhi3(p,H,n_sim,varargin)
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
p_indiv=arrayfun(@(x)[p(x.EtaIndex) p(x.EtaIndex)*p(x.Index)],H.IndividualParams,'UniformOutput',false);  % Obtain indexes for each individually-varying parameter
p_indiv=cell2mat(p_indiv')';                                              % Convert the ind indexes into a matrix
try
p_cell=arrayfun(@(x)H.CellIndx*([p(x.EtaIndex) p(x.EtaIndex)*p(x.Index)])',H.CellParams,'UniformOutput',false);         % Obtain indexes for each individually-varying parameter
p_cell=cell2mat(p_cell);                                                % Convert the ind indexes into a matrix
catch
    p_cell = [];
end

if size(p_indiv,1)<n_sim
        p_indiv=p_indiv';
end

if size(p_cell,1)<length(H.CellParams(1).Index)
        p_cell=p_cell';
end

cellEtaIndx = [H.CellParams.EtaIndex];
indivEtaIndx = [H.IndividualParams.EtaIndex];

phi = repmat(p(H.PopulationParams),n_sim,1);
try
    phi(:,cellEtaIndx) = p_cell;
catch 
end
try
    phi(:,indivEtaIndx) = p_indiv;
catch 
end
if ~isempty(par.initialValue)
    phi(:,end+1:end+size(par.initialValue,2)) = par.initialValue;
else
end