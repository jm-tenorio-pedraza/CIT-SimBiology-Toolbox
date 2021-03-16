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
p_indiv=arrayfun(@(x)p(x.EtaIndex)*p(x.Index),H.IndividualParams,...  % Obtain indexes for each individually-varying parameter
    'UniformOutput',false);
p_indiv=cell2mat(p_indiv')';                                              % Convert the ind indexes into a matrix
try
    p_cell=arrayfun(@(x)H.CellIndx*(p(x.Index)*p(x.EtaIndex))',...        % Obtain indexes for each individually-varying parameter
        H.CellParams,'UniformOutput',false); 
    p_cell=cell2mat(p_cell);                                                % Convert the ind indexes into a matrix
catch
    p_cell = [];
end

try
    p_resp=arrayfun(@(x)H.RespIndx*(p(x.Index)*p(x.EtaIndex))',...        % Obtain indexes for each individually-varying parameter
        H.RespParams,'UniformOutput',false); 
    p_resp=cell2mat(p_resp);                                                % Convert the ind indexes into a matrix
catch
    p_resp = [];
end


if size(p_indiv,1)<n_sim
        p_indiv=p_indiv';
end

if size(p_cell,1)<length(H.CellParams(1).Index)
        p_cell=p_cell';
end

if size(p_resp,1)<length(H.RespParams(1).Index)
        p_resp=p_resp';
end
% cellEtaIndx = [H.CellParams.EtaIndex];
% indivEtaIndx = [H.IndividualParams.EtaIndex];

phi = repmat(p(H.PopulationParams),n_sim,1);
try
    phi(:,[H.CellParams.EtaIndex]) = p_cell;
catch 
end
try
    phi(:,[H.IndividualParams.EtaIndex]) = p_indiv;
catch 
end
try
    phi(:,[H.RespParams.EtaIndex]) = p_resp;
catch 
end

if ~isempty(par.initialValue)
    phi(:,end+1:end+size(par.initialValue,2)) = par.initialValue;
else
end