function phi=getPhi3(p,H,n_sim)
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

if nargin<3
    error('getPhi requires 3 inputs')
end
if size(p,1)>size(p,2)
    p=p';
end
p_indiv=arrayfun(@(x)p(x.Index),H.IndividualParams,'UniformOutput',false);
p_indiv=cellfun(@(x) randn(n_sim)*p(x(2))+p(x(1)), p_indiv);
p_indiv=cell2mat(p_indiv);
try
    phi=[repmat(p(H.PopulationParams),n_sim,1),p_indiv'];
catch 
    phi=[repmat(p(H.FixedParams),n_sim,1),p_indiv'];
end