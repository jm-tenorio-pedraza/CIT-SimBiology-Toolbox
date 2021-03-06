function phi=getPhi(p,H,n_sim)
% Input:    1) parameter vector of dimension (px1)
%           2) H, structure array with k fields, one for each variant +1 ... 
%               for joint sigma parameters:
%               SigmaParams: Indexes of parameters that quantify variable
%                               measurement error standard deviation
%               Each k field contains:
%               PopulationParams: Indexes of population-wise parameters in
%                                   p
%               IndividualParams: Indexes of parameters that vary at the...
%                                   individual level
%               
% Output:   Matrix (n_simxp) of parameters
if nargin<3
    error('getPhi requires 3 inputs')
end
if size(p,1)>size(p,2)
    p=p';
end
% p_indiv=structfun(@(x)p(x),H.IndividualParams,'UniformOutput',false);
try
    phi=[repmat(p(H.PopulationParams),n_sim,1),p(H.IndividualParams)'];
catch 
    phi=[repmat(p(H.FixedParams),n_sim,1),p(H.IndividualParams)'];
end