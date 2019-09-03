function sigma=getCovariance(p,H)
% Input:    1) parameter vector of dimension (px1)
%           2) H, structure array with fields:
%               PopulationParams: Indexes of population-wise parameters in
%               p
%               IndividualParams: Indexes of parameters that vary at the...
%                                   individual level
%               SigmaParams: Indexes of parameters that quantify variable
%               measurement error standard deviation
% Output:   Matrix (1xn_var) of sigmas
sigma=p(setdiff(H.SigmaParams, [H.IndividualParams.OmegaIndex]));

end