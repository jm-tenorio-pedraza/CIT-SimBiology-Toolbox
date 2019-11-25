function logP=likelihood_censored(log_data_value,log_mean,sigma,direction,varargin)
if nargin<3
    error('likelihood_censored:toofewinputs','likelihood_censored requires atleast 4 inputs.')
end

if strcmp(direction,'right')
    logP=log(1-1/2*(1+erf((log_data_value-log_mean)/(sqrt(2)*sigma))));
elseif strcmp(direction, 'left')
    logP=log(1/2*(1+erf((log_data_value-log_mean)/(sqrt(2)*sigma))));
else
    logP = 0;
end
