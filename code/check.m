function [X, p_X] = check(X,p_X_mean,p_X,varargin)
inputs = inputParser;
inputs.addParameter('N', size(X,2));
inputs.addParameter('d', size(X,1));
inputs.addParameter('dxN', true);
inputs.parse(varargin{:});
inputs=inputs.Results;

if inputs.dxN
else
    X = X';
end
p_X_iqr = quantile(p_X_mean, [0.25 0.75]); % 1 x 2 vector 
p_X_indx = p_X_mean<(p_X_iqr(1) - 2*(p_X_iqr(2)-p_X_iqr(1)));

if any(p_X_indx)
    p_X_max_indx = ismember(p_X,max(p_X));
    X_max = X(:,p_X_max_indx);
    p_X_max = p_X(p_X_max_indx);
    X(:,p_X_indx) = repmat(X_max(:,1),1,sum(p_X_indx));
    p_X(p_X_indx,1) = repelem(p_X_max(1), sum(p_X_indx),1);
else
end
if inputs.dxN
        
else
    X = X';
end

return