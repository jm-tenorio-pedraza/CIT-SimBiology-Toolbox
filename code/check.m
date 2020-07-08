function [X, p_X] = check(X,p_X_mean,p_X)

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
return