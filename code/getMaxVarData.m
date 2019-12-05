function [train_indx] = getMaxVarData(data,time,n,varargin)
p=inputParser;
p.addParameter('criteria', 'maxVar')
p.parse(varargin{:})
p=p.Results;
n_data = length(data);

K = cellfun(@(x)max(x),time,'UniformOutput', true);
min_K = min(K);
dataValues_K = nan(n_data,1);
for i=1:n_data
   dataValues_K(i) = data{i}(time{i}==min_K,1);
end
mean_data = mean(dataValues_K);
deviations = dataValues_K - mean_data;
[~, indx] = sort(deviations);
train_indx = repelem({'nan'}, n_data,1);
    lower = floor(n/2);
    upper = ceil(n/2);
    train_indx(indx(1:lower)) = repelem({'train'}, lower,1);
    train_indx(indx(end-upper+1:end)) = repelem({'train'}, upper,1);
    train_indx(indx(lower+1:end-upper)) = repelem({'test'}, n_data-n,1);
return
