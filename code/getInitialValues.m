function x_0 = getInitialValues(groups, initialstruct)

x_0 = nan(length(groups),1);
groups = cellfun(@(x) x(1:regexp(x,'_')-1),groups, 'UniformOutput',false);

for i = 1:size(initialstruct,1)
    cell_indx = ismember(groups, initialstruct(i).name);
    x_0(cell_indx,1) = initialstruct(i).initialValue;
end
return