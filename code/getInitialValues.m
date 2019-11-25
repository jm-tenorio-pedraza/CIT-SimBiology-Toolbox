function [x_0, variant] = getInitialValues(groups, initialstruct,varargin)
par= inputParser;
par.addParameter('parameters', {})
par.parse(varargin{:})
par=par.Results;
x_0 = nan(length(groups),1+length(par.parameters));
variant = repelem({nan},length(groups),1);

groups = cellfun(@(x) x(1:regexp(x,'_')-1),groups, 'UniformOutput',false);

for i = 1:size(initialstruct,1)
    
    cell_indx = ismember(groups, initialstruct(i).name);
    x_0(cell_indx,1) = initialstruct(i).initialValue;
    
    if ~isempty(par.parameters)
        for j=1:length(par.parameters)
            par_i = par.parameters(j);
            variant_i = get(initialstruct(i).variant, 'Content');
            values_i = cellfun(@(x)strcmp(x{2},par_i),variant_i);
            values_i = variant_i{values_i}{4};
            x_0(cell_indx,1+j)=values_i;
        end
    end
end
return