%% Generate cross-validation datasets from original data
MetaPI = getPICrossValData('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat',...
    stateVar,groups_subset,observables,'zeroAction', 'omit','mergePhenotypes', true);

%% Add fields to each PI structure


for l=1:length(MetaPI)
    tv = arrayfun(@(x) x.dataValue(~isnan(x.dataValue(:,1)),1),MetaPI(l).PI.data,'UniformOutput',false);
    [MetaPI(l).PI.data(1:end).TV] = tv{:,:};
    MetaPI(l).PI.variableUnits=PI.variableUnits;
    MetaPI(l).PI.normIndx = 7:8;
    MetaPI(l).PI.x_0 = PI.x_0;
    MetaPI(l).PI.variants =PI.variants;
    MetaPI(l).PI.u = PI.u;
    MetaPI(l).PI.H =PI.H;
end
%% Loop for optimization

