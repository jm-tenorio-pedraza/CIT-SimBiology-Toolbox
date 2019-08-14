function par_struct = getParamStruct(simFun,Pop_indx, Ind_indx,n_sim, p,sigmaNames)


paramNames=[simFun.Parameters.Name(Pop_indx); repelem(simFun.Parameters.Name(Ind_indx),n_sim,1);...
    sigmaNames];

par_struct=struct('name', paramNames,'minValue',num2cell(p(:,1)), 'maxValue',...
    num2cell(p(:,3)),'startValue', num2cell(p(:,2)),'mu_prior',  num2cell(p(:,4)),'sigma_prior',  num2cell(p(:,5)));
return