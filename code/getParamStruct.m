function par_struct = getParamStruct(simFun,Pop_indx, Ind_indx,n_sim, p)


paramNames=[simFun.Parameters.Name(Pop_indx); repelem(simFun.Parameters.Name(Ind_indx),n_sim);...
    {'b_TV'; 'b_CD8'; 'b_CD107a'; 'b_DC'; 'b_MDSC'; 'b_PDL1'}];

par_struct=struct('Name', paramNames,'minValue',num2cell(p(:,1)), 'maxValue',...
    num2cell(p(:,3)),'startValue', num2cell(p(:,2)));
return