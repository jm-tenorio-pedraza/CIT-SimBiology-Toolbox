%% Search paths
warning on
clear all
if ispc
    addpath(genpath('\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox'))
    cd('\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\output\CIM\PI\CIM23')
    out = sbioloadproject('\Users\jmten\OneDrive\Dokumente\GitHub\sbio-projects\CIM_5.sbproj');
    
else
    addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
    cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/CIM/PI/CIM23')
    out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/sbio-projects/CIM_5.sbproj');
    
end


%% Obtain data, simulation function and dose table

    
    if ispc
    Clavijo=load('\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\PI_Clavijo.mat');
    Morisada=load('\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\PI_Morisada.mat');
    
    else
        Clavijo=load('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat');
        Morisada=load('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Morisada.mat');
    end
    
    
    %% Clavijo analysis
    vars = {'Tumor' 'CD8' 'Treg' 'DC' 'GMDSC' 'CD107a_Rel' 'Tumor_PDL1_Rel' 'Myeloid_PDL1_Rel'};
    varsIndx = ismember(Clavijo.PI.stateVar, vars);
    n_Var = sum(varsIndx);

    groups = {Clavijo.PI.data(:).Group};
    unique_groups = unique({Clavijo.PI.data(:).Group});
    N_groups = length(unique_groups);
    Clavijo.summary = [];
    for i=1:N_groups
        groupIndx = ismember(groups, unique_groups(i));
        N_group_i = sum(groupIndx);
        
        time = unique(cell2mat({Clavijo.PI.data(groupIndx).dataTime}'));
        n_Time = length(time);
        dataMean = nan(n_Time, n_Var);
        dataVar = nan(n_Time, n_Var);
        dataTimePoints = nan(n_Time, n_Var);
        dataValue = cellfun(@(x) x(:,varsIndx),{Clavijo.PI.data(groupIndx).dataValue}, 'UniformOutput', false);
        dataTime = {Clavijo.PI.data(groupIndx).dataTime}';
        for k = 1:n_Var
            dataset_k = nan(n_Time, N_group_i);
            for j=1:N_group_i
                timeIndx_j = ismember(time, dataTime{j,:});
                dataset_k(timeIndx_j, j) = dataValue{j}(:,k);
            end
            
            dataMean(:,k) = mean(dataset_k,2,'omitnan');
            dataVar(:,k) = var(dataset_k,[],2, 'omitnan');
            dataTimePoints(:,k) = sum(~isnan(dataset_k),2);
            
            if ismember(unique_groups(i), {'MOC1_Control_Mean' 'MOC2_Control_Mean'})
                dataSD = cellfun(@(x) x(:,varsIndx),{Clavijo.PI.data(groupIndx).SD}, 'UniformOutput', false);
                
                dataVar(:,k) = dataSD{j}(:,k).^2;
            end
        end
        Clavijo.summary(i).Group = unique_groups(i);
        Clavijo.summary(i).Mean = dataMean;
        Clavijo.summary(i).Var = dataVar;
        Clavijo.summary(i).N = N_group_i;
        Clavijo.summary(i).dataTimePoints = dataTimePoints;
        
        if ismember(unique_groups(i), {'MOC1_Control_Mean' 'MOC2_Control_Mean'})
            Clavijo.summary(i).N = 5;
            Clavijo.summary(i).dataTimePoints =dataTimePoints*5;
        end
    end
    
    %% Morisada
    vars = {'Tumor' 'CD8' 'Treg' 'DC' 'GMDSC' 'CD107a_Rel' 'Tumor_PDL1_Rel' 'Myeloid_PDL1_Rel'};
    varsIndx = ismember(Morisada.PI.stateVar, vars);
    n_Var = sum(varsIndx);

    groups = {Morisada.PI.data(:).Group};
    unique_groups = unique({Morisada.PI.data(:).Group});
    N_groups = length(unique_groups);
    Morisada.summary = [];
    for i=1:N_groups
        groupIndx = ismember(groups, unique_groups(i));
        N_group_i = sum(groupIndx);
        
        time = unique(cell2mat({Morisada.PI.data(groupIndx).dataTime}'));
        n_Time = length(time);
        dataMean = nan(n_Time, n_Var);
        dataVar = nan(n_Time, n_Var);
        dataTimePoints = nan(n_Time, n_Var);
        dataValue = cellfun(@(x) x(:,varsIndx),{Morisada.PI.data(groupIndx).dataValue}, 'UniformOutput', false);
        dataTime = {Morisada.PI.data(groupIndx).dataTime}';
        for k = 1:n_Var
            dataset_k = nan(n_Time, N_group_i);
            for j=1:N_group_i
                timeIndx_j = ismember(time, dataTime{j,:});
                dataset_k(timeIndx_j, j) = dataValue{j}(:,k);
            end
            
            dataMean(:,k) = mean(dataset_k,2,'omitnan');
            dataVar(:,k) = var(dataset_k,[],2, 'omitnan');
            dataTimePoints(:,k) = sum(~isnan(dataset_k),2);
            
            if ismember(unique_groups(i), {'MOC1_Control_Mean' 'MOC2_Control_Mean'})
                dataSD = cellfun(@(x) x(:,varsIndx),{Morisada.PI.data(groupIndx).SD}, 'UniformOutput', false);
                
                dataVar(:,k) = dataSD{j}(:,k).^2;
            end
        end
        Morisada.summary(i).Group = unique_groups(i);
        Morisada.summary(i).Mean = dataMean;
        Morisada.summary(i).Var = dataVar;
        Morisada.summary(i).N = N_group_i;
        Morisada.summary(i).dataTimePoints = dataTimePoints;
        
        if ismember(unique_groups(i), {'MOC1_Control_Mean' 'MOC2_Control_Mean'})
            Morisada.summary(i).N = 5;
            Morisada.summary(i).dataTimePoints =dataTimePoints*5;
        end
    end
   %% Plot
   summary = [Clavijo.summary Morisada.summary];
   N_groups = length(summary);
   ncol = ceil(sqrt(n_Var));
   nrow = ceil(n_Var/ncol);
   colors = linspecer(N_groups);
   for i=1:n_Var
       subplot(nrow, ncol,i)
       hold on
       h=arrayfun(@(x) plot(x.Mean(:,i), x.Var(:,i),'-+'), summary);
       for j=1:N_groups
           h(j).Color = colors(j,:);
           h(j).MarkerEdgeColor = colors(j,:);
           h(j).MarkerFaceColor = colors(j,:);
       end
       X = cell2mat(cellfun(@(x) x(:,i)', {summary(:).Mean}, 'UniformOutput',false));
       Y = cell2mat(cellfun(@(x) x(:,i)', {summary(:).Var}, 'UniformOutput',false));
       zeroIndx = or(X ==0, Y == 0);
       X = X(~zeroIndx);
       Y = Y(~zeroIndx);
        
       mdl = fitlm(log(X),log(Y));
       alpha = mdl.Coefficients.Estimate(1);
       beta = mdl.Coefficients.Estimate(2);
       legend(strjoin({'\alpha = ' num2str(alpha), ',' '\beta =', num2str(beta)},''))
       title(vars(i))
       xlabel('Mean')
       ylabel('Var')
       set(gca,'XScale', 'log', 'YScale', 'log');
       ax = gca;
       minL =  min([ax.XLim(1) ax.YLim(1)]);
       maxL = max([ax.XLim(2) ax.YLim(2)]);
       plot(minL:.01:maxL, minL:.01:maxL, '-k')
       plot(minL:.01:maxL, exp(alpha+beta*log(minL:.01:maxL)),'-r')
       plot(minL:.01:maxL, (minL:.01:maxL).^2/.9,'-g')
   end
   %legend(unique_groups)
       