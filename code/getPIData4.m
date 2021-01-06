function PI = getPIData4(data_ext, stateVar, groups_subset, varargin)

% Inputs:
%   data_ext:   full path extension of data files
%   stateVar:   state variables to sample from the complete dataset
%   groups_subset:  Subset of treatments to choose from
%   Optional:
%       mergePhenotypes: false (default) keeps mean groups with different ...
%                       names but same group variable as separate entries
%       output: 'mean' (default) calculate mean over the
%       individually-observed data
%       grouping: 'original' (default) groups observations by the original
%       'group' variable in the PI. 'response' groups by responders vs
%       non-responders.
%       kineticGrouping: true groups by slow and fast growing tumors.
p=inputParser;
p.addParameter('mergePhenotypes', false)
p.addParameter('output', 'mean')
p.addParameter('responseGrouping', false)
p.addParameter('kineticGrouping', false)
p.addParameter('zeroHandling', 'rm')
p.parse(varargin{:})
p=p.Results;

PI_struct = cellfun(@(x)load(x,'PI'), data_ext);                           % load data structure previously processed
groups = {};
data_subset = {};
dataTime={};
SD_subset = {};
data_indx = 1;

for i=1:length(PI_struct)
    if size(PI_struct(i).PI.data,1)>size(PI_struct(i).PI.data,2)           % transpose data array for following processing
    else
        PI_struct(i).PI.data = PI_struct(i).PI.data';
    end
    [~,varindx]= (ismember(stateVar,PI_struct(i).PI.stateVar));            % identify and select variables of interest
    if any(varindx==0)
        data_subset_i = arrayfun(@(x)nan(size(x.dataValue,1), length(stateVar)),...
            PI_struct(i).PI.data, 'UniformOutput', false);
        SD_subset_i = arrayfun(@(x)nan(size(x.SD,1), length(stateVar)),...
            PI_struct(i).PI.data, 'UniformOutput', false);
        varindx = varindx~=0;
        for j=1:length(data_subset_i)
            data_subset_i{j,1}(:,varindx) = PI_struct(i).PI.data(j).dataValue(:, varindx);
            SD_subset_i{j,1}(:,varindx) = PI_struct(i).PI.data(j).SD(:, varindx);
        end
    else
        data_subset_i = arrayfun(@(x) x.dataValue(:,varindx),...           % Remove variables not to be considered
            PI_struct(i).PI.data,'UniformOutput', false)';
        SD_subset_i = arrayfun(@(x) x.SD(:,varindx),...                    % Remove variables not to be considered
            PI_struct(i).PI.data,'UniformOutput', false)';
    end
    
    groups_i = {PI_struct(i).PI.data(:).Group}';                             % extract groups cell array
    groupIndx = ismember(groups_i, groups_subset);                           % Identify and select the therapy groups
    groups_i = groups_i(groupIndx);
    
    data_subset_i = data_subset_i(groupIndx);                               
    SD_subset_i = SD_subset_i(groupIndx);                  % Select the measured SDs if they are available
    dataTime_i = {PI_struct(i).PI.data(groupIndx).dataTime}';                % extract time-points cell array
   
    data_subset(data_indx:(data_indx+length(data_subset_i)-1)) = data_subset_i;
    SD_subset(data_indx:(data_indx+length(SD_subset_i)-1)) = SD_subset_i;
    dataTime(data_indx:(data_indx+length(dataTime_i)-1)) = dataTime_i;
    groups(data_indx:(data_indx+length(groups_i)-1)) = groups_i;
    data_indx = data_indx+length(data_subset_i);
end

dataTimeZeroIndx = cellfun(@(x) x(1)==0, dataTime);
for i=1:length(dataTime)
    if dataTimeZeroIndx(i)
        
        data_subset{i} = data_subset{i}(2:end,:);
        dataTime{i} = dataTime{i}(2:end);
        SD_subset{i} = SD_subset{i}(2:end,:);
    end
end
    data = [];
    struct_indx = 1;
    unique_groups = unique(groups);                                        % Determine the groups found in the dataset

if  ~p.responseGrouping && ~(p.kineticGrouping)
    case_i = 1;
elseif p.responseGrouping && ~(p.kineticGrouping)
    case_i = 2;
else
    if ~p.responseGrouping && (p.kineticGrouping)
        case_i = 3;
    else
        case_i = 4;
    end
end

for i=1:length(unique_groups)
    group_i = ismember(groups', unique_groups{i});
    data_i = data_subset(group_i)';
    dataTime_i = dataTime(group_i)';
    SD_indx = SD_subset(group_i);                                    % Array of standard deviations cells if dealing with mean data, empty otherwise

    try
        SD_i = cellfun(@(x) x(:,varindx),SD_indx, 'UniformOutput', false);
    catch
        SD_i = cellfun(@(x) nan(length(x),1), dataTime_i, 'UniformOutput', false);
    end
        
    time = unique(cell2mat(dataTime_i));                                        % Count how many time-resolved data points there are for each individual
    n_i = length(data_i);
    mean_i_p = NaN(length(time), length(stateVar));                         % Allocate space for mean of all variables (Non-responders, fast growing tumors matrix)
    sd_i_p = NaN(length(time), length(stateVar));
    tvIndx = ismember(stateVar, {'Tumor' 'TV'});
    TV_i = cellfun(@(x) x(:,tvIndx),data_i, 'UniformOutput',false);              % Extract tumor volumes to determine classification
    
    
    switch case_i
        case 1
        case 3
            time_end =  cellfun(@(x) x(end), dataTime_i);
            kinIndx = time_end<median(time_end);
            
            mean_i_f = NaN(length(time),length(stateVar));
            sd_i_f = NaN(length(time),length(stateVar));                     % Allocate space for mean of all variables (fast growing tumors matrix)
            mean_i_s = NaN(length(time),length(stateVar));
            sd_i_s = NaN(length(time),length(stateVar));                     % Allocate space for mean of all variables (slow growing tumors matrix)
            
            kinetic = repelem({'NaN'}, n_i,1);
            kinetic(kinIndx) = {'Fast'};
            kinetic(~kinIndx) = {'Slow'};
            
        case 2
            respIndx = cellfun(@(x) x(end)>0.1 || all(isnan(x)) || x(end-1)>0.1, TV_i);
            
            mean_i_r = NaN(length(time),length(stateVar));
            sd_i_r = NaN(length(time),length(stateVar));                     % Allocate space for mean of all variables (Responders, slow growing tumors matrix)
            
            response = repelem({'NaN'}, n_i,1);
            response(respIndx) = {'Non-responder'};
            response(~respIndx) =  {'Responder'};
            
        case 4
            temp = [];
             nan_indx = cellfun(@(x) isnan(x), TV_i,'UniformOutput', false);
            [temp(1:length(TV_i)).time] = dataTime_i{:,:};
            [temp(1:length(TV_i)).nanIndx] = nan_indx{:,:};
            [temp(1:length(TV_i)).TV] = TV_i{:,:};

            time_an = arrayfun(@(x) x.time(~x.nanIndx), temp, 'UniformOutput', false);
            tv_an = arrayfun(@(x) x.TV(~x.nanIndx), temp, 'UniformOutput', false);
            try
            respIndx = cellfun(@(x) x(end)>0.01 || isnan(x(end)) || all(isnan(x)), tv_an);
            catch
                respIndx = true;
            end
            try
            time_end_p =  cellfun(@(x) x(end), time_an(respIndx));
            time_end_r =  cellfun(@(x) x(end), time_an(~respIndx));

            catch
                time_an = arrayfun(@(x) x.time, temp, 'UniformOutput', false);
                time_end_p =  cellfun(@(x) x(end), time_an(respIndx));
                time_end_r =  cellfun(@(x) x(end), time_an(~respIndx));
            end
            endValueIndx_r = cellfun(@(x) x(end)==0, TV_i(~respIndx));

            kinIndx_p = (time_end_p<=median(time_end_p));
            kinIndx_r = or(time_end_r<median(time_end_r), endValueIndx_r');
            
            mean_i_f_p = NaN(length(time),length(stateVar));
            sd_i_f_p = NaN(length(time),length(stateVar));
            
            mean_i_s_p = NaN(length(time),length(stateVar));
            sd_i_s_p = NaN(length(time),length(stateVar));
            
            
            mean_i_f_r = NaN(length(time),length(stateVar));
            sd_i_f_r = NaN(length(time),length(stateVar));
            
            mean_i_s_r = NaN(length(time),length(stateVar));
            sd_i_s_r = NaN(length(time),length(stateVar));
            
            kinetic = repelem({'NaN'}, n_i,1);
            kinetic(kinIndx_p) = {'Fast'};
            kinetic(kinIndx_r) = {'Fast'};
            kinetic(~kinIndx_p) = {'Slow'};
            kinetic(~kinIndx_r) = {'Slow'};
            
            response = repelem({'NaN'}, n_i,1);
            response(respIndx) = {'Non-responder'};
            response(~respIndx) =  {'Responder'};
    end
    
    if strcmp(p.output, 'mean')
        for j=1:length(stateVar)
            data_ij = cellfun(@(x) x(:,j), data_i, 'UniformOutput', false);
            try
                sd_ij = cellfun(@(x) x(:,j), SD_i, 'UniformOutput', false);
            catch
                sd_ij = cellfun(@(x) nan(length(x),1), dataTime_i, 'UniformOutput', false);
            end
            mat_ij = nan(length(time), n_i);
            matSD_ij = nan(length(time), n_i);
            for k=1:n_i
                timeIndx_ij = ismember(time, dataTime_i{k});
                mat_ij(timeIndx_ij, k) = data_ij{k};
                if n_i==1
                matSD_ij(timeIndx_ij, k) = sd_ij{k};
                end
            end
            
            switch case_i
                case 1
                    mean_i_p(:,j) = mean(mat_ij,2, 'omitnan');
                    sd_i_p(:,j) = std(mat_ij,[],2, 'omitnan');
                    if n_i==1
                        sd_i_p(:,j) = (matSD_ij{:,:});
                    else
                        sd_i_p(:,j) = std(mat_ij(:,:),[],2, 'omitnan');
                    end
                   
                case 2
                    mean_i_p(:,j) = mean(mat_ij(:,respIndx),2, 'omitnan');
                    mean_i_r(:,j) = mean(mat_ij(:,~respIndx),2, 'omitnan');
                    if n_i==1
                        sd_i_p(:,j) = (matSD_ij);
                    else
                        sd_i_p(:,j) = std(mat_ij(:,respIndx),[],2, 'omitnan');
                        try
                        sd_i_r(:,j) = std(mat_ij(:,~respIndx),[],2, 'omitnan');
                        catch 
                        end
                    end
                    
                    
                case 3
                    mean_i_f(:,j) = mean(mat_ij(:,kinIndx),2, 'omitnan');
                    mean_i_s(:,j) = mean(mat_ij(:,~kinIndx),2, 'omitnan');
                    
                    sd_i_f(:,j) =  std(mat_ij(:,kinIndx),[],2, 'omitnan');
                    sd_i_s(:,j) = std(mat_ij(:,~kinIndx),[],2, 'omitnan');
                    if n_i==1
                        sd_i_f(:,j) = (matSD_ij);
                    else
                        sd_i_f(:,j) = std(mat_ij(:,kinIndx),[],2, 'omitnan');
                        try
                        sd_i_s(:,j) = std(mat_ij(:,~kinIndx),[],2, 'omitnan');
                        catch
                        end
                    end
                    
                case 4
                    nonRespIndx = find(respIndx);
                    RespIndx = find(~respIndx);
                    mean_i_f_p(:,j) = mean(mat_ij(:,nonRespIndx(kinIndx_p)),2, 'omitnan');
                    mean_i_f_r(:,j) = mean(mat_ij(:,RespIndx(kinIndx_r)),2, 'omitnan');
                    mean_i_s_p(:,j) = mean(mat_ij(:,nonRespIndx(~kinIndx_p)),2, 'omitnan');
                    mean_i_s_r(:,j) = mean(mat_ij(:,RespIndx(~kinIndx_r)),2, 'omitnan');
                    if n_i == 1
                        sd_i_f_p(:,j) = matSD_ij;
                    else
                        sd_i_f_p(:,j) = std(mat_ij(:,kinIndx_p),[],2, 'omitnan');
                        try
                        sd_i_f_r(:,j) = std(mat_ij(:,kinIndx_r),[],2, 'omitnan');
                        sd_i_s_p(:,j) = std(mat_ij(:,~kinIndx_p),[],2, 'omitnan');
                        sd_i_s_r(:,j) = std(mat_ij(:,~kinIndx_r),[],2, 'omitnan');
                        catch
                        end
                    end
                    
                  
            end
        end
        
        switch case_i
            case 1
                 data(struct_indx).Name = unique_groups{i};
                    data(struct_indx).dataTime = time;
                    data(struct_indx).dataValue = mean_i_p;
                    data(struct_indx).SD = sd_i_p;
                    data(struct_indx+1).Group = unique_groups(i);
                    struct_indx = struct_indx +1;
                                        

            case 2
                    data(struct_indx).Name = strjoin({unique_groups{i}, 'Progressor'}, '_');
                    data(struct_indx).dataTime = time;
                    data(struct_indx).dataValue = mean_i_p;
                    data(struct_indx).SD = sd_i_p;
                    data(struct_indx).Group = unique_groups(i);
                    data(struct_indx).Response = 'Progressor';
                    
                    data(struct_indx+1).Name = strjoin({unique_groups{i}, 'Responders'}, '_');
                    data(struct_indx+1).dataTime = time;
                    data(struct_indx+1).dataValue = mean_i_r;
                    data(struct_indx+1).SD = sd_i_r;
                    data(struct_indx+1).Group = unique_groups(i);
                    data(struct_indx+1).Response = 'Responder';
                    struct_indx = struct_indx +2;

            case 3
                    data(struct_indx).Name = strjoin({unique_groups{i}, 'FastGrowth'}, '_');
                    data(struct_indx).dataTime = time;
                    data(struct_indx).dataValue = mean_i_f;
                    data(struct_indx).SD = sd_i_f;
                    data(struct_indx).Group = unique_groups(i);
                    data(struct_indx).Kinetic = 'FastGrowth';
                    
                    data(struct_indx+1).Name = strjoin({unique_groups{i}, 'SlowGrowth'}, '_');
                    data(struct_indx+1).dataTime = time;
                    data(struct_indx+1).dataValue = mean_i_s;
                    data(struct_indx+1).SD = sd_i_s;
                    data(struct_indx+1).Group = unique_groups(i);
                    data(struct_indx+1).Kinetic = 'SlowGrowth';
                    struct_indx = struct_indx +2;

                    
            case 4
                  data(struct_indx).Name = strjoin({unique_groups{i}, 'FastGrowth' 'Progressor'}, '_');
                    data(struct_indx).dataTime = time;
                    data(struct_indx).dataValue = mean_i_f_p;
                    data(struct_indx).SD = sd_i_f_p;
                    data(struct_indx).Group = unique_groups(i);
                    data(struct_indx).Response = 'Progressor';
                    data(struct_indx).Kinetic = 'FastGrowth';
                    
                    data(struct_indx+1).Name = strjoin({unique_groups{i}, 'FastGrowth' 'Responder'}, '_');
                    data(struct_indx+1).dataTime = time;
                    data(struct_indx+1).dataValue = mean_i_f_r;
                    data(struct_indx+1).SD = sd_i_f_r;
                    data(struct_indx+1).Group = unique_groups(i);
                    data(struct_indx+1).Response = 'Responder';
                    data(struct_indx+1).Kinetic = 'FastGrowth';
                    
                    data(struct_indx+2).Name = strjoin({unique_groups{i}, 'SlowGrowth' 'Progressor'}, '_');
                    data(struct_indx+2).dataTime = time;
                    data(struct_indx+2).dataValue = mean_i_s_p;
                    data(struct_indx+2).SD = sd_i_s_p;
                    data(struct_indx+2).Group = unique_groups(i);
                    data(struct_indx+2).Response = 'Progressor';
                    data(struct_indx+2).Kinetic = 'SlowGrowth';
                    
                    data(struct_indx+3).Name = strjoin({unique_groups{i}, 'SlowGrowth' 'Responder'}, '_');
                    data(struct_indx+3).dataTime = time;
                    data(struct_indx+3).dataValue = mean_i_s_r;
                    data(struct_indx+3).SD = sd_i_s_r;
                    data(struct_indx+3).Group = unique_groups(i);
                    data(struct_indx+3).Response = 'Responder';
                    data(struct_indx+3).Kinetic = 'SlowGrowth';
                    struct_indx = struct_indx +4;

        end
        
        
        
    else
        [data(struct_indx:struct_indx+n_i-1).dataTime] = dataTime_i{:,:};
        [data(struct_indx:struct_indx+n_i-1).dataValue] = data_i{:,:};
        [data(struct_indx:struct_indx+n_i-1).SD] = SD_i{:,:};
        [data(struct_indx:struct_indx+n_i-1).Group] = groups{group_i,:};
        
        if p.responseGrouping
            [data(struct_indx:struct_indx+n_i-1).Response] = response{:,:};
        end
        
        if p.kineticGrouping
            [data(struct_indx:struct_indx+n_i-1).Growth] = kinetic{:,:};
        end
         struct_indx = struct_indx +length(data_i);

    end
end

PI.data = data;
nanIndx = arrayfun(@(x) all(all(isnan(x.dataValue))), data);                % Identify and eliminate individuals where all obs are nan
PI.data = data(:,~nanIndx)';

try
PI.data(ismember([PI.data(:).Group], 'MOC1_Control_Mean')).Group ...
    = {'MOC1_Control'};
PI.data(ismember([PI.data(:).Group], 'MOC2_Control_Mean')).Group ...
    = {'MOC2_Control'};
catch
end

%% Control for data with non-admissible values (x=0)

zero_indx = arrayfun(@(x) logical(sum(x.dataValue(:,1)==0,2)),PI.data,...
    'UniformOutput', false);                                                % Identify indexes of tumor volume equal to 0
[PI.data(1:end).zero_indx] = zero_indx{:,:};                                % Add index array to data array
if strcmp(p.zeroHandling, 'rm')
    dataValue = arrayfun(@(x) x.dataValue(~x.zero_indx,:), PI.data,...          % Subset out the zero values and their respective time points
        'UniformOutput', false);
    dataTime = arrayfun(@(x) x.dataTime(~x.zero_indx,:), PI.data,...
        'UniformOutput', false);
        SD = arrayfun(@(x) x.SD(~x.zero_indx,:), PI.data,...          % Subset out the zero values and their respective time points
        'UniformOutput', false);
    [PI.data(1:end).dataValue] = dataValue{:,:};
[PI.data(1:end).dataTime] = dataTime{:,:};
[PI.data(1:end).SD] = SD{:,:};

elseif strcmp(p.zeroHandling, 'imput')
    [minValues, minValuesIndx] = arrayfun(@(x) min(x.dataValue(~x.zero_indx,1), [], 1,'omitnan'), PI.data,...
    'UniformOutput', false);
    [PI.data(1:end).minValues] = minValues{:,:};
    for i = 1:length(PI.data)
%         if PI.data(i).zero_indx(1)==0
            try
                PI.data(i).dataValue((PI.data(i).zero_indx),1) = minValues{i};
                PI.data(i).SD(PI.data(i).zero_indx,1) = PI.data(i).SD(minValuesIndx{i}, 1);
            catch
            end
%         else
%             PI.data(i).dataValue = PI.data(i).dataValue(2:end,:);
%             PI.data(i).dataTime = PI.data(i).dataTime(2:end,:);
%             PI.data(i).SD = PI.data(i).SD(2:end,:);
%         end
    end
end



PI.tspan = unique(cat(1,PI.data(:).dataTime));
PI.n_data=sum(cellfun(@(x)sum(sum(~isnan(x))),{PI.data.dataValue},...
    'UniformOutput',true));

if size(PI.data,1)>size(PI.data,2)
else
    PI.data=PI.data';
end

return
