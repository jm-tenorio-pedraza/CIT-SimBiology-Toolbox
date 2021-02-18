function PI = getPIData3(data_ext, stateVar, groups_subset, varargin)

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
p.parse(varargin{:})
p=p.Results;

load(data_ext,'PI')                                                         % load data structure previously processed
if size(PI.data,1)>size(PI.data,2)                                          % transpose data array for following processing
else
    PI.data = PI.data';
end

[~,varindx]= (ismember(stateVar,PI.stateVar));                              % identify and select variables of interest
varindx = varindx(~varindx==0);

data_subset = arrayfun(@(x) x.dataValue(:,varindx), PI.data,...
    'UniformOutput', false)';
try
groups = {PI.data(:).Group}';                                               % extract groups cell array
groupIndx = ismember(groups, groups_subset);
catch
    groups = [PI.data(:).Group]';                                               % extract groups cell array
groupIndx = ismember(groups, groups_subset);
end
groups = groups(groupIndx);
data_subset = data_subset(groupIndx);
try
    SD_subset = {PI.data(groupIndx).SD};
catch
    SD_subset = cellfun(@(x) nan(length(x),length(PI.stateVar)), data_subset,'UniformOutput', false);
end
dataTime = {PI.data(groupIndx).dataTime}';                                          % extract time-points cell array
unique_groups = unique(groups);
data = [];
struct_indx = 1;

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
    group_i = ismember(groups, unique_groups{i});
    data_i = data_subset(group_i);
    dataTime_i = dataTime(group_i);
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
            kinIndx_r = or(time_end_r<=median(time_end_r), endValueIndx_r);
            
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
                        sd_i_p(:,j) = (matSD_ij(:,:));
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

zero_indx = arrayfun(@(x) sum(x.dataValue(:,1)==0,2),PI.data,...
    'UniformOutput', false);                                                % Identify indexes of tumor volume equal to 0
[PI.data(1:end).zero_indx] = zero_indx{:,:};                                % Add index array to data array
dataValue = arrayfun(@(x) x.dataValue(~x.zero_indx,:), PI.data,...          % Subset out the zero values and their respective time points
    'UniformOutput', false);
dataTime = arrayfun(@(x) x.dataTime(~x.zero_indx,:), PI.data,...
    'UniformOutput', false);

[PI.data(1:end).dataValue] = dataValue{:,:};
[PI.data(1:end).dataTime] = dataTime{:,:};
try
    SD = arrayfun(@(x) x.SD(~x.zero_indx,:), PI.data,...          % Subset out the zero values and their respective time points
        'UniformOutput', false);
    [PI.data(1:end).SD] = SD{:,:};
catch
end

PI.tspan = unique(cat(1,PI.data(:).dataTime));
PI.n_data=sum(cellfun(@(x)sum(sum(~isnan(x))),{PI.data.dataValue},...
    'UniformOutput',true));

if size(PI.data,1)>size(PI.data,2)
else
    PI.data=PI.data';
end

return
