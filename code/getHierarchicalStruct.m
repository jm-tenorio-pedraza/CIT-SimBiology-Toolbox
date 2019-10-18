function H = getHierarchicalStruct(params, PI, varargin)
p=inputParser;
p.addParameter('n_sigma',1,@isnumeric);
p.addParameter('rand_indx',0,@isnumeric);
p.addParameter('n_indiv',1,@isnumeric);
p.addParameter('cell_indx',[],@isnumeric);
p.parse(varargin{:})
p=p.Results;

rand_indx = p.rand_indx;                                                       % Indexes of the parameters to vary at the individual level
cell_indx = p.cell_indx;                                                    

n_rand = length(rand_indx);
n_cell = length(cell_indx);
n_indiv = size(PI.data,1);
n_sigma=p.n_sigma;
try
cell_groups = cellfun(@(x)x(1:find(x=='_',1)-1),[PI.data(:).Group],...
    'UniformOutput',false);
catch
    cell_groups = cellfun(@(x)x(1:find(x=='_',1)-1),{PI.data(:).Group},...
    'UniformOutput',false);
end
celltypes = unique(cell_groups);                                                % Extract the cell types in the dataset
n_celltypes = length(celltypes);
%% Basic Hierarchical structure
H.PopulationParams = 1:length(params);                                      % Create population parameters index
[H.IndividualParams(1:n_rand).name] = params{rand_indx};                    % Add names of the individual parameters
[H.CellParams(1:n_cell).name] = params{cell_indx};
CellIndex=mat2cell(reshape(length(params)+1:length(params)...               % Create indexes for each parameter varying at the cell level
    +n_cell*(n_celltypes),...                                                             % starting at the end of the population params
    [],n_rand)',ones(n_rand,1));

IndivIndex=mat2cell(reshape(length(params)+1+n_cell*n_celltypes:length(params)+...
    n_cell*(n_celltypes)+n_indiv*n_rand,...                                               % Create indexes of individually varying parameters 
    [],n_rand)',ones(n_rand,1));                                            % for each simulation starting at the end of the cell params


IndivEtaIndex = mat2cell(reshape(rand_indx,n_rand,[]),ones(n_rand,1));
CellEtaIndex = mat2cell(reshape(cell_indx,n_cell,[]),ones(n_cell,1));

[H.IndividualParams(1:n_rand).Index] = IndivIndex{:,:};
[H.CellParams(1:n_cell).Index] = CellIndex{:,:};
%% Add cellular parameter indexes if available
[H.CellParams(1:n_cell).EtaIndex] = CellEtaIndex{:,:};

try
    CellOmegaIndex = mat2cell(reshape(IndivIndex{end,end}(end)+1:...       % Create indexes for the variance parameters of
        IndivIndex{end,end}(end)+n_cell,n_cell,[]),ones(n_cell,1));         % individually-varying params starting at the end of
                                                                            % the individual params                                                                            
    H.SigmaParams = [CellOmegaIndex{:,:}];                                 % Add omega indexes to vector of variance indexes
    CellIndx = cellfun(@(x)ismember(cell_groups,x),celltypes,'UniformOutput', false);
    CellIndx = reshape([CellIndx{:,:}],n_indiv,n_celltypes);
    H.CellIndx = CellIndx;
catch
    CellOmegaIndex ={};                                                    % If there are no individual parameters    
end

[H.CellParams(cell_indx).OmegaIndex] = CellOmegaIndex{:,:};                % Add omega indexes of cell to H

%% Add individual parameter indexes if available
[H.IndividualParams(1:n_rand).EtaIndex] = IndivEtaIndex{:,:};
try
    IndivOmegaIndex = mat2cell(reshape(IndivIndex{end,end}(end)+1+n_cell:...       % Create indexes for the variance parameters of
        IndivIndex{end,end}(end)+n_rand+n_cell,n_rand,[]),ones(n_rand,1));         % individually-varying params starting at the end of
                                                                            % the individual params
                                                                            
    H.SigmaParams = [H.SigmaParams [IndivOmegaIndex{:,:}]];                                       % Add noise variance indexes to the vector of variance indexes
    H.SigmaParams = [H.SigmaParams H.SigmaParams(end)+1:H.SigmaParams(end)+n_sigma];
catch
    IndivOmegaIndex ={};                                                    % If there are no individual parameters
    H.SigmaParams = length(params)+1:length(params)+n_sigma;
    
end


[H.IndividualParams(rand_indx).OmegaIndex] = IndivOmegaIndex{:,:};          % Add omega indexes to H



end