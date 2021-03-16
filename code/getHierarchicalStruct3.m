function H = getHierarchicalStruct3(params, PI, varargin)
p=inputParser;
p.addParameter('n_sigma',1,@isnumeric);
p.addParameter('rand_indx',[],@isnumeric);
p.addParameter('n_indiv',1,@isnumeric);
p.addParameter('cell_indx',[],@isnumeric);
p.addParameter('resp_indx',[],@isnumeric);
p.addParameter('CellField', 'Group');
p.addParameter('ResponseField', 'Response');
p.parse(varargin{:})
p=p.Results;

rand_indx = p.rand_indx;                                                       % Indexes of the parameters to vary at the individual level
cell_indx = p.cell_indx;                                                    
resp_indx = p.resp_indx;

n_rand = length(rand_indx);
n_cell = length(cell_indx);
n_resp = length(resp_indx);
n_indiv = length(PI.data);
n_sigma=p.n_sigma;
try
    cell_groups = cellfun(@(x)x(1:find(x=='_',1)-1),[PI.data(:).(p.CellField)],...
        'UniformOutput',false);
catch
    cell_groups = cellfun(@(x)x(1:find(x=='_',1)-1),{PI.data(:).(p.CellField)},...
    'UniformOutput',false);
end
response_groups = ({PI.data(:).(p.ResponseField)});
if ischar(response_groups) 
        response_groups = ([PI.data(:).(p.ResponseField)]);     
end    
celltypes = unique(cell_groups);                                                % Extract the cell types in the dataset
n_celltypes = length(celltypes);

resptypes = unique(response_groups);
n_resptypes = length(resptypes);

%% Basic Hierarchical structure
H.PopulationParams = 1:length(params);                                      % Create population parameters index

cellIndexStart = length(params) +1;
cellIndexEnd = cellIndexStart + n_cell*n_celltypes-1;
CellIndex=mat2cell(reshape(cellIndexStart:cellIndexEnd,...               % Create indexes for each parameter varying at the cell level
    [],n_cell)',ones(n_cell,1));
CellEtaIndex = mat2cell(reshape(cell_indx,n_cell,[]),ones(n_cell,1));

indivIndexStart = length(params)+1+n_cell*n_celltypes;
indivIndexEnd = indivIndexStart+n_indiv*n_rand-1;
IndivIndex=mat2cell(reshape(indivIndexStart:indivIndexEnd,...                                               % Create indexes of individually varying parameters 
    [],n_rand)',ones(n_rand,1));                                            % for each simulation starting at the end of the cell params
IndivEtaIndex = mat2cell(reshape(rand_indx,n_rand,[]),ones(n_rand,1));

respIndexStart = length(params) + 1 + n_cell*n_celltypes + n_indiv*n_rand;
respIndexEnd = respIndexStart+n_resp*n_resptypes-1;
RespIndex = mat2cell(reshape(respIndexStart:respIndexEnd, [], n_resp)',...
    ones(n_resp,1));
RespEtaIndex = mat2cell(reshape(resp_indx, n_resp, []), ones(n_resp,1));

[H.CellParams(1:n_cell).name] = params{cell_indx};
[H.CellParams(1:n_cell).Index] = CellIndex{:,:};

[H.IndividualParams(1:n_rand).name] = params{rand_indx};                    % Add names of the individual parameters
[H.IndividualParams(1:n_rand).Index] = IndivIndex{:,:};

[H.RespParams(1:n_resp).name] = params{resp_indx};
[H.RespParams(1:n_resp).Index] = RespIndex{:,:};

%% Add cellular parameter indexes if available
[H.CellParams(1:n_cell).EtaIndex] = CellEtaIndex{:,:};

try
    if isempty(IndivIndex) && isempty(RespIndex)
        CellOmegaIndex = mat2cell(reshape(CellIndex{end,end}(end)+1:...       % Create indexes for the variance parameters of
            CellIndex{end,end}(end)+n_cell,n_cell,[]),ones(n_cell,1));
    elseif isempty(RespIndex)
        CellOmegaIndex = mat2cell(reshape(IndivIndex{end,end}(end)+1:...       % Create indexes for the variance parameters of
            IndivIndex{end,end}(end)+n_cell,n_cell,[]),ones(n_cell,1));         % individually-varying params starting at the end of
                                                                                % the individual params    
    else
        CellOmegaIndex =  mat2cell(reshape(IndivIndex{end,end}(end)+1:...       % Create indexes for the variance parameters of
            IndivIndex{end,end}(end)+n_cell,n_cell,[]),ones(n_cell,1));         % individually-varying params starting at the end of
    end                                                                         % the response params
    
catch
    CellOmegaIndex ={};                                                    % If there are no cell parameters    
end
CellIndx = cellfun(@(x)ismember(cell_groups,x),celltypes,'UniformOutput', false);
CellIndx = reshape([CellIndx{:,:}],n_indiv,n_celltypes);

RespIndx = cellfun(@(x)ismember(response_groups,x), resptypes, 'UniformOutput', false);
RespIndx = reshape([RespIndx{:,:}], n_indiv, n_resptypes);

H.CellIndx = CellIndx;
H.CellTypes = celltypes;
H.RespIndx = RespIndx;
H.RespTypes = resptypes;


H.SigmaParams = [CellOmegaIndex{:,:}];                                 % Add omega indexes to vector of variance indexes
[H.CellParams(1:n_cell).OmegaIndex] = CellOmegaIndex{:,:};          % Add omega indexes to H

%% Add individual parameter indexes if available
[H.IndividualParams(1:n_rand).EtaIndex] = IndivEtaIndex{:,:};

    if isempty(RespIndex) && ~isempty(IndivIndex)
        indivStartIndex = IndivIndex{end,end}(end)+1+n_cell;
        indivEndIndex =   indivStartIndex+n_rand-1;
        IndivOmegaIndex = mat2cell(reshape(indivStartIndex:...       % Create indexes for the variance parameters of
            indivEndIndex,n_rand,[]),ones(n_rand,1));         % individually-varying params starting at the end of
        
        
        H.SigmaParams = [H.SigmaParams [IndivOmegaIndex{:,:}]];                                       % Add noise variance indexes to the vector of variance indexes

    elseif ~isempty(RespIndex) && ~isempty(IndivIndex)
        indivStartIndex = RespIndex{end,end}(end)+1+n_cell;
        indivEndIndex =   indivStartIndex+n_rand-1;
        IndivOmegaIndex = mat2cell(reshape(indivStartIndex:...       % Create indexes for the variance parameters of
            indivEndIndex,n_rand,[]),ones(n_rand,1));         % individually-varying params starting at the end of
        
        
        H.SigmaParams = [H.SigmaParams [IndivOmegaIndex{:,:}]];                                       % Add noise variance indexes to the vector of variance indexes
    else
            IndivOmegaIndex ={};                                                    % If there are no individual parameters

    end
[H.IndividualParams(1:n_rand).OmegaIndex] = IndivOmegaIndex{:,:};          % Add omega indexes to H


%% Add response parameter indexes if available
[H.RespParams(1:n_resp).EtaIndex] = RespEtaIndex{:,:};
 if ~isempty(RespIndex) && isempty(IndivIndex)
        respStartIndex = RespIndex{end,end}(end)+1+n_cell;
        respEndIndex =   respStartIndex+n_resp-1;
        RespOmegaIndex = mat2cell(reshape(respStartIndex:...           % Create indexes for the variance parameters of
            respEndIndex,n_resp,[]),ones(n_resp,1));                    % response-varying params starting at the end of cell params
        
    elseif ~isempty(RespIndex) && ~isempty(IndivIndex)
        respStartIndex = IndivIndex{end,end}(end)+1+n_cell;
        respEndIndex =   indivStartIndex+n_resp-1;
        RespOmegaIndex = mat2cell(reshape(respStartIndex:...       % Create indexes for the variance parameters of
            respEndIndex,n_resp,[]),ones(n_resp,1));         % individually-varying params starting at the end of
        
        
    else
            RespOmegaIndex ={};                                                    % If there are no individual parameters

 end
H.SigmaParams = [H.SigmaParams [RespOmegaIndex{:,:}]];        % Add noise variance indexes to the vector of variance indexes

[H.RespParams(1:n_resp).OmegaIndex] =RespOmegaIndex{:,:};          % Add omega indexes to H

%% Add observables variances
try
 H.SigmaParams = [H.SigmaParams H.SigmaParams(end)+1:H.SigmaParams(end)+n_sigma];
catch
    H.SigmaParams = [H.PopulationParams(end)+1: H.PopulationParams(end)+n_sigma];
end