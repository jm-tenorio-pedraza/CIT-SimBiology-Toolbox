function H = getHierarchicalStruct(params,  varargin)
p=inputParser;
p.addParameter('n_sigma',1,@isnumeric);
p.addParameter('n_rand',0,@isnumeric);
p.addParameter('n_indiv',1,@isnumeric);
p.parse(varargin{:})
p=p.Results;

n_sigma=p.n_sigma;
n_rand = p.n_rand;
n_indiv = p.n_indiv;

H.PopulationParams = 1:length(params);
 [H.IndividualParams(1:n_rand).name] = params{1:n_rand};
    Index=mat2cell(reshape(length(params)+1:length(params)+n_indiv*n_rand,[],n_rand)',ones(n_rand,1));
    EtaIndex = mat2cell(reshape(1:n_rand,n_rand,[]),ones(n_rand,1));
    [H.IndividualParams(1:n_rand).Index] = Index{:,:};
    [H.IndividualParams(1:n_rand).EtaIndex] = EtaIndex{:,:};
try
   
    OmegaIndex = mat2cell(reshape(Index{end,end}(end)+1:Index{end,end}(end)+n_rand,n_rand,[]),ones(n_rand,1)); 
    H.SigmaParams = [OmegaIndex{:,:}];
    H.SigmaParams = [H.SigmaParams H.SigmaParams(end)+1:H.SigmaParams(end)+n_sigma];
    
catch
    OmegaIndex ={};
    H.SigmaParams = length(params)+1:length(params)+n_sigma;
    
end
    [H.IndividualParams(1:n_rand).OmegaIndex] = OmegaIndex{:,:};

end