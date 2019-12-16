function [x, p_x,accept,pCR] = dreamH(X0,likelihood,prior,N,T,d,varargin)
% Inputs
% prior = handle to generate initial random samples
% pdf = handle to log-posterior function
% N = numer of chains (walkers)
% T = numger of generations
% d = parameter dimension
% Outputs
% x = samples from posterior distribution
% p_x = log-posterio probabilities for each sample in x
% Default parameters
p=inputParser;
p.addParameter('ThinChain',1,@isnumeric);
p.addParameter('ProgressBar',true,@islogical);
p.addParameter('delta',3,@isnumeric); %addParamValue is chose for compatibility with octave. Still Untested.
p.addParameter('c',0.1,@isnumeric); %addParamValue is chose for compatibility with octave. Still Untested.
p.addParameter('c_star',1e-12,@isnumeric); %addParamValue is chose for compatibility with octave. Still Untested.
p.addParameter('n_CR',3,@isnumeric); %addParamValue is chose for compatibility with octave. Still Untested.
p.addParameter('p_g',0.2,@isnumeric); %addParamValue is chose for compatibility with octave. Still Untested.
p.addParameter('StepSize',2.38,@isnumeric); %addParamValue is chose for compatibility with octave. Still Untested.
p.addParameter('H',[],@isstruct); % Hierarchical structure with Individual, Population and Sigma parameters
p.addParameter('BurnIn',0,@(x)(x>=0));
p.addParameter('Parallel',true, @islogical);
p.parse(varargin{:});
p=p.Results;

if p.ProgressBar
    progress=@textprogress;
else
    progress=@noaction;
end
BurnIn=p.BurnIn;
delta = p.delta;                                            % Default parameter values
c = p.c;
c_star = p.c_star;
n_CR = p.n_CR;
p_g = p.p_g;
stepSize=p.StepSize;
posterior = @(x)likelihood(x) + prior(x);
x=nan(d,N,T); p_x = nan(T,N);                               % Preallocate chains and density
R = nan(N,N-1); p_X = nan(N,1);

[J,n_id] = deal(zeros(1,n_CR));                             % Variables selection prob. crossover
for i=1:N, R(i, 1:N-1) = setdiff(1:N,i); end                % R-matrix: index of chains for DE
CR = (1:n_CR)/n_CR; pCR = ones(1,n_CR)/n_CR;                % Crossover values and select. prob.

X = X0;                                                     % Create initial  population
for i = 1:N, p_X(i,1) = posterior(X(i, :)); end          % Compute density of initial population
x(1:d, 1:N, 1) = X'; p_x(1,1:N) = p_X';                     % Store initial states and density

H = p.H;
try
    ind_index = [H.CellParams(:).Index H.IndividualParams(:).Index];
    param_index = setdiff(1:d, ind_index);
catch
    param_index = 1:d;
end


D = reshape(randsample(1:delta, N*T, true),T,N);            % Select delta (equal selection probability)
z = rand(N*T,d);                                            % Draw d values from U[0,1]
U_pop = log(rand(T,N));                                     % Draw N values from U(0,1)
U_ind = log(rand(T,N));                                     % Draw N values from U(0,1)
lambda = unifrnd(-c, c, T,N);                               % Draw N lambda values

accept_pop = zeros(N,1);
accept_ind = zeros(N,1);

for t = 2:T
    [~, draw] = sort(rand(N-1,N));                          % Permute [1, ..., N-1] N times
    dX = zeros(N,d);                                        % Set N jump vectors to zero
    Xp = X;
    std_X = std(X);                                         % Compute std each dimension
    id = randsample(1:n_CR, N, true, pCR);                  % Select index of crossover value
    % Compute ith proposal for the individual parameters
    for i=1:N
        a = R(i, draw(1:D(t,i), i));                        % Calculate step size and update only those columns in A
        b = R(i, draw(D(t,i)+1:2*D(t,i),i));                % Extract vectors a and b whose entries are indexes not equal to i
        A = find(z((t-1)*N+i,:) < CR(id(i)));               % Derive subset A with selected dimensions to be updated
        d_star = numel(A);                                  % How many dimensions sampled?
        
        if d_star == 0, [~, A] = min(z((t-1)*N+i,:));
            d_star = 1; end                                 % A must contain at least 1 value
        
        gamma_d = stepSize/sqrt(2*D(t,i)*d_star);           % Calculate jump rate
        g=randsample([gamma_d 1], 1, true, [1-p_g p_g]);    % Select gamma: 80/20 mix [default 1]
        
        dX(i,A) = c_star*randn(1, d_star) + ...
            (1+lambda(t,i))*g*sum(X(a,A)-X(b,A),1);         % Compute ith jump diff. evol.
       
        if ~isempty(H.IndividualParams)                                     % Check if there are individual parameters to estimate
            Xp(:,ind_index) = Xp(:,ind_index)...                            % Evaluate proposal of new individual parameters
                + dX(:,ind_index);
            for k=1:N
                [X(k,:), p_X(k,1), accept_i]= mcmcstep(X(k,:),...
                    Xp(k,:), p_X(k,1), likelihood, prior,...
                    U_ind(t,k), accept_ind(k));                             % Calculate pdf for ith proposal from individual parameters
                if accept_i > accept_ind(k)                                 % MH criterion
                    accept_ind(k) = accept_ind(k) + 1;
                else
                    dX(k,ind_index) = 0;                                    % Set jump back to 0 for pCR for the individual parameters
                end
            end
        end
         Xp(i,param_index) = Xp(i,param_index)...
            + dX(i,param_index);
        
        [X(i,:), p_X(i,1), accept_i]= mcmcstep(X(i,:),...
            Xp(i,:), p_X(i,1), likelihood, prior,...
            U_pop(t,i), accept_pop(i));
        if accept_i > accept_pop(i)                                         % MH criterion for the population parameters
            accept_pop(i) = accept_pop(i) + 1;
        else
            dX(i,param_index) = 0;                                          % Set jump back to 0 for pCR for population parameters
        end
        J(id(i)) = J(id(i)) + sum((dX(i,:)./std_X).^2);                     % Update jump distance crossover idx
        n_id(id(i)) = n_id(id(i)) + 1;
    end
    totcount = N*(t);
   if ~isempty(H.IndividualParams)
        accept = mean([sum(accept_pop),sum(accept_ind)])/totcount;
    else
        accept = sum(accept_pop)/totcount;
   end
    
    progress((t-1)/T,mean(X)',...
        accept)                                                 % Print out progress status
    
    x(1:d, 1:N, t) = X'; p_x(t, 1:N) = p_X';                    % Append current X and density
    
    if BurnIn>t*N
        if (sum(J)>0), pCR = J./n_id;
            pCR = pCR/sum(pCR); end                             % update selection prob. crossover
        [X, p_X] = check(X, mean((p_x(ceil(t/2):t,1:N))),p_X);  % Outlier detection and correction
        
        if accept<0.2
            stepSize = stepSize*0.67;
        elseif accept>0.4
            stepSize = stepSize*1.33;
        end
    end
end


function textprogress(pct,curm,acceptpct)
persistent lastNchar lasttime starttime
if isempty(lastNchar)||pct==0
    lasttime=cputime-10;starttime=cputime;lastNchar=0;
    pct=1e-16;
end
if pct==1
    fprintf('%s',repmat(char(8),1,lastNchar));lastNchar=0;
    return
end
if (cputime-lasttime>0.1)
    
    ETA=datestr((cputime-starttime)*(1-pct)/(pct*60*60*24),13);
    progressmsg=[183-uint8((1:40)<=(pct*40)).*(183-'*') ''];
    curmtxt=sprintf('% 9.3g\n',curm(1:min(end,20),1));
    progressmsg=sprintf('\nDREAM %5.1f%% [%s] %s\n%3.0f%% accepted\n%s\n',...
        pct*100,progressmsg,ETA,acceptpct*100,curmtxt);
    
    fprintf('%s%s',repmat(char(8),1,lastNchar),progressmsg);
    drawnow;lasttime=cputime;
    lastNchar=length(progressmsg);
end

function noaction(varargin)
