function [x, p_x,accept] = dreamHParallel(X0,likelihood,prior,N,T,d,varargin)
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
Xp = nan(N,d);
R = nan(N,N-1); p_X = nan(N,1);

[J,n_id] = deal(zeros(1,n_CR));                             % Variables selection prob. crossover
for i=1:N, R(i, 1:N-1) = setdiff(1:N,i); end                % R-matrix: index of chains for DE
CR = (1:n_CR)/n_CR; pCR = ones(1,n_CR)/n_CR;                % Crossover values and select. prob.

X = X0;                                                     % Create initial  population
for i = 1:N, p_X(i,1) = posterior(X(i, :)); end          % Compute density of initial population
x(1:d, 1:N, 1) = X'; p_x(1,1:N) = p_X';                     % Store initial states and density

H = p.H;
try
    ind_index = [H.IndividualParams(:).Index];
    param_index = setdiff(1:d, H.IndividualParams.Index);
catch
    param_index = 1:d;
end


D = reshape(randsample(1:delta, N*T, true),T,N);            % Select delta (equal selection probability)
z = rand(N*T,d);                                            % Draw d values from U[0,1]
U_ind = log(rand(T,N));                                     % Draw N values from U(0,1)
U_pop = log(rand(T,N));                                     % Draw N values from U(0,1)
lambda = unifrnd(-c, c, T,N);                               % Draw N lambda values
Draw = rand(N-1,N*T);
accept_pop = zeros(N,1);
accept_ind = zeros(N,1);

for t = 2:T
    z_t = z((t-1)*N+1:(t-1)*N+N,:);
    [~, draw] = sort(Draw(:,(t-1)*N+1:(t-1)*N+N));          % Permute [1, ..., N-1] N times
    dX = zeros(N,d);                                        % Set N jump vectors to zero
    std_X = std(X);                                         % Compute std each dimension
    Xp = X;
    id = randsample(1:n_CR, N, true, pCR);                  % Select index of crossover value
    A = z_t < CR(id);                       % Derive subset A with selected dimensions to be updated
    d_star = sum(A,2);                                  % How many dimensions sampled?
    
    if sum(d_star == 0)
        A(d_star==0,:) = ismember(z_t(d_star==0,:),min(z_t(d_star==0,:)));
        d_star(d_star==0) = 1;
    end                                 % A must contain at least 1 value
    gamma_d = stepSize./sqrt(2*D(t,:)'.*d_star);           % Calculate jump rate
    
    for i = 1:N                                             % Create proposals for individual parameters
        a = R(i, draw(1:D(t,i), i));
        b = R(i, draw(D(t,i)+1:2*D(t,i),i));                % Extract vectors a and b whose entries are indexes not equal to i
        % Calculate step size and update only those columns in A
        g=randsample([gamma_d(i) 1], 1, true, [1-p_g p_g]);    % Select gamma: 80/20 mix [default 1]
        
        dX(i,A(i,:)) = c_star*randn(1, d_star(i)) + ...
            (1+lambda(t,i))*g*sum(X(a,A(i,:))-X(b,A(i,:)),1);         % Compute ith jump diff. evol.
        J(id(i)) = J(id(i)) + sum((dX(i, 1:d)./std_X).^2);  % Update jump distance crossover idx
        n_id(id(i)) = n_id(id(i)) + 1;                      % How many times idx crossover used
    end
    
    Xp(:,ind_index) = Xp(:,ind_index)...
        + dX(:,ind_index);                          % Compute ith proposal for the individual parameters
    if ~isempty(H.IndividualParams)
        parfor i=1:N
            % Check if there are individual parameters to estimate
            % Evaluate proposal of new individual parameters
            [X(i,:), p_X(i,1), accept_i]= mcmcstep(X(i,:),...
                Xp(i,:), p_X(i,1), likelihood, prior,...
                U_ind(i), accept_ind(i));                      % Calculate pdf for ith proposal from individual parameters
            if accept_i > accept_ind(i)                        % MH criterion
                accept_ind(i) = accept_ind(i) + 1;
            else
                dX(i,ind_index) = 0;    % Set jump back to 0 for pCR for the individual parameters
            end
        end
    end
    
    Xp(:,param_index) = Xp(:,param_index)...
        + dX(:,param_index);      % Compute ith proposal for the individual parameters
    parfor i=1:N
        [X(i,:), p_X(i,1), accept_i]= mcmcstep(X(i,:),...
            Xp(i,:), p_X(i,1), likelihood, prior,...
            U_pop(i), accept_pop(i));
        
        if accept_i > accept_pop(i)                            % MH criterion for the population parameters
            accept_pop(i) = accept_pop(i) + 1;
        else
            dX(i,param_index) = 0;                          % Set jump back to 0 for pCR for population parameters
        end
    end
    totcount = N*(t);
    
    progress((t-1)/T,mean(X)',...
        mean([sum(accept_pop),sum(accept_ind)])/totcount)             % Print out progress status
    
    x(1:d, 1:N, t) = X'; p_x(t, 1:N) = p_X';                % Append current X and density
    
    if BurnIn>t*N
        if (sum(J)>0), pCR = J./n_id;
            pCR = pCR/sum(pCR); end                         % update selection prob. crossover
        [X, p_X] = check(X, mean((p_x(ceil(t/2):t,1:N))),p_X);       % Outlier detection and correction
    end
end

accept = mean([sum(accept_pop),sum(accept_ind)])/totcount;

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
