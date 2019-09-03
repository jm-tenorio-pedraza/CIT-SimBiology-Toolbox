function [x, p_x,accept] = dream(X0,pdf,N,T,d,varargin)
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

p.addParameter('BurnIn',0,@(x)(x>=0));

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

x=nan(d,N,T); p_x = nan(T,N);                               % Preallocate chains and density
Xp = nan(N,d); p_Xp = nan(N,1);
[J,n_id] = deal(zeros(1,n_CR));                             % Variables selection prob. crossover
for i=1:N, R(i, 1:N-1) = setdiff(1:N,i); end                % R-matrix: index of chains for DE
CR = (1:n_CR)/n_CR; pCR = ones(1,n_CR)/n_CR;                % Crossover values and select. prob.

X = X0;                                                     % Create initial  population
for i = 1:N, p_X(i,1) = pdf(X(i, 1:d)); end                 % Compute density of initial population
x(1:d, 1:N, 1) = X'; p_x(1,1:N) = p_X';                      % Store initial states and density

accept = 0;
for t = 2:T
    [~, draw] = sort(rand(N-1,N));                          % Permute [1, ..., N-1] N times
    dX = zeros(N,d);                                        % Set N jump vectors to zero
    lambda = unifrnd(-c, c, N,1);                           % Draw N lambda values
    std_X = std(X);                                         % Compute std each dimenstion
    for i = 1:N                                             % Create proposals and accept/reject
        D = randsample(1:delta, 1, true);                   % Select delta (equal selection probability)
        a = R(i, draw(1:D, i)); b = R(i, draw(D+1:2*D,i));  % Extract vectors a and b whose entries are indexes not equal to i
        id = randsample(1:n_CR, 1, true, pCR);              % Select index of crossover value
        z = rand(1,d);                                      % Draw d values from U[0,1]
        A = find(z < CR(id));                               % Derive subset A selected dimensions to be updated
        d_star = numel(A);                                  % How many dimensions sampled?
        if d_star == 0, [~, A] = min(z); d_star = 1; end    % A must contain at least 1 value
        gamma_d = stepSize/sqrt(2*D*d_star);                     % Calculate jump rate
        g=randsample([gamma_d 1], 1, true, [1-p_g p_g]);    % Select gamma: 80/20 mix [default 1]
        dX(i,A) = c_star*randn(1, d_star) + ...             
            (1+lambda(i))*g*sum(X(a,A)-X(b,A),1);           % Compute ith jump diff. evol.
        Xp(i,1:d) = X(i,1:d) + dX(i,1:d);                   % Compute ith proposal
        p_Xp(i,1) = pdf(Xp(i, 1:d));                        % Compute density of ith proposal
        p_acc = min(0,p_Xp(i,1)-p_X(i,1));                  % Compute acceptance prob.
        if p_acc > log(randn)                               % MH criterion
            X(i,1:d) = Xp(i,1:d); p_X(i,1) = p_Xp(i,1);
            accept = accept +1;
        else
            dX(i,1:d) = 0;                                  % Set jump back to 0 for pCR
        end
        J(id) = J(id) + sum((dX(i, 1:d)./std_X).^2);        % Update jump distance crossover idx
        n_id(id) = n_id(id) + 1;                            % How many times idx crossover used
        totcount = N*(t-1)+i;

    end
    progress((t-1)/T,Xp(i,:)',(accept)/(totcount))            % Print out progress status

    x(1:d, 1:N, t) = X'; p_x(t, 1:N) = p_X'; % Append current X and density
    if BurnIn>t*N
        if (sum(J)>0), pCR = J./n_id; pCR = pCR/sum(pCR); end           % update selection prob. crossover

        [X, p_X] = check(X, mean((p_x(ceil(t/2):t,1:N))),p_X);       % Outlier detection and correction
    else
    end
end


accept = accept/totcount;

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
 