function [x, p_x,accept,pCR,stepSize,J, n_id] = par_H_dream(X0,prior,likelihood,N,T,d,varargin)
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
p.addParameter('stepSize',2.38);
p.addParameter('H', []);
p.addParameter('BurnIn',0);
p.addParameter('pCR',[]);
p.addParameter('J',[]);
p.addParameter('n_id',[]);
p.parse(varargin{:});
p=p.Results;

if p.ProgressBar
    progress=@textprogress;
else
    progress=@noaction;
end

[delta, c,c_star,n_CR, p_g,stepSize, BurnIn,H] = ...
    deal(p.delta, p.c, p.c_star, p.n_CR, p.p_g, p.stepSize,p.BurnIn, p.H);
if isempty(p.pCR)
    pCR = ones(1,n_CR)/n_CR;                % Crossover values and select. prob.
    [J,n_id] = deal(zeros(1,n_CR));                             % Variables selection prob. crossover

else
    [pCR,J,n_id] = deal(p.pCR, p.J, p.n_id);
end
x=nan(d,N,T); p_x = nan(T,N);                               % Preallocate chains and density
Xp = nan(N,d); p_Xp = nan(N,1);p_X = nan(N,1);
CR = (1:n_CR)/n_CR;                 % Crossover values and select. prob.
R = nan(N,N-1);
for i=1:N, R(i, 1:N-1) = setdiff(1:N,i); end                % R-matrix: index of chains for DE

X = X0; % dxN                                                     % Create initial  population
for i = 1:N, p_X(i,1) = likelihood(X(:,i))+prior(X(:,i)); end                 % Compute density of initial population
x(1:d, 1:N,1) =X0; p_x(1,1:N) = p_X';    % Store initial states and density

% Defining static inputs outside of loop
D = reshape(randsample(1:delta,N*T, true),T,N);                   % Select delta (equal selection probability)
z = rand(N*T,d);
U = log(rand(T,N));
accept = zeros(N,1);
lambda = unifrnd(-c, c,T, N);                           % Draw N lambda values
try
    if ~isempty(H.CellParams(1).Index) 
    popParamIndx = [H.PopulationParams H.CellParams.Index H.SigmaParams];
    indivParamIndx = [H.IndividualParams.Index];
    else
        popParamIndx = [H.PopulationParams  H.SigmaParams];
    indivParamIndx = [H.IndividualParams.Index];
    end
catch
    popParamIndx = 1:d;
    indivParamIndx = [];
end
for t = 2:T
    [~, draw] = sort(rand(N-1,N));                          % Permute [1, ..., N-1] N times
    dX = zeros(d,N);                                        % Set N jump vectors to zero
    std_X = std(X,[],2);                                         % Compute std each dimenstion
    id = randsample(1:n_CR,N, true, pCR);               % Select index of crossover value
    for i = 1:N                                             % Create proposals and accept/reject
        a = R(i, draw(1:D(t,i), i)); 
        b = R(i, draw(D(t,i)+1:2*D(t,i),i));  % Extract vectors a and b not equal to i
        A = find(z((t-1)*N+i,:) < CR(id(i)));                               % Derive subset A selected dimensions
        d_star = numel(A);                                  % How many dimensions sampled?
        if d_star == 0, [~, A] = min(z((t-1)*N+i,:)); d_star = 1; end    % A must contain at least 1 value
        gamma_d = stepSize/sqrt(2*D(t,i)*d_star);                     % Calculate jump rate
        g = randsample([gamma_d 1], 1, true, [1-p_g p_g]);    % Select gamma: 80/20 mix [default 1]
        dX(A,i) = c_star*randn(1, d_star)' + ...             
            (1+lambda(t,i))*g*sum(X(A,a)-X(A,b),2);           % Compute ith jump diff. evol.
    end
    Xp = X;
    Xp(popParamIndx,:) = X(popParamIndx,:) + dX(popParamIndx,:);                   % Compute ith proposal
    
    parfor i = 1:N
        prop_p = Xp(:,i);
        p_Xp_i= prior(prop_p);
        if isinf(p_Xp_i) || ~isreal(p_Xp_i) || isnan(p_Xp_i)
            dX(popParamIndx,i) = 0;
        else
           propL = likelihood(prop_p)+p_Xp_i;                        % Compute density of ith proposal
           r = propL-p_X(i,1);
           if U(t,i)<=min(r,0)
               p_X(i,1) = propL;
               X(popParamIndx,i) = prop_p(popParamIndx);
               accept(i,1)=accept(i,1)+1;
           else 
               dX(popParamIndx,i) = 0;
           end
        end
    end
    if ~isempty(indivParamIndx)
        Xp(popParamIndx,:) = X(popParamIndx,:);
        Xp(indivParamIndx,:) = X(indivParamIndx,:) + dX(indivParamIndx,:);                   % Compute ith proposal
        parfor i = 1:N
            prop_p = Xp(:,i);
            p_Xp_i= prior(prop_p);
            if isinf(p_Xp_i) || ~isreal(p_Xp_i) || isnan(p_Xp_i)
                dX(indivParamIndx,i) = 0;
            else
                propL = likelihood(prop_p)+p_Xp_i;                        % Compute density of ith proposal
                r = propL-p_X(i,1);
                if U(t,i)<=min(r,0)
                    p_X(i,1) = propL;
                    X(indivParamIndx,i) = prop_p(indivParamIndx);
                    accept(i,1)=accept(i,1)+1;
                else
                    dX(indivParamIndx,i) = 0;
                end
            end
        end
        
        
    end
for i=1:N 
    J(id(i)) = J(id(i)) +sum(dX(:,i)./std_X).^2;
    n_id(id(i)) = n_id(id(i))+1;
end
    x(1:d,1:N,t) = X; p_x(t, 1:N) = p_X'; % Append current X and density
    if t*N<BurnIn, pCR =J./n_id; 
        pCR = pCR/sum(pCR); 
        if mod(t,ceil(BurnIn/N/40))==0
        if mean(accept)/t<0.2
            stepSize=stepSize*(1-1/N);
        elseif mean(accept)/t>0.3
            stepSize=stepSize*(1+1/N);
        end
        end

    end           % update selection prob. crossover
    [X, p_X] = check(X, mean((p_x(ceil(t/2):t,1:N))),p_X);       % Outlier detection and correction
    
    progress((t-1)/T,mean(Xp,2),(mean(accept)/(t)))            % Print out progress status
    if mod(t,20)==0
    plot(p_x)
    hold on
    plot(mean(p_x,2), 'LineWidth', 4)
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
 