function [x, p_x, stepSize, J, n_id,accept] = par_dream_1(X,prior, likelihood, N, T,d, varargin)
par = inputParser;
par.addParameter('delta', 3);
par.addParameter('c', 0.1);
par.addParameter('c_star', 1e-12);
par.addParameter('n_CR', 3);
par.addParameter('p_g', 0.2);
par.addParameter('stepSize', 2.38);
par.addParameter('burnIn', 0);
par.addParameter('ProgressBar', true);
par.addParameter('J',[] );
par.addParameter('n_id',[]);
par.addParameter('plotLogL',false);
par.addParameter('H', []);
par.parse(varargin{:});
par=par.Results;

[delta, c, c_star, n_CR, p_g, stepSize,burnIn,H] = deal(par.delta, par.c, par.c_star,...
    par.n_CR, par.p_g, par.stepSize, par.burnIn,par.H);
if par.ProgressBar
    progress=@textprogress;
else
    progress=@noaction;
end

x = nan(d,N,T); p_x = nan(T,N);
CR = [1:n_CR]/n_CR; 
if isempty(par.J)
[J, n_id] = deal(ones(1,n_CR));
p_CR = ones(1,n_CR)/n_CR;
else
    J = par.J;
    n_id = par.n_id;
    p_CR = (J./n_id)/sum(J./n_id);
end
for i=1:N, R(i, 1:N-1) = setdiff(1:N,i); end
if ~isempty(H)
%     popParams = setdiff([H.PopulationParams H.SigmaParams], [H.CellParams.EtaIndex...
%         H.IndividualParams.EtaIndex H.RespParams.EtaIndex H.CellParams.Index H.IndividualParams.Index H.RespParams.Index]);
    IndivParamsIndx = ismember(1:d, [H.CellParams.EtaIndex...
        H.IndividualParams.EtaIndex H.RespParams.EtaIndex H.CellParams.Index H.IndividualParams.Index H.RespParams.Index]);
    PopParamsIndx= ismember(1:d, [H.PopulationParams H.SigmaParams]);
else
    PopParamsIndx = ismember(1:d,1:d);
    IndivParamsIndx = ismember(1:d,1:d);
end

%%X = prior(N,d);
for i=1:N, p_X(i, 1) = prior(X(i, 1:d))+likelihood(X(i, 1:d)); end
x(1:d, 1:N,1) = X'; p_x(1,1:N) = p_X';
X_p = X;
accept = zeros(N,1);
p=parpool('local');
pctRunOnAll warning off
for t = 2:T
    [~, draw] = sort(rand(N-1, N));
    dX = zeros(N, d);
    lambda=(unifrnd(-c, c, N,1));
    std_X = std(X);
    id = randsample(1:n_CR, N, true, p_CR);
    D = randsample([1:delta], N, true);
    z = rand(N, d);
    A = (z < CR(id)');
    walkerIndx = randsample([1:N], N,false);
    A(walkerIndx(1:sum(IndivParamsIndx)),:) = repmat(IndivParamsIndx,sum(IndivParamsIndx),1);
    A(walkerIndx(sum(IndivParamsIndx)+1:N),:) = A(walkerIndx(sum(IndivParamsIndx)+1:N),:).*PopParamsIndx;

    for i=1:N
%         if any(A(i,:).*IndivParamsIndx)
%             A(i,IndivParamsIndx) = 1;
%         end
        d_star = sum(A(i,:));
        if (d_star == 0)
            [~,minZindx]= min(z(i,:));
            A(i,:) = ismember(1:d,minZindx);
            d_star = 1;
        end
        gamma_d = stepSize/sqrt(2*D(i).*d_star');
        g = randsample([gamma_d 1], 1, true, [1-p_g, p_g]);
        a = R(i, draw(1:D(i), i)); b = R(i, draw(D(i)+1:2*D(i), i));
        dX(i,A(i,:)) = c_star*randn(1, d_star) + (1+lambda(i))*g.*sum(X(a,A(i,:))-X(b,A(i,:)),1);
        X_p(i,1:d) = X(i, 1:d) + dX(i,1:d);
    end
    logU = log(rand(1,N));
    try
    parfor i=1:N
        p_Xp(i,1) = prior(X_p(i, 1:d));
        if isinf(p_Xp(i,1)) || isnan(p_Xp(i,1)) || ~isreal(p_Xp(i,1))
            dX(i,:) =0;
            continue
        else
            p_Xp(i,1) = p_Xp(i,1) + likelihood(X_p(i,1:d));
            p_acc = min(0, p_Xp(i,1)-p_X(i,1));
      
             if p_acc > logU(i)
                 X(i, :) = X_p(i, 1:d); p_X(i,1) = p_Xp(i,1);
                 accept(i,:) = accept(i,:) + 1;
             else
                 dX(i,:) = 0;
             end
        end
    end
    catch
        delete(p)
        pause(30*60)
        p=parpool('local');
        pctRunOnAll warning off
    end

    x(1:d,1:N,t) = X'; p_x(t, 1:N) = p_X';
    if  t*N<burnIn
        [X, p_X] = check(X, mean((p_x(ceil(t/2):t,1:N))),p_X,'dxN',false);
        for i=1:N
            J(id(i)) = J(id(i)) + sum((dX(i,1:d)./std_X).^2);
            n_id(id(i)) = n_id(id(i)) + 1;
        end
        p_CR = J./n_id; p_CR = p_CR/sum(p_CR); 
%         if mod(t, 100)==0
%             if min(accept)/(t) < 0.2
%                 stepSize = stepSize*0.9;
%             elseif max(accept)/(t) > 0.65
%                 stepSize = stepSize*1.1;
%             end
%         end
    end
    progress((t-1)/T,mean(X)',(sum(accept)/(N*t)))            % Print out progress status
    if mod(t,100)==0 && par.plotLogL
    plot(p_x)
    hold on
    plot(mean(p_x,2), 'LineWidth', 4)
    ylim([quantile(p_x(t,:), 0.02) quantile(p_x(t,:), 0.97)])
    end
    
end
delete(p)

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