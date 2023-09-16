function [x, p_x, stepSize, J, n_id,accept,Z] = dream_zs(Z0,prior, likelihood, N, T,d, varargin)
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
par.addParameter('m0',size(Z0,1));
par.addParameter('k',10);
par.addParameter('p_z',.1);

par.parse(varargin{:});
par=par.Results;

[delta, c, c_star, n_CR, p_g, stepSize,burnIn,H,m0,k,p_z] = deal(par.delta, par.c, par.c_star,...
    par.n_CR, par.p_g, par.stepSize, par.burnIn,par.H,par.m0,par.k,par.p_z);
if par.ProgressBar
    progress=@textprogress;
else
    progress=@noaction;
end

x = nan(d,N,T); p_x = nan(T,N);
Z = nan(m0+N*floor(T/k),d);
Z(1:m0,:) = Z0;
Zindx= m0+1;
m=m0;
X = Z0(end-N+1:end,:);
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
p_X=nan(N,1);
for i=1:N, p_X(i, 1) = prior(X(i, 1:d))+likelihood(X(i, 1:d)); end
x(1:d, 1:N,1) = X'; p_x(1,1:N) = p_X';
X_p = X;
accept = zeros(N,1);

Delta = kron(ones(delta,1),diag(ones(N,1)));
IDrows = ones(n_CR,N).*[1:n_CR]';
for t = 2:T
    draw = (randsample(1:m,N*delta*2,false));
    std_X = std(X);
    id = randsample(1:n_CR, N, true, p_CR);
    ID= (IDrows == id)';
%     dX = zeros(N, d);
    zprob=randsample([0,1],1,true, [1-p_z, p_z]); % Prob. of a snooker update
    if zprob==0
        lambda=(unifrnd(-c, c, N,1));
        z = rand(N, d);
        A = (z < CR(id)');
        d_star = sum(A,2);
        if any(d_star == 0)
            [~,minZindx]= min(z(d_star == 0,:));
            A(d_star == 0,minZindx) = 1;
            d_star(d_star == 0) = 1;
        end
            gamma_d = stepSize./sqrt(2*delta*d_star');
            g = ones(N,1);
            jumpProb = randsample([0 1], N, true, [1-p_g, p_g]);
            g(jumpProb<1) = gamma_d(jumpProb<1);
            Za = Z(draw(1:delta*N),:)'*Delta;
            Zb = Z(draw(delta*N+1:end),:)'*Delta;
            dX = A.*(c_star*randn(N,d) + (1+lambda).*g.*(Za-Zb)');

    else
            g = rand(N,1)*(2.2-1.2)+1.2;
            za= X - Z(draw(1:N),:);
            za_norm =za./diag(za*za'); 
            zb = Z(draw(N+1:N*2),:) - diag(za*Z(draw(N+1:N*2),:)').*za_norm;
            zc = Z(draw(N*2+1:N*3),:) - diag(za*Z(draw(N*2+1:N*3),:)').*za_norm;
            dX = c_star*randn(N, d) + g.*(zb-zc);

    end
    X_p = X + dX;
    logU = log(rand(1,N));
    try
    for i=1:N
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
    end

    x(1:d,1:N,t) = X'; p_x(t, 1:N) = p_X';
    dXstd = sum((dX./std_X).^2,2);
    if  t*N<burnIn
        J = J+(sum((dX./std_X).^2,2)'*ID);
        n_id = n_id+sum(ID);
%         for i=1:N
%             J(id(i)) = J(id(i)) + sum((dX(i,1:d)./std_X).^2);
%             n_id(id(i)) = n_id(id(i)) + 1;
%         end
        p_CR = J./n_id; p_CR = p_CR/sum(p_CR); 

    end
    if mod(t,k)==0
        Z(Zindx:Zindx+N-1,:) = X;
        Zindx=Zindx+N;
        m=Zindx-1;
    end
    progress((t-1)/T,mean(X)',(sum(accept)/(N*t)))            % Print out progress status
    if mod(t,100)==0 && par.plotLogL
    plot(p_x)
    hold on
    plot(mean(p_x,2), 'LineWidth', 4)
%     ylim([quantile(p_x(t,:), 0.02) quantile(p_x(t,:), 1)])
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