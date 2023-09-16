function [x,pdf_x] = hmcmc(logpdf,q_curr,N,varargin)
par=inputParser;
par.addParameter('epsilon',0.01);
par.addParameter('L',1);
par.addParameter('delta',0.000001);
par.addParameter('ProgressBar', true);
par.addParameter('plotLogL',true);
par.parse(varargin{:});
par=par.Results;
if par.ProgressBar
    progress=@textprogress;
else
    progress=@noaction;
end
[epsilon,L,delta] = deal(par.epsilon,par.L,par.delta);
d = length(q_curr);
U_curr = logpdf(q_curr);
P = randn(N,d);
g_up = nan(1,d);
g_down = nan(1,d);

logUnif = log(rand(N*L,1));
x = nan(N,d);
pdf_x = nan(N,1);
accept=0;
for i=1:N
    q=q_curr;
    p = P(i,:);
    p_curr = p;
    q_delta_up = exp(q) +delta*diag(ones(d,1)).*exp(q);
    q_delta_down = exp(q)-delta*diag(ones(d,1)).*exp(q);
    parfor k=1:d
        g_up(k) = logpdf(log(q_delta_up(k,:)));
        g_down(k) = logpdf(log(q_delta_down(k,:)));
    end
        g = (g_up-g_down)./((diag(q_delta_up)'-q)*2);
        p = p-epsilon*g/2;
    
    %     for j=1:L
    %     q = q + epsilon*p;
    %     U = logpdf(q);
    %         q_delta_up = exp(q) +delta*diag(ones(d,1)).*exp(q);
    %         q_delta_down = exp(q)-delta*diag(ones(d,1)).*exp(q);
    %
    %         if j<L
    %             parfor k=1:d
    %                 g_up(k) = logpdf(log(q_delta_up(k,:)));
    %                 g_down(k) = logpdf(log(q_delta_down(k,:)));
    %
    %             end
    %     g = (g_up-g_down)./((diag(q_delta_up)'-q)*2);
    %             p = p-epsilon*g;
    %         end
    %     end
    %
    %     parfor k=1:d
    %         g_up(k) = logpdf(log(q_delta_up(k,:)))
    %         g_down(k) = logpdf(log(q_delta_down(k,:)))
    %     end    
    
  
    q = q + epsilon*p;
    U = logpdf(q);
    p = -p;
    
    K_curr = sum(p_curr.^2)/2;
    U_prop = U;
    K_prop = sum(p.^2)/2;
    
    p_acc = min(0, U_prop-U_curr +K_prop-K_curr);
    if logUnif(i)<(p_acc)
        
        pdf_x(i) = U_prop;
        U_curr= U_prop;
        
        accept = accept+1;
    else
        q = q_curr;
    end
    x(i,:) = q;
    pdf_x(i) = U_curr;
    if mod(i,20)==0 && par.plotLogL
        plot(pdf_x(1:i))
        hold on
        ylim([quantile(pdf_x(1:i), 0.02) quantile(pdf_x(1:i), 0.97)])
    end
    progress(i/N,x(i,:)',(accept/(i)))            % Print out progress status
    
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


