function [g] = getGradient(p,obj_fun,varargin)
par=inputParser;
par.addParameter('delta',1e-6)
par.addParameter('alpha',10)
par.addParameter('parallel',true)

par.parse(varargin{:})
par=par.Results;

[delta,~] = deal(par.delta,par.alpha);
n=length(p);
if size(p,1)> size(p,2)
else
    p=p';
end
p_delta_up = (exp(p) +delta*diag(ones(n,1)).*exp(p));
p_delta_down = (exp(p) -delta*diag(ones(n,1)).*exp(p));
if par.parallel
    g_delta_up = nan(n,1);
    g_delta_down = nan(n,1);
    try
    parfor i=1:n
        g_delta_up(i) = obj_fun(log(p_delta_up(:,i)));
        g_delta_down(i) = obj_fun(log(p_delta_down(:,i)));
    end
    catch
        pause(60*30)
        parfor i=1:n
        g_delta_up(i) = obj_fun(log(p_delta_up(:,i)));
        g_delta_down(i) = obj_fun(log(p_delta_down(:,i)));
        end
    end
else
g_delta_up = rowfun(obj_fun,table(log(p_delta_up')),'OutputFormat','uniform');
g_delta_down = rowfun(obj_fun,table(log(p_delta_down')),'OutputFormat','uniform');
end
g = (g_delta_up-g_delta_down)'./(log(diag(p_delta_up)) - log(diag(p_delta_down)))';
% g=nan(1,n);
% for i=1:n
%     g(1,i) = ( fun_val-obj_fun(p_delta(i,:)))/delta;
% end


