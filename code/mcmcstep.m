function [curr_p, curr_L, accept]=mcmcstep(curr_p, prop_p, curr_L, logL_fun, prior_fun,U, accept)
prior_L=prior_fun(prop_p);
if isinf(prior_L) || ~isreal(prior_L) || isnan(prior_L)
   
else
    prop_L=logL_fun(prop_p)+prior_L;
    r=(prop_L-curr_L);
    if U<=min(r,0)
        curr_L=prop_L;
        curr_p=prop_p;
        accept=accept+1;
    else
    end
        
end
end
