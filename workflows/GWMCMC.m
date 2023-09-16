p=parpool('local');
pctRunOnAll warning off
w0 = x12(:,:,end);
[x1, p_x1] = gwmcmc_edit(w0,{prior_fun likelihood_fun},3e5,'Parallel',true,'ThinChain',1);
