function PI=mcmcCI(PI, p_hat, logP, prob, varargin)
if nargin<3
    error('GWMCMC:toofewinputs','AMCMC requires atleast 3 inputs.')
end
p=inputParser;
p.addParameter('method','HPD');

p.parse(varargin{:});
p=p.Results;

n_p = size(PI.par,2);
n_samples= size(p_hat,1);
p_sorted=sortrows([logP p_hat],1);
p_sorted=p_sorted(end:-1:1,:);
indx=cumsum(p_sorted(:,1))/sum(p_sorted(:,1))>=prob;
int =1:n_samples;
indx=int(indx);
indx=indx(1);
for i=1:n_p
    modeHPD = mode(p_hat(:,i));
    PI.par(i).posterior_mode=p_sorted(1,i+1);
    if strcmp(p.method, 'symmetric')
        LB=quantile(p_hat(:,i),(1-prob)/2);
        UB=quantile(p_hat(:,i),0.5+prob/2);

        PI.par(i).posterior_mode=modeHPD;
        PI.par(i).LB=LB;
        PI.par(i).UB=UB;
    elseif strcmp(p.method, 'HPD')
        
        p_credible_set=p_sorted(p_sorted(:,1)>=p_sorted(indx,1),i+1);
        PI.par(i).LB=min(p_credible_set);
        PI.par(i).UB=max(p_credible_set);
        
    else
    end
end
end
