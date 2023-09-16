function DIC=getDIC(postSamples,logL,likelihood, varargin)
if nargin<3
    error('GWMCMC:toofewinputs','AMCMC requires atleast 3 inputs.')
end
p=inputParser;
p.addParameter('method','Gelman');

p.parse(varargin{:});
p=p.Results;
theta_hat=mean(postSamples);
DIC_theta_hat=2*likelihood(theta_hat);
DIC_hat= mean(2*logL);
if ismember(p.method,'Gelman')
p_D=var(2*logL)/2;
else
p_D=DIC_hat-DIC_theta_hat;
end
DIC=DIC_hat+p_D;
end
