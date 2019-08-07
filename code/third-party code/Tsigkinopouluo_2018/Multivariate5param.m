function [Keqfin, KAfin, KPfin, V1fin, V2fin] = Multivariate5param(muKeq, sigmaKeq, mukmf, sigmakmf, mukmr, sigmakmr, muvmaxf, sigmavmaxf, muvmaxr, sigmavmaxr, ParamNo,Eo)
% The function generates a multivariate distribution for reactions with   
% Michaelis-Menten kinetics with 5 parameters. The user must provide the  
% mu and sigma for the 5 distributions in the correct order,as well as the  
% number of parameter samples that are needed (ParamNo).

%        k1    k3
% E + A <-> X <-> P + E
%        k2    k4
%        KD1   KD2

% STAGE 1: Parameter consistency
% This stage finds the parameter with the largest geometric coeeficient of
% variation and assigns it as dependent on the remaining 4 parameters. The
% distribution for the dependent parameter is calculated via the Haldane
% equation from the distributions of the 4 independent parameters.

% Calculate geometric coefficients of variation:
GCV1=exp(sigmaKeq)-1;
GCV2=exp(sigmakmf)-1;
GCV3=exp(sigmakmr)-1;
GCV4=exp(sigmavmaxf)-1;
GCV5=exp(sigmavmaxr)-1;

% Choose the parameter with the largest coefficient of variation:
A=[GCV1 GCV2 GCV3 GCV4 GCV5];
M=max(A);

% Set dependent parameter based on the largest GCV and 
% calculate new mu and sigma for the dependent distribution:

if M==GCV2  % kmf dependent
   mukmf=(muvmaxf+mukmr)-(muKeq+muvmaxr);
   sigmakmf=sqrt(sigmavmaxf^2+sigmakmr^2+sigmaKeq^2+sigmavmaxr^2);
   
elseif M==GCV1 % Keq dependent
   muKeq=(muvmaxf+mukmr)-(mukmf+muvmaxr);
   sigmaKeq=sqrt(sigmavmaxf^2+sigmakmr^2+sigmakmf^2+sigmavmaxr^2);
   
elseif M==GCV3 % kmr dependent
   mukmr=(muKeq+muvmaxr+mukmf)-muvmaxf;
   sigmakmr=sqrt(sigmaKeq^2+sigmavmaxr^2+sigmakmf^2+sigmavmaxf^2);

elseif M==GCV4 % vmaxf dependent
   muvmaxf=(muKeq+muvmaxr+mukmf)-mukmr;
   sigmavmaxf=sqrt(sigmaKeq^2+sigmavmaxr^2+sigmakmf^2+sigmakmr^2);

elseif M==GCV5 % vmaxr dependent
   muvmaxr=(muvmaxf+mukmr)-(muKeq+mukmf);
   sigmavmaxr=sqrt(sigmavmaxf^2+sigmakmr^2+sigmaKeq^2+sigmakmf^2);
end

% STAGE 1: Parameter correlation
% This stage calculates the correlation between the external
% michaelis-menten parameters, by using the internal mass action kinetics.

% Design mass action kinetics' distributions:
% k1=(V1+V2)./(KA.*Eo);
sigmak1=sqrt(log(((exp(2*muvmaxf+sigmavmaxf^2)*(exp(sigmavmaxf^2)-1)+exp(2*muvmaxr+sigmavmaxr^2)*(exp(sigmavmaxr^2)-1))/((exp(muvmaxf+(sigmavmaxf^2)/2)+exp(muvmaxr+(sigmavmaxr^2)/2))^2))+1)+sigmakmf^2);
muk1=log(exp(muvmaxf+(sigmavmaxf^2)/2)+exp(muvmaxr+(sigmavmaxr^2)/2))-(log(((exp(2*muvmaxf+sigmavmaxf^2)*(exp(sigmavmaxf^2)-1)+exp(2*muvmaxr+sigmavmaxr^2)*(exp(sigmavmaxr^2)-1))/((exp(muvmaxf+(sigmavmaxf^2)/2)+exp(muvmaxr+(sigmavmaxr^2)/2))^2))+1))/2-(mukmf+log(Eo));

% k2=V2./Eo;
muk2=muvmaxr-log(Eo);
sigmak2=sigmavmaxr;

% k3=V1./Eo;
muk3=muvmaxf-log(Eo);
sigmak3=sigmavmaxf;

% k4=(V1+V2)./(KP.*Eo);
sigmak4=sqrt(log(((exp(2*muvmaxf+sigmavmaxf^2)*(exp(sigmavmaxf^2)-1)+exp(2*muvmaxr+sigmavmaxr^2)*(exp(sigmavmaxr^2)-1))/((exp(muvmaxf+(sigmavmaxf^2)/2)+exp(muvmaxr+(sigmavmaxr^2)/2))^2))+1)+sigmakmr^2);
muk4=log(exp(muvmaxf+(sigmavmaxf^2)/2)+exp(muvmaxr+(sigmavmaxr^2)/2))-(log(((exp(2*muvmaxf+sigmavmaxf^2)*(exp(sigmavmaxf^2)-1)+exp(2*muvmaxr+sigmavmaxr^2)*(exp(sigmavmaxr^2)-1))/((exp(muvmaxf+(sigmavmaxf^2)/2)+exp(muvmaxr+(sigmavmaxr^2)/2))^2))+1))/2-(mukmr+log(Eo));

%KD1=k2./k1;
muKD1=muk2-muk1;
sigmaKD1=sqrt(sigmak2^2+sigmak1^2);

%KD2=k4./k3;
muKD2=muk4-muk3;
sigmaKD2=sqrt(sigmak4^2+sigmak3^2);

% Generate random values from the distributions of KD1, KD2, k2 and k3:
pd6=makedist('Lognormal',muKD1,sigmaKD1);
KD1=random(pd6,1000000,1);
pd7=makedist('Lognormal',muKD2,sigmaKD2);
KD2=random(pd7,1000000,1);
pd9=makedist('Lognormal',muk2,sigmak2);
k2=random(pd9,1000000,1);
pd10=makedist('Lognormal',muk3,sigmak3);
k3=random(pd10,1000000,1);

% Calculate the values for k1 and k3 through KD1 & k2 and KD2 & k4 respectively:
k1=k2./KD1;
k4=k3.*KD2;

% Calculate the new Keq through the mass action kinetics' parameters: 
Keqnew=(k1.*k3)./(k2.*k4);

% Calculate the new Michaelis-Menten parameters:
V1new=k3.*Eo;
V2new=k2.*Eo;
KAnew=(k2+k3)./k1;
KPnew=(k2+k3)./k4;

% Generate multivariate distribution by using the covariance matrix defined by the new distributions:
B=[Keqnew KAnew KPnew V1new V2new];
CorrMat=corrcoef(B);
MuB    = [muKeq mukmf mukmr muvmaxf muvmaxr];
SigmaB = [sigmaKeq sigmakmf sigmakmr sigmavmaxf sigmavmaxr];
% Force column vectors
MuB     = MuB(:);
SigmaB  = SigmaB(:);
% Calculate the covariance structure
sigma_down = repmat( SigmaB' , numel(SigmaB), 1            );
sigma_acrs = repmat( SigmaB  , 1           , numel(SigmaB) );
covv = log( CorrMat .* sqrt(exp(sigma_down.^2)-1) .* ...
                       sqrt(exp(sigma_acrs.^2)-1) + 1 );
% The Simulation
x = exp( mvnrnd( MuB , covv , ParamNo ));
Keqfin=x(:,1);
KAfin=x(:,2);
KPfin=x(:,3);
V1fin=x(:,4);
V2fin=x(:,5);

end