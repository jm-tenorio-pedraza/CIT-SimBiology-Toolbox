function [Keqfin, KAfin, KBfin, KPfin, KQfin, V1fin, V2fin] = MultivariateBiBi(muKeq, sigmaKeq, mukA, sigmakA, mukB, sigmakB, mukP, sigmakP, mukQ, sigmakQ, muvmaxf, sigmavmaxf, muvmaxr, sigmavmaxr, ParamNo)
% The function generates a multivariate distribution for Bi-Bi reactions with   
% Michaelis-Menten kinetics. The user must provide the  
% mu and sigma for the 7 distributions in the correct order,as well as the  
% number of parameter samples that are needed (ParamNo).

%Define lognormal distributions:
% Keq 
pd1=makedist('Lognormal',muKeq,sigmaKeq);
Keq=random(pd1,1000000,1);

% kA 
pd2=makedist('Lognormal',mukA, sigmakA);
kA=random(pd2,1000000,1);

% kB 
pd3=makedist('Lognormal',mukB, sigmakB);
kB=random(pd3,1000000,1);

% kP 
pd4=makedist('Lognormal',mukP, sigmakP);
kP=random(pd4,1000000,1);

% kQ 
pd5=makedist('Lognormal',mukQ, sigmakQ);
kQ=random(pd5,1000000,1);

% V1 
pd6=makedist('Lognormal',muvmaxf, sigmavmaxf);
V1=random(pd6,1000000,1);

% V2 
pd7=makedist('Lognormal',muvmaxr, sigmavmaxr);
V2=random(pd7,1000000,1);

% Calculate geometric coefficients of variation:
GCV1=exp(sigmaKeq)-1;
GCV2=exp(sigmakA)-1;
GCV3=exp(sigmakB)-1;
GCV4=exp(sigmakP)-1;
GCV5=exp(sigmakQ)-1;
GCV6=exp(sigmavmaxf)-1;
GCV7=exp(sigmavmaxr)-1;

% Choose the parameter with the largest coefficient of variation:
A=[GCV1 GCV2 GCV3 GCV4 GCV5 GCV6 GCV7];
M=max(A);

% Set dependent parameter based on the largest GCV and 
% calculate new mu and sigma for the dependent distribution:

if M==GCV2  % kA dependent
   kA=(kP.*kQ.*(V1.^2))./(Keq.*kB.*(V2.^2));
   mukA=(mukP+mukQ+2*muvmaxf)-(muKeq+mukB+2*muvmaxr);
   sigmakA=sqrt(sigmaKeq^2+sigmakB^2+sigmakP^2+sigmakQ^2+(2*sigmavmaxr)^2+(2*sigmavmaxf)^2);
   
elseif M==GCV1 % Keq dependent
   Keq=(kP.*kQ.*(V1.^2))./(kA.*kB.*(V2.^2));
   muKeq=(2*muvmaxf+mukP+mukQ)-(mukA+mukB+2*muvmaxr);
   sigmaKeq=sqrt(sigmakA^2+sigmakB^2+sigmakP^2+sigmakQ^2+(2*sigmavmaxr)^2+(2*sigmavmaxf)^2);
   
elseif M==GCV3 % kB dependent
   kB=(kP.*kQ.*(V1.^2))./(kA.*Keq.*(V2.^2));
   mukB=(mukP+mukQ+2*muvmaxf)-(muKeq+mukA+2*muvmaxr);
   sigmakB=sqrt(sigmakA^2+sigmaKeq^2+sigmakP^2+sigmakQ^2+(2*sigmavmaxr)^2+(2*sigmavmaxf)^2);
   
elseif M==GCV4 % kP dependent
   kP=(kA.*kB.*Keq.*(V2.^2))./(kQ.*(V1.^2));
   mukP=(muKeq+mukA+mukB+2*muvmaxr)-(mukQ+2*muvmaxf);
   sigmakP=sqrt(sigmakA^2+sigmaKeq^2+sigmakB^2+sigmakQ^2+(2*sigmavmaxr)^2+(2*sigmavmaxf)^2);
   
elseif M==GCV5 % kQ dependent
   kQ=(kA.*kB.*Keq.*(V2.^2))./(kP.*(V1.^2));
   mukQ=(muKeq+mukA+mukB+2*muvmaxr)-(mukP+2*muvmaxf);
   sigmakQ=sqrt(sigmakA^2+sigmaKeq^2+sigmakB^2+sigmakP^2+(2*sigmavmaxr)^2+(2*sigmavmaxf)^2);

elseif M==GCV6 % vmaxf dependent
   V1=sqrt((kA.*kB.*Keq.*(V2.^2))./(kQ.*kP));
   muvmaxf=(0.5*muKeq+0.5*mukA+0.5*mukB+muvmaxr)-(0.5*mukP+0.5*mukQ);
   sigmavmaxf=sqrt(sigmakA^2+sigmaKeq^2+sigmakB^2+sigmakP^2+sigmakQ^2+(2*sigmavmaxr)^2);

elseif M==GCV6 % vmaxr dependent
   V2=sqrt((kP.*kQ.*(V1.^2))./(kA.*Keq.*kB));
   muvmaxr=(0.5*mukP+0.5*mukQ+muvmaxf)-(0.5*muKeq+0.5*mukA+0.5*mukB);
   sigmavmaxr=sqrt(sigmakA^2+sigmaKeq^2+sigmakB^2+sigmakP^2+sigmakQ^2+(2*sigmavmaxf)^2);
end

if M==GCV1 % Keq dependent
   % Generate multivariate distribution:
   B=[Keq kA kB kP kQ V1];
   CorrMat=corrcoef(B);
   Mu    = [muKeq mukA mukB mukP mukQ muvmaxf];
   Sigma = [sigmaKeq sigmakA sigmakB sigmakP sigmakQ sigmavmaxf];
   % Force column vectors
   Mu     = Mu(:);
   Sigma  = Sigma(:);
   % Calculate the covariance structure
   sigma_down = repmat( Sigma' , numel(Sigma), 1            );
   sigma_acrs = repmat( Sigma  , 1           , numel(Sigma) );
   covv = log( CorrMat .* sqrt(exp(sigma_down.^2)-1) .* ...
                       sqrt(exp(sigma_acrs.^2)-1) + 1 );
   % The Simulation
   y = exp( mvnrnd( Mu , covv , ParamNo ));
   Keqfin=y(:,1);
   KAfin=y(:,2);
   KBfin=y(:,3);
   KPfin=y(:,4);
   KQfin=y(:,5);
   V1fin=y(:,6);
   V2fin=sqrt((KPfin.*KQfin.*(V1fin.^2))./(KAfin.*Keqfin.*KBfin));
   
   else
   % Generate multivariate distribution:
   B=[kA kB kP kQ V1 V2];
   CorrMat=corrcoef(B);
   Mu    = [mukA mukB mukP mukQ muvmaxf muvmaxr];
   Sigma = [sigmakA sigmakB sigmakP sigmakQ sigmavmaxf sigmavmaxr];
   % Force column vectors
   Mu     = Mu(:);
   Sigma  = Sigma(:);
   % Calculate the covariance structure
   sigma_down = repmat( Sigma' , numel(Sigma), 1            );
   sigma_acrs = repmat( Sigma  , 1           , numel(Sigma) );
   covv = log( CorrMat .* sqrt(exp(sigma_down.^2)-1) .* ...
                       sqrt(exp(sigma_acrs.^2)-1) + 1 );
   % The Simulation
   y = exp( mvnrnd( Mu , covv , ParamNo ));
   KAfin=y(:,1);
   KBfin=y(:,2);
   KPfin=y(:,3);
   KQfin=y(:,4);
   V1fin=y(:,5);
   V2fin=y(:,6);
   Keqfin=(KPfin.*KQfin.*(V1fin.^2))./(KAfin.*KBfin.*(V2fin.^2));
   
end 