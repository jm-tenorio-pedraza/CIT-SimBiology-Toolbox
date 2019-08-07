function [Keqfin, KAfin, KPfin, KiQfin, V1fin, V2fin] = MultivariateUniBi(muKeq, sigmaKeq, mukA, sigmakA, mukP, sigmakP, mukiQ, sigmakiQ, muvmaxf, sigmavmaxf, muvmaxr, sigmavmaxr, ParamNo)
% The function generates a multivariate distribution for Uni-Bi reactions with   
% Michaelis-Menten kinetics. The user must provide the  
% mu and sigma for the 6 distributions in the correct order,as well as the  
% number of parameter samples that are needed (ParamNo).

%Define lognormal distributions:
% Keq 
pd1=makedist('Lognormal',muKeq,sigmaKeq);
Keq=random(pd1,1000000,1);

% kA 
pd2=makedist('Lognormal',mukA, sigmakA);
kA=random(pd2,1000000,1);

% kP 
pd3=makedist('Lognormal',mukP, sigmakP);
kP=random(pd3,1000000,1);

% kiQ 
pd4=makedist('Lognormal',mukiQ, sigmakiQ);
kiQ=random(pd4,1000000,1);

% V1 
pd5=makedist('Lognormal',muvmaxf, sigmavmaxf);
V1=random(pd5,1000000,1);

% V2 
pd6=makedist('Lognormal',muvmaxr, sigmavmaxr);
V2=random(pd6,1000000,1);

% Calculate geometric coefficients of variation:
GCV1=exp(sigmaKeq)-1;
GCV2=exp(sigmakA)-1;
GCV3=exp(sigmakP)-1;
GCV4=exp(sigmakiQ)-1;
GCV5=exp(sigmavmaxf)-1;
GCV6=exp(sigmavmaxr)-1;

% Choose the parameter with the largest coefficient of variation:
A=[GCV1 GCV2 GCV3 GCV4 GCV5 GCV6];
M=max(A);

% Set dependent parameter based on the largest GCV and 
% calculate new mu and sigma for the dependent distribution:

if M==GCV2  % kA dependent
   kA=(V1.*kiQ.*kP)./(V2.*Keq);
   mukA=(mukP+mukiQ+muvmaxf)-(muKeq+muvmaxr);
   sigmakA=sqrt(sigmaKeq^2+sigmakP^2+sigmakiQ^2+sigmavmaxr^2+sigmavmaxf^2);
   
elseif M==GCV1 % Keq dependent
   Keq=(kP.*kiQ.*V1)./(kA.*V2);
   muKeq=(muvmaxf+mukP+mukiQ)-(mukA+muvmaxr);
   sigmaKeq=sqrt(sigmakA^2+sigmakP^2+sigmakiQ^2+sigmavmaxr^2+sigmavmaxf^2);
   
elseif M==GCV3 % kP dependent
   kP=(kA.*Keq.*V2)./(kiQ.*V1);
   mukP=(muKeq+mukA+muvmaxr)-(mukiQ+muvmaxf);
   sigmakP=sqrt(sigmaKeq^2+sigmakA^2+sigmakiQ^2+sigmavmaxr^2+sigmavmaxf^2);
   
elseif M==GCV4 % kiQ dependent
   kiQ=(kA.*Keq.*V2)./(kP.*V1);
   mukiQ=(muKeq+mukA+muvmaxr)-(mukP+muvmaxf);
   sigmakiQ=sqrt(sigmaKeq^2+sigmakA^2+sigmakP^2+sigmavmaxr^2+sigmavmaxf^2);
   
elseif M==GCV5 % vmaxf dependent
   V1=(kA.*Keq.*V2)./(kiQ.*kP);
   muvmaxf=(muKeq+mukA+muvmaxr)-(mukiQ+mukP);
   sigmavmaxf=sqrt(sigmaKeq^2+sigmakA^2+sigmakiQ^2+sigmavmaxr^2+sigmakP^2);
   
elseif M==GCV6 % vmaxr dependent
   V2=(V1.*kiQ.*kP)./(kA.*Keq);
   muvmaxr=(mukP+mukiQ+muvmaxf)-(muKeq+mukA);
   sigmavmaxr=sqrt(sigmaKeq^2+sigmakP^2+sigmakiQ^2+sigmakA^2+sigmavmaxf^2);
end

if M==GCV1 % Keq dependent
   % Generate multivariate distribution:
   B=[Keq kA kP kiQ V1];
   CorrMat=corrcoef(B);
   Mu    = [muKeq mukA mukP mukiQ muvmaxf];
   Sigma = [sigmaKeq sigmakA sigmakP sigmakiQ sigmavmaxf];
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
   KPfin=y(:,3);
   KiQfin=y(:,4);
   V1fin=y(:,5);
   V2fin=(V1fin.*KiQfin.*KPfin)./(KAfin.*Keqfin);

   else
   % Generate multivariate distribution:
   B=[kA kP kiQ V1 V2];
   CorrMat=corrcoef(B);
   Mu    = [mukA mukP mukiQ muvmaxf muvmaxr];
   Sigma = [sigmakA sigmakP sigmakiQ sigmavmaxf sigmavmaxr];
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
   KPfin=y(:,2);
   KiQfin=y(:,3);
   V1fin=y(:,4);
   V2fin=y(:,5);
   Keqfin=(KPfin.*KiQfin.*V1fin)./(KAfin.*V2fin);
   
end  
   
   