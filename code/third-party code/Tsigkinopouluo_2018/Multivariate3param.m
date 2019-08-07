function [KD1, kon1, koff1] = Multivariate3param(muKD, sigmaKD, mukon, sigmakon, mukoff, sigmakoff, ParamNo)
% The function generates a multivariate distribution for mass action reactions 
% with 3 parameters. The user must provide the mu and sigma for the three 
% distributions in the correct order, as well as the number of parameter 
% samples that are needed (ParamNo).

%Define lognormal distributions:
% KD 
pd1=makedist('Lognormal',muKD,sigmaKD);
KD=random(pd1,1000000,1);

% kon 
pd2=makedist('Lognormal',mukon,sigmakon);
kon=random(pd2,1000000,1);

% koff 
pd3=makedist('Lognormal',mukoff,sigmakoff);
koff=random(pd3,1000000,1);

% Calculate coefficients of variation:

GCV1=exp(sigmaKD)-1;
GCV2=exp(sigmakon)-1;
GCV3=exp(sigmakoff)-1;

% Choose the parameter with the largest coefficient of variation:
A=[GCV1 GCV2 GCV3];
M=max(A);

% Set dependent parameter based on the largest GCV and 
% calculate new mu and sigma for the dependent distribution:

if M==GCV1 %KD dependent
   KD=koff./kon;
   muKD=mukoff-mukon;
   sigmaKD=sqrt(sigmakoff^2+sigmakon^2);
     
elseif M==GCV2 % kon dependent
   kon=koff./KD; 
   mukon=mukoff-muKD;
   sigmakon=sqrt(sigmakoff^2+sigmaKD^2);

elseif M==GCV3 % koff dependent
   koff=kon.*KD;
   mukoff=mukon+muKD;
   sigmakoff=sqrt(sigmakon^2+sigmaKD^2);
end
   
if M==GCV2 % kon dependent
   % Generate multivariate distribution:
   B=[kon koff];
   CorrMat=corrcoef(B);
   Mu    = [mukon mukoff];
   Sigma = [sigmakon sigmakoff];
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
   kon1=y(:,1);
   koff1=y(:,2);
   KD1=koff1./kon1;
   
else % koff or KD dependent
   % Generate multivariate distribution:
   B=[KD koff];
   CorrMat=corrcoef(B);
   Mu    = [muKD mukoff];
   Sigma = [sigmaKD sigmakoff];
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
   KD1=y(:,1);
   koff1=y(:,2);
   kon1=koff1./KD1;

end

