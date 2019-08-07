function [Mode, Spread] = CalcModeSpread(V)

% This function generates the mode and confidence interval factor for a lognormal
% distribution, given a matrix P which contains 4 columns. 1)The parameter values 
% retrieved from literature 2)the error of measurements in parameter values, 
% 3)the weights of the parameter values and 4) a definition of the error as 
% multiplicative or additive. If the error is additive (Value +/- Error) input 0 
% and if it is multiplicative (Value */÷ Error) input 1 in the last column.
% Sample format of the table:
% V=[10000 NaN 2 0;100 0 8 0; 0.01 5 10 0; 160 NaN 4 0; 1 0.1 4 1];
% The input is given in the format: 
% [Mode, CIfact] = MCIcalc(V)--> input table V.

if nargin > 1
  error('Too many inputs'); 
end

D = [ ];
W = [ ];
    
for i=1:length(V(:,1))
    if V(i,4)==0 %logtransform additive SD
       lnE(i)=sqrt(log(1+(V(i,2).^2./V(i,1).^2)));
       if isnan(V(i,2)) %if SD is NaN use 10% multiplicative SD
       lnP(i)=log(V(i,1))-1/2.*log(1.1)^2;
       else
       lnP(i)=log(V(i,1))-1/2.*lnE(i).^2; 
       end
    else %logtransform multiplicative SD
        lnP(i)=log(V(i,1)); 
        lnE(i)=log(V(i,2));
    end   
end

V(:,1)= lnP;
V(:,2)= lnE;

% Sort table from smallest to largest parameter value
A=sortrows(V,1);

% Split table columns into separate vectors
P=A(:,1);
E=A(:,2);
Wo=A(:,3);

if any(Wo < 0.0001)
    error('The weights cannot have values smaller than 0.0001.');
end

for i=1:length(P)
    if isnan(E(i)) %if SD is NaN assign default 10% SD
       mu=P(i);  
       sigma=log(1.1);
       nbins = 1000; %generate bins within the parameter distribution range
       binEdges = linspace(mu-5*sigma,mu+5*sigma,nbins+1);
       aj = binEdges(1:end-1);
       aj=aj(:);
       
       bj = binEdges(2:end);
       bj=bj(:);
       cj = ( aj + bj ) ./ 2; %find centre of bins
       Pj = exp(-(cj-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
       Wj=Wo(i)*Pj.*(bj - aj); %calculate weight at the centre of bins
               
    elseif E(i)~=0 %if SD is not NaN or 0 use logtransformed SD
           mu=P(i);
           sigma=E(i);
           nbins = 1000;
           binEdges = linspace(mu-5*sigma,mu+5*sigma,nbins+1);
           aj = binEdges(1:end-1);
           aj=aj(:);
           bj = binEdges(2:end);
           bj=bj(:);
           cj = ( aj + bj ) ./ 2;
           Pj = exp(-(cj-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
           Wj=Wo(i)*Pj.*(bj - aj);          
    else 
        cj=P(i); %if SD is 0 do not assign SD but keep only the single value
        Wj=Wo(i);
            
    end
    
%if the value is not the minimum, it's not a single value and it does not overlap
%with the previous value, generate additional bins between twice the distance of P(i) and 
%P(i-1), otherwise do nothing   
     if P(i)~=min(P) && min(cj)>P(i-1) && length(cj)~=1 
          nbins2 = 1000; 
          binEdges2 = linspace(min(cj)-2*abs(min(cj)-P(i-1)),min(cj),nbins2+1);
          ajad = binEdges2(1:end-1);
          ajad=ajad(:);
          bjad = binEdges2(2:end);
          bjad=bjad(:);
          cjad = ( ajad + bjad ) ./ 2;
          Pjad = exp(-(cjad-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
          Wjad=Wo(i)*Pjad.*(bjad - ajad);
          
    else
       cjad =[]; 
       Wjad =[];
       
    end

%if the value is not the maximum, it's not a single value and it does not overlap
%with the next value, generate additional bins between twice the distance of P(i) and 
%P(i+1), otherwise do nothing       
       if P(i)~=max(P)&& max(cj)<P(i+1) && length(cj)~=1
          nbins3 = 1000;
          binEdges3 = linspace(max(cj),max(cj)+2*abs(P(i+1)-max(cj)),nbins3+1);
          ajad2 = binEdges3(1:end-1);
          ajad2=ajad2(:);
          bjad2 = binEdges3(2:end);
          bjad2=bjad2(:);
          cjad2 = ( ajad2 + bjad2 ) ./ 2;
          Pjad2 = exp(-(cjad2-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
          Wjad2=Wo(i)*Pjad2.*(bjad2 - ajad2);
       else 
           cjad2 =[]; 
           Wjad2 =[];
           
       end
    
    D = [D ; cj ; cjad ; cjad2]; %add the centres of all bins in a vector
    W = [W ; Wj ; Wjad ; Wjad2]; %add the weights of all bin centres in a vector
    
end

wMed = weightedMedian(D,W); %calculate the weighted median of values in vector D

S = std(D,W); %calculate the weighted standard deviation of values in vector D

Mode=exp(wMed); %calculate Mode
Spread=exp(S); %calculate Spread
