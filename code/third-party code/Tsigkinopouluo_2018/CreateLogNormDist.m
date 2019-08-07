function [mu, sigma, Xmin, Xmax] = CreateLogNormDist(Mode, Spread, Percentage)

% This function generates the location and scale parameters for a lognormal
% distribution, given a mode, a Spread and a percentage of values that will 
% be contained within the defined boundaries. The input can be given in two ways: 
% [mu, sigma,Xmin, Xmax] = CreateLogNormDist(Mode, Spread)--> only given Mode and  
% Spread. In this case, the percentage is set to default value (0.6827).
% [mu, sigma, Xmin, Xmax] = CreateLogNormDist(Mode, Spread, Percentage)--> given 
% Mode, Spread and percentage. The function then generates the location and 
% scale parameters mu and sigma of the lognormal distribution, as well as  
% the Xmin(=Mode/Spread) and Xmax(=Mode*Spread).

% Error checking
if nargin < 2
    error('Must have at least 2 input arguments');
end

if nargin > 3
  error('Too many inputs'); 
end

% Assign default percentage if not given
if nargin == 2
  Percentage=0.6827;
end

% Calculate Xmin and Xmax of the range of values from the Mode and Spread

Xmin=Mode/Spread;
Xmax=Mode*Spread;

% Calculate mu and sigma of the lognormal distribution
syms s
eqn = 1/2+1/2*erf((log(Xmax)-(log(Mode)+s^2))/(sqrt(2)*s))-(1/2+1/2*erf((log(Xmin)-(log(Mode)+s^2))/(sqrt(2)*s)))==Percentage;
sigma = vpasolve(eqn,s);
mu=log(Mode)+sigma^2;

%Plot final distribution
mu=double(mu);
sigma=double(sigma);
pd=makedist('Lognormal',mu,sigma);
Dist=random(pd,10000,1);
figure()
histfit(Dist,80,'lognormal'); 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','c','facealpha',0.3);
htitle = title('Distribution for generated parameter values'); hxlabel=xlabel('Parameter values'); hylabel=ylabel('Kernel Density'); 
set(gca,'FontSize',18) % font size of axis values
set(htitle,'FontSize',26) % font size of title
set(hxlabel,'FontSize',20) % font size of x axis
set(hylabel,'FontSize',20) % font size of y axis
end