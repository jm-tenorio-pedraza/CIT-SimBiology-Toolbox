%% Search paths
warning off
clear all
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/UTG/PI')
%% Defining groups and data
groups_subset = {'MOC1_Control', 'MOC2_Control'};
observables={'TV'};
stateVar={'Tumor'};
PI=getPIData('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat',...
    stateVar,groups_subset,observables,'zeroAction', 'omit','mergePhenotypes', false,...
    'output', 'individual');

censored = repelem({'right'},length(PI.data),1);
[PI.data(1:end).censored] = censored{:,:};
%% Parameters
Tmax = 3100;
MOC1_T_0 = 0.0076;
MOC2_T_0 = 1.5e-4;
kel = 0.01;

kpro = [repelem({0.5736},12,1);repelem({1.0397},7,1)];
E = repelem({0.8},length(PI.data),1);
T_0 = [repelem({MOC1_T_0},12,1);repelem({MOC2_T_0},7,1)];

kill = repelem({0.1},19,1);

[PI.data(1:end).E] = E{:,:};
[PI.data(1:end).kpro] = kpro{:,:};
[PI.data(1:end).T_0] = T_0{:,:};
[PI.data(1:end).kill] = kill{:,:};
%% Model
model = @(kill, kpro, E, T_0,t)(Tmax./(kpro*(1/(kpro-kel-kill*E)+((Tmax-T_0)/T_0)*exp(-(kpro-kel-kill*E)*t))));

%% Non-linear LSQ
residuals_fn  = @(p)getResidualsKill(exp(p),PI,model,'free','kill'); 

lb = log([repelem(0.001,length(PI.data),1)]);
ub = log([repelem(10,length(PI.data),1)]);
options_fminsearch=optimset('Display','iter','MaxFunEvals', 1e4, 'MaxIter',1e4, 'TolFun', 1e-4);
p0=log([PI.data(:).kill]);
[phat, resnorm, res] = lsqnonlin(residuals_fn,p0, lb,ub, options_fminsearch);
%% Checking results
phat = num2cell(exp(phat_anneal));

%[PI.data(1:end).E] = phat{3:end};
kill=phat(1:length(PI.data));
[PI.data(1:end).kill] = kill{:};
simValues = arrayfun(@(x)model(x.kill,x.kpro,x.E,x.T_0,x.dataTime),PI.data,'UniformOutput',false);

[PI.data(1:end).simValues] = simValues{:,:};
figure
hold on
arrayfun(@(x) plot(x.dataTime, x.dataValue,'+', 'Color', x.colors),PI.data)
arrayfun(@(x) plot(x.dataTime, x.simValues,'-', 'Color', x.colors),PI.data)

%% MLE
options_anneal.Verbosity=2;
options_anneal.InitTemp=100;

likelihood_fun = @(x)getLikelihoodKill(exp(x),PI,model,'free','kill');
invgamma_prior= @(x,a,b)sum(log(b.^a./gamma(a).*x.^(-a-1).*exp(-b./x)));
unif_prior=@(x,indx)log((prod(and(x(indx)>=exp(lb(indx)),x(indx)<=exp(ub(indx))))));
lognorm_prior=@(x,m,s)sum(log(exp(-0.5 * ((log(x) - m)./s).^2) ./ (x .* sqrt(2*pi) .* s)));

prior_fun = @(x)unif_prior(exp(x'),1:19)+invgamma_prior(exp(x(end))^2,0.01,0.01);
p0 = [log([phat{1:19}]) log(0.1)];

obj_fun = @(p)(likelihood_fun((p))*(-1)+prior_fun((p))*(-1));
[phat_anneal, fval_anneal]=anneal(obj_fun,p0,options_anneal);

%% MCMC
finalValues = phat_anneal;
N = length(finalValues);

X0 =[ (finalValues); randn(100,length(finalValues))*0.005 + finalValues];
logL=rowfun(obj_fun,table(X0));
[L,I]=sort(logL{:,:});

w0=X0(I(1:N),:);

h.IndividualParams=[];
tic
[x, p_x,accept,pCR] = dreamH(w0,likelihood_fun,prior_fun,...
    size(w0,1),ceil(1e6/size(w0,1)), length(finalValues), 'BurnIn', 2e5,'StepSize',...
    1.38,'H', h);
toc
names = [repelem({'MOC1_{kill}'},12,1); repelem({'MOC2_{kill}'},7,1); '\sigma'];
plotMCMCDiagnostics(x, p_x,'name',names,'model', 'unperturbed tumor growth model')
postSamples =x(:,:,2e4:300:end);
logP_thinned = p_x(2e4:300:end,:);
plotMCMCDiagnostics(postSamples,logP_thinned,'name',...
    names,'model', 'unperturbed tumor growth model');

% Population Parameters
postSamples = postSamples(:,:)';
plotBivariateMarginals_2((postSamples),...
    'names',names)
kill = mat2cell(exp(postSamples(:,1:end-1))',repelem(1,length(PI.data),1));
[PI.data(1:end).kill]= kill{:,:};
% Killing rates in MOC1 tumors
figure('Renderer', 'painters', 'Position', [10 10 1500 600])

subplot(1,2,1)
hold on
arrayfun(@(x)histogram(x.kill,'Normalization','probability','FaceColor', x.colors),PI.data)
legend({PI.data(1:end).Name},'interpreter', 'none')
title('Killing rate distribution in MOC1 and MOC2 SyMMs')
xlabel('kill [10^{-6}cells \bullet day^{-1}]')
ylabel('Probability')


%% Predictions

simMat = arrayfun(@(x)nan(size(postSamples,1),length(x.dataTime)),PI.data,'UniformOutput',false);
[PI.data(1:end).postSim] = simMat{:,:};
[PI.data(1:end).postPred] = simMat{:,:};

for i=1:length(PI.data)

    simValues = rowfun(@(x)model(x,PI.data(i).kpro,PI.data(i).E,PI.data(i).T_0,...
        PI.data(i).dataTime)',table(exp(postSamples(:,i))));
    predValues = rowfun(@(x)exp(log(x)+randn(1,length(x))*exp(postSamples(i,end))), simValues);
    
    PI.data(i).postSim = simValues{:,:};
    PI.data(i).postPred = predValues{:,:};
        
end


postMean = arrayfun(@(x)mean(x.postSim),PI.data,'UniformOutput',false);
postMeanLB = arrayfun(@(x)quantile(x.postSim, 0.025),PI.data,'UniformOutput',false);
postMeanUB = arrayfun(@(x)quantile(x.postSim, 0.975),PI.data,'UniformOutput',false);
postPredLB = arrayfun(@(x)quantile(x.postPred, 0.025),PI.data,'UniformOutput',false);
postPredUB = arrayfun(@(x)quantile(x.postPred, 0.975),PI.data,'UniformOutput',false);

[PI.data(1:end).postMean] = postMean{:,:};
[PI.data(1:end).postMeanLB] = postMeanLB{:,:};
[PI.data(1:end).postMeanUB] = postMeanUB{:,:};
[PI.data(1:end).postPredLB] = postPredLB{:,:};
[PI.data(1:end).postPredUB] = postPredUB{:,:};

%% Plotting Posterior predictions results
% All plots in one figure
subplot(1,2,2)
hold on
arrayfun(@(x) plot(x.dataTime, x.dataValue,'--d', 'Color', x.colors,'MarkerSize',5,'MarkerEdgeColor', 'k'),PI.data)
arrayfun(@(x) plot(x.dataTime, x.postMean,'-', 'LineWidth', 2,'Color', x.colors),PI.data)
arrayfun(@(x) patch('XData', [x.dataTime', x.dataTime(end:-1:1)'],...
    'YData', [x.postPredUB x.postPredLB(end:-1:1)], 'FaceColor', x.colors,...
    'FaceAlpha', 0.2,'EdgeColor', 'none'),PI.data)
arrayfun(@(x) patch('XData', [x.dataTime', x.dataTime(end:-1:1)'],...
    'YData', [x.postMeanUB x.postMeanLB(end:-1:1)], 'FaceColor', x.colors,...
    'FaceAlpha', 0.3,'EdgeColor', 'none'),PI.data)

title('Posterior predictions for unperturbed tumor growth in SyMM')
ylabel('Volume [mL]')
xlabel('Time [days]')
legend('off')
% set(gca,'YScale', 'log')

% Single plot for each dataentry
ncol = ceil(sqrt(length(PI.data)));
nrow = ceil(length(PI.data)/ncol);
figure('Renderer', 'painters', 'Position', [10 10 1500 600])
for i=1:length(PI.data)
    subplot(nrow,ncol,i)
    hold on
    plot(PI.data(i).dataTime, PI.data(i).dataValue,'--d', 'Color', PI.data(i).colors,...
        'MarkerSize',5,'MarkerEdgeColor', 'k')
    plot(PI.data(i).dataTime, PI.data(i).postMean,'-', 'LineWidth', 2,'Color', PI.data(i).colors)
    patch('XData', [PI.data(i).dataTime', PI.data(i).dataTime(end:-1:1)'],...
    'YData', [PI.data(i).postPredUB PI.data(i).postPredLB(end:-1:1)], 'FaceColor', PI.data(i).colors,...
    'FaceAlpha', 0.2,'EdgeColor', 'none')
    patch('XData', [PI.data(i).dataTime', PI.data(i).dataTime(end:-1:1)'],...
        'YData', [PI.data(i).postMeanUB PI.data(i).postMeanLB(end:-1:1)],...
        'FaceColor', PI.data(i).colors,'FaceAlpha', 0.3,'EdgeColor', 'none')
    title(PI.data(i).Name, 'interpreter', 'none')
    ylabel('Volume [mL]')
xlabel('Time [days]')

end
%% CI estimates

[PI.par(1:size(postSamples,2)).name] = names{:,:};
PI = mcmcCI(PI,exp(postSamples),reshape(logP_thinned,[],1),0.95);
plotCI(PI,'unperturbed tumor growth in SyMM')
