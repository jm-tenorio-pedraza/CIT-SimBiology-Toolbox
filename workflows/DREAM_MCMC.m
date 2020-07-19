%% DREAM MCMC
finalValues = log([PI.par(:).finalValue]);
N = length(finalValues);

X0 =[ (finalValues); randn(100,length(finalValues))*0.1 + finalValues];
logL=rowfun(obj_fun,table(X0));
[L,I]=sort(logL{:,:});
w0 = X0(I(1:N),:);
h.IndividualParams=[];
p=parpool('local');
tic
[x1, p_x1,accept1,pCR1,stepSize1, J1, n_id1] = par_dream(w0',likelihood_fun,prior_fun,...
    N,ceil(1e6/N), N, 'BurnIn', ...
    2e5,'stepSize',2.38);
toc

w0 = x1(:,:,end);
tic
[x2, p_x2,accept2,pCR2,stepSize2,J2,n_id2] = par_dream(w0,likelihood_fun,prior_fun,...
    N,ceil(1e6/N), N, 'BurnIn', ...
    0,'StepSize',stepSize1,'pCR', pCR1, 'J', J1, 'n_id', n_id1);
toc
w0 = x2(:,:,end);

tic
[x3, p_x3,accept3,pCR3,stepSize3,J3,n_id3] = par_dream(w0,likelihood_fun,prior_fun,...
    N,ceil(1e6/size(w0,1)), N, 'BurnIn', ...
    0,'StepSize',stepSize2,'pCR', pCR2,'J', J2, 'n_id', n_id2);
toc

w0 = x3(:,:,end);
tic
[x4, p_x4,accept4,pCR4,stepSize4,J4, n_id4] = dreamHParallel(w0',likelihood_fun,...
    prior_fun,size(w0,1),ceil(1e6/size(w0,1)), length(finalValues),...
    'BurnIn', 5e5,'StepSize',stepSize3,'H', h,'pCR', pCR3, 'J', J3, 'n_id', n_id3);
toc

w0 = x4(:,:,end);
tic
[x5, p_x5,accept5,pCR5,stepSize5,J5, n_id5] = dreamHParallel(w0',likelihood_fun,prior_fun_MCMC,...
    size(w0,1),ceil(3e5/size(w0,1)), length(finalValues), 'BurnIn', ...
    3e5,'StepSize',stepSize4,'H', h, 'pCR', pCR4, 'J', J4, 'n_id', n_id4);
toc

w0 = x5(:,:,end);
tic
[x6, p_x6,accept6,pCR6,stepSize6,J6, n_id6] = dreamHParallel(w0',likelihood_fun,prior_fun_MCMC,...
    size(w0,1),ceil(3e5/size(w0,1)), length(finalValues), 'BurnIn', ...
    3e5,'StepSize',stepSize5,'H', h, 'pCR', pCR5, 'J', J5, 'n_id', n_id5);
toc

w0 = x6(:,:,end);
tic
[x7, p_x7,accept7,pCR7,stepSize7,J7, n_id7] = dreamHParallel(w0',likelihood_fun,prior_fun_MCMC,...
    size(w0,1),ceil(3e5/size(w0,1)), length(finalValues), 'BurnIn', ...
    3e5,'StepSize',.4,'H', h, 'pCR', pCR6, 'J', J6, 'n_id', n_id6);
toc

w0 = x7(:,:,end);
tic
[x8, p_x8,accept8,pCR8, stepSize8, J8, n_id8] = dreamHParallel(w0',likelihood_fun,prior_fun_MCMC,...
    size(w0,1),ceil(3e5/size(w0,1)), length(finalValues), 'BurnIn', ...
    3e5,'StepSize',stepSize7,'H', h, 'pCR', pCR7, 'J', J7, 'n_id', n_id7);
toc

w0 = x8(:,:,end);
tic
[x9, p_x9,accept9,pCR9, stepSize9, J9, n_id9] = dreamHParallel(w0',likelihood_fun,prior_fun_MCMC,...
    size(w0,1),ceil(3e5/size(w0,1)), length(finalValues), 'BurnIn', ...
    3e5,'StepSize',stepSize8,'H', h, 'pCR', pCR8, 'J', J8, 'n_id', n_id8);
toc
x1=cat(3,x1,x2,x3,x4,x5,x6, x7, x8, x9);
p_x1 = [p_x1; p_x2; p_x3; p_x4; p_x5; p_x6; p_x7; p_x8; p_x9];

%% 

covM = cov(postSamples);
R = chol(covM);
X0 = [finalValues; (randn(100, length(finalValues))*0.1)*R + finalValues];
logL=rowfun(obj_fun,table(X0));
[L,I]=sort(logL{:,:});
w0=X0(I(1:N),:);

tic
[x10, p_x10,accept10,pCR10, stepSize10, J10, n_id10] = dreamHParallel(w0,likelihood_fun,prior_fun_MCMC,...
    size(w0,1),ceil(3e5/size(w0,1)), length(finalValues), 'BurnIn', ...
    3e5,'StepSize',2.38^2,'H', h);
toc

w0=x10(:,:,end)';

tic
[x11, p_x11,accept11,pCR11, stepSize11, J11, n_id11] = dreamHParallel(w0,likelihood_fun,prior_fun_MCMC,...
    size(w0,1),ceil(3e5/size(w0,1)), length(finalValues), 'BurnIn', ...
    3e5,'StepSize',stepSize10,'H', h, 'pCR', pCR10, 'J', J10, 'n_id', n_id10);
toc


w0=x11(:,:,end)';

tic
[x12, p_x12,accept12,pCR12, stepSize12, J12, n_id12] = dreamHParallel(w0,likelihood_fun,prior_fun_MCMC,...
    size(w0,1),ceil(3e5/size(w0,1)), length(finalValues), 'BurnIn', ...
    0,'StepSize',stepSize10,'H', h, 'pCR', pCR11, 'J', J11, 'n_id', n_id11);
toc

w0=x12(:,:,end)';

tic
[x13, p_x13,accept13,pCR13, stepSize13, J13, n_id13] = dreamHParallel(w0,likelihood_fun,prior_fun_MCMC,...
    size(w0,1),ceil(3e5/size(w0,1)), length(finalValues), 'BurnIn', ...
    0,'StepSize',stepSize12,'H', h, 'pCR', pCR12, 'J', J12, 'n_id', n_id12);
toc


x1 = cat(3,x10, x11, x12, x13);
p_x1 = [p_x10; p_x11; p_x12;p_x13];

%%
covM = cov(postSamples);
R = chol(covM);
X0 = [finalValues; (randn(100, length(finalValues))*0.1)*R + finalValues];
logL=rowfun(obj_fun,table(X0));
[L,I]=sort(logL{:,:});
w0=X0(I(1:N),:);

tic
[x14, p_x14,accept14,pCR14, stepSize14, J14, n_id14] = dreamHParallel(w0,likelihood_fun,prior_fun_MCMC,...
    size(w0,1),ceil(2.5e5/size(w0,1)), length(finalValues), 'BurnIn', ...
    3e5,'StepSize',2.38^2,'H', h);
toc

w0=x14(:,:,end)';

tic
[x15, p_x15,accept15, pCR15, stepSize15, J15, n_id15] = dreamHParallel(w0,likelihood_fun,prior_fun_MCMC,...
    size(w0,1),ceil(3e5/size(w0,1)), length(finalValues), 'BurnIn', ...
    3e5,'StepSize',stepSize14,'H', h, 'pCR', pCR14, 'J', J14,'n_id', n_id14);
toc

w0=x15(:,:,end)';

tic
[x16, p_x16, accept16, pCR16, stepSize16, J16, n_id16] = dreamHParallel(w0,likelihood_fun,prior_fun_MCMC,...
    size(w0,1),ceil(3e5/size(w0,1)), length(finalValues), 'BurnIn', ...
    3e5,'StepSize',stepSize15,'H', h, 'pCR', pCR15, 'J', J15, 'n_id', n_id15);
toc


w0=x16(:,:,end)';

tic
[x17, p_x17, accept17, pCR17, stepSize17, J17, n_id17] = dreamHParallel(w0,likelihood_fun,prior_fun_MCMC,...
    size(w0,1),ceil(3e5/size(w0,1)), length(finalValues), 'BurnIn', ...
    3e5,'StepSize',stepSize16,'H', h, 'pCR', pCR16, 'J', J16, 'n_id', n_id16);
toc
w0=x17(:,:,end)';

tic
[x18, p_x18, accept18, pCR18, stepSize18, J18, n_id18] = dreamHParallel(w0,likelihood_fun,prior_fun_MCMC,...
    size(w0,1),ceil(3e5/size(w0,1)), length(finalValues), 'BurnIn', ...
    0,'StepSize',stepSize17,'H', h, 'pCR', pCR17, 'J', J17, 'n_id', n_id17);
toc
x1 = cat(3,x14, x15, x16,x17, x18);
p_x1 = [p_x14; p_x15; p_x16; p_x17; p_x18];

%% Save results
save(strjoin({cd '/PI_CIM21_Control_11_2_DREAM_MCMC_x.mat'},''), 'x')
save(strjoin({cd '/PI_CIM21_Control_11_2_DREAM_MCMC_p_x.mat'},''), 'p_x')

load(strjoin({cd '/CIM_red2_DREAM_MCMC_x.mat'},''))
load(strjoin({cd '/CIM_red2_DREAM_MCMC_p_x.mat'},''))
