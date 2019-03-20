% Runs nested particle filter on supplied model and parameters
% Represent parameters and state variables as particles and update their
% distributions each iteration

% In this script we generate observations to perform TWIN EXPERIMENT to
% validate our method for time-varying parameter estimation
clearvars;
close all;
%% set up parameters, in this case we use functions to generate time varying parameters
ind = [2];
% functions to generate time dependent parameter
% === uncomment Fpara1 if want to estimate Vp as well
% Fpara1 = @(t) 10*sin(1/10^5/2*t) + 100; % Vp
Fpara2 = @(t) 1*cos(1/10^5*t) +  5; % S
pmin = [3]; %particles originally uniformly distributed b/w pmin and pmax
pmax = [9];

% stochastic bimodal model
model = 'FAST'; % which model to use, here we test on Kerry's 2017 FAST model
mdim=2; %state dimension

%% critical parameters
N = 1000; %number of particles in filter
M = 20; % number of state particle per parameter particle
o_sigma=0.4; %observational error std
resamp_thresh=1; %threshhold for resampling
wiggle = 0.05; %noise added on resampling for simple resampling scheme

sigma=0;%0.05; %standard deviation of model noise

alpha = 0.98; beta = 0.5; %beta = [0.01;0.01];% jittering parameter

%% numerical integration parameters
T=20*2*pi*10^5; %final time, three periods
dt=2.5*10^3;% time step for `truth' integration
obs_step=8*10^4; %time interval between observations
m = 1;
mdt = m*dt; %`model' time step
t_truth = 0:dt:T; % truth times
tvals = 0 : obs_step : T; %analysis times
tdim = length(tvals);
xt = [50; 0.8]; %truth - initial condition
xp = xt+ o_sigma.*randn(mdim,1); %initial condition for model state
xp = repmat(xp,[1,M,N]); %forecast associated with each particle - dimension mdim x N

%% model & model parameters

% parameters for FAST, param contains the parameters that are not estimated
Cd = 1.2*10^-3;
h = 1400; % in m
epsilon = 0.33;
kappa = 0.1;
hm = 50;
Vp = 100;
S = 5;
Gamma = 2.5;
uT = 5;

% parameters in the truth time resolution
theta_true = [Vp*ones(size(t_truth)); S*ones(size(t_truth)); Gamma*ones(size(t_truth)); hm*ones(size(t_truth)); ...
    uT*ones(size(t_truth)); Cd*ones(size(t_truth)); h*ones(size(t_truth)); epsilon*ones(size(t_truth)); kappa*ones(size(t_truth))];
theta_true(ind,:) = [Fpara2(t_truth)];
true_para = theta_true(ind,:);
% parameters in the analysis time resolution
theta = [Vp*ones(size(tvals)); S*ones(size(tvals)); Gamma*ones(size(tvals)); hm*ones(size(tvals)); ...
    uT*ones(size(tvals)); Cd*ones(size(tvals)); h*ones(size(tvals)); epsilon*ones(size(tvals)); kappa*ones(size(tvals))];
theta(ind,:) = [Fpara2(tvals)];
theta = repmat(theta, [1, 1, N]); % extend to N particles

%% truth and observations
truth=model;

obs_op=@(x)x; % observational operator
obsdim=2; %dimension of obervations
R=o_sigma^2*eye(obsdim); %observation covariance matrix

%% generate particle distribution
pdim = length(ind); %dimension of parameters
particle = samp_param(pmin,pmax,N);
Wpara = 1/N*ones(N,1); %initial weights
Wstate = 1/M*ones(M,1);

%% generate truth and observations to run twin experiments to validate the method

xt = EulerM_hist_theta_time(truth,tvals,dt,xt,sigma,theta_true); % truth
obshist = obs_op(xt) + repmat(o_sigma,[2,tdim]).*randn(mdim,tdim); % observations

phist = zeros(pdim,N,tdim); %to store parameters particles
pstd = zeros(pdim, tdim);
Whist = zeros(N,tdim); % to store parameters particles weights

phist(:,:,1) = particle; % initialize
Whist(:,1) = Wpara; % initialize
para_est = zeros(pdim, tdim); % estimation of parameters
para_est(:,1) = particle*Wpara; 
x_est = zeros(mdim, tdim); % estimation of state variables
x_est(:,1) = reshape(xp(:,:,1),[mdim,M])*Wstate;
xi = zeros(mdim, tdim, N); % mean in state variables amoung M state particles for each N parameter particles
resampcount=0; %count how many times the PF resampled

%% do particle filter
for tau=1:tdim-1
    obs = obshist(:,tau+1);
    % jittering, important step!! see Liu & West
    particle = alpha*particle + (1-alpha)*repmat(particle*Wpara, 1, N) + repmat(beta,[1,N]).*randn(size(particle));
    % update state particle for each theta_ii, ii = 1...N
    theta(ind, tau, :) = particle;
    for ii = 1:N
       % prediction
       xp(:,:,ii) = EulerM_theta(model,tvals(tau:tau+1),mdt,xp(:,:,ii),sigma,theta(:, tau, ii)); 
       
       xi(:,tau,ii) = 1/M*sum(xp(:,:,ii),2);
       % update weights of state particles
       innov_state = abs(repmat(obs,1,M) - obs_op(xp(:,:,ii)));
       Wstatetmp = -0.5*sum((innov_state.^2)./R(1),1);
       Wstatemax=max(Wstatetmp);
       Wstatetmp  = Wstatetmp-Wstatemax;
       Wstate = Wstate.*exp(Wstatetmp');
       Wstate=Wstate/sum(Wstate);
       
       % resample in x_i^j, for j = 1...M. NOTE no resampling condition

        a1 = 3/4; a2=(sqrt(13)+1)/8; a3 = -(sqrt(13)-1)/8;
        sampIndex1 = ResampSimp(Wstate,M);
        sampIndex2 = ResampSimp(Wstate,M);
        sampIndex3 = ResampSimp(Wstate,M);
        xp(:,:,ii) = a1*xp(:,sampIndex1,ii)+...
            a2*xp(:,sampIndex2,ii)+...
            a3*xp(:,sampIndex3,ii);
        Wstate = 1/M*ones(M,1);
       
    end
    
    % update parameter particles' weights
    innov = abs(repmat(obs,1,N) - obs_op(reshape(xi(:,tau,:), [mdim, N])));
    Wtmp = -0.5*sum(innov.^2./R(1),1);
    Wmax = max(Wtmp);
    Wtmp = Wtmp - Wmax;
    Wpara = Wpara.*exp(Wtmp');
    Wpara = Wpara/sum(Wpara);
    
    % resample if weights are concentrated
    if 1/sum(Wpara.^2)/N < resamp_thresh
        resampcount = resampcount+1;
        
        a1 = 3/4; a2=(sqrt(13)+1)/8; a3 = -(sqrt(13)-1)/8;
             sampIndex1 = ResampSimp(Wpara,N);
             sampIndex2 = ResampSimp(Wpara,N);
             sampIndex3 = ResampSimp(Wpara,N);
             particle = a1*particle(:,sampIndex1)+...
                        a2*particle(:,sampIndex2)+...
                        a3*particle(:,sampIndex3);
              xp = a1*xp(:,:,sampIndex1)+...
                   a2*xp(:,:,sampIndex2)+...
                   a3*xp(:,:,sampIndex3);

        Wpara = 1/N*ones(N,1);
    end
    
    % compute the time evolved estimation
    para_est(:, tau+1) = particle*Wpara;
    % compute the estimation in state variable
    x_temp = zeros(mdim, N); % temporary variable to store the mean of state variables amoung state particles
    for ii = 1:N
        x_temp(:, ii) = reshape(xp(:,:,ii),[mdim,M])*Wstate;
    end
    
    x_est(:, tau+1) = x_temp*Wpara;

    % compute standard deviation
    pstd(1,tau) = std(particle(1,:), Wpara); 
end

% plotting
figure;
plot(tvals,para_est(1,:),'color','b')
hold on;
plot(tvals,Fpara2(tvals),'or','MarkerFaceColor','r')
hold on;
plot(tvals,para_est(1,:)+2*pstd(1,:), '--', tvals,para_est(1,:)-2*pstd(1,:), '--')
xlabel('$$t$$','interpreter','latex','fontsize',28);
set(gca,'fontsize', 25);
title('$$S$$','interpreter','latex', 'fontsize', 25);

figure;
subplot(2,1,1);
title('State solutions comparison', 'fontsize', 25);
plot(tvals, xt(1,:),'r.', tvals, x_est(1,:),'b.');
set(gca,'fontsize', 25);
xlabel('$$t$$','interpreter','latex','fontsize',28);
ylabel('$$V(t)$$','interpreter','latex','fontsize',30);
legend('Truth', 'Estimation');
subplot(2,1,2);
plot(tvals, xt(2,:),'r.', tvals, x_est(2,:),'b.');
set(gca,'fontsize', 25);
xlabel('$$t$$','interpreter','latex','fontsize',28);
ylabel('$$m(t)$$','interpreter','latex','fontsize',30);
legend('Truth', 'Estimation');

figure;

subplot(2,1,1);
plot(tvals, abs(xt(1,:) - x_est(1,:))./xt(1,:), 'o-');
title('Errors in states', 'fontsize', 25);
xlabel('$$t$$','interpreter','latex','fontsize',30);
ylabel('$$V(t)$$','interpreter','latex','fontsize',30);
set(gca,'fontsize', 28);
subplot(2,1,2);
plot(tvals, abs(xt(2,:) - x_est(2,:))./xt(2,:), 'o-');
xlabel('$$t$$','interpreter','latex','fontsize',30); 
ylabel('$$m(t)$$','interpreter','latex','fontsize',30);
set(gca,'fontsize', 28);
figure;
subplot(2,1,1); plot(tvals,xt(1,:),'linewidth',3); 
xlabel('$$t$$','interpreter','latex','fontsize',30); ylabel('$$V(t)$$','interpreter','latex','fontsize',30);
set(gca,'fontsize',28);
title('True solutions','fontsize',26);
subplot(2,1,2);
plot(tvals,xt(2,:),'linewidth',3); 
xlabel('$$t$$','interpreter','latex','fontsize',30); ylabel('$$m(t)$$','interpreter','latex','fontsize',30);
set(gca,'fontsize',28);
