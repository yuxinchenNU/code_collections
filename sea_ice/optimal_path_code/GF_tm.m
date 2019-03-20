function [t, optimal_path] = GF_tm(tmIn)
% solves the gradient flow given an initial guesses
% tmIn is the time when the step function we take as the initial guess
% transitions

global A eps beta alpha D timesteps tm
global E_u E_s2 E_s1
global ft f fx fxx sigma
global E0_stable2 E0_stable1 tOrbit
% model parameters
alpha = 0.15;
beta = 0;
eps = 0.4;
A = 0.7;
sigma = 0.3;
D = sigma^2/2;
tm = tmIn;
% directory I used to save my data - change accordingly
dir = ['../results/eps0p' num2str(eps*1000)];
mkdir(dir);
dataname = [dir '/opAlpha0p' num2str(alpha*100) 'sigma0p' num2str(sigma*1000)];
figurename = [dir '/allPathsAlpha' num2str(alpha*100) 'sigma0p' num2str(sigma*1000)];
% years I take to get the optimal path
max_year = 5;
% number of time steps in each year
M0 = 2^9;

% functional terms to define F
U = @(x,t) -1/eps*(1/2*x.^2 - 1/4*x.^4 - 1/3*beta*x.^3 + alpha*x + A*cos(2*pi*t).*x);
f = @(x,t) 1/eps*(x - x.^3 - beta*x.^2 + alpha + A*cos(2*pi*t));
fx = @(x,t) 1/eps*(1 - 3*x.^2 - 2*beta*x);
fxx = @(x,t) 1/eps*(-6*x - 2*beta);
ft = @(x,t) 1/eps*A*2*pi*(-sin(2*pi*t));


% First I need to get IC for period orbits
years = 10;
% lower stable orbit
[~,X1]=ODE45SolverStable(-1,A,beta,alpha,eps,years);
% upper stable orbit
[~,X2]=ODE45SolverStable(1,A,beta,alpha,eps,years);
% unstable orbit
[~,X3]=ODE45SolverUnstable(0,A,beta,alpha,eps,years);

E0_stable1 = X1(end);
E0_stable2 = X2(end);
E0_unstable = X3(end);
% fprintf('Initial conditions are %3.12g, %3.12g, and %3.12g\n',X1(end),X2(end),X3(end));

% get deterministic orbits
[ tOrbit,E_u,E_s2,E_s1 ] = getDetOrbits(eps,A,beta,alpha,M0,max_year,E0_unstable,E0_stable1,E0_stable2);
tOrbit = tOrbit'; % convert to same dimension as E_s1

startTimeInd = M0*(max_year-2)+1+0*round(M0/3); % the start point of t_f domain
endTimeInd = M0*(max_year-1)+M0+1-0*round(M0/2); % the end point of t_f domain

IC = 1; % starting in the first period

tip = []; % vector to store tip time
Prob = [];
timesteps = 4; % time steps jump in picking t_f
TI = endTimeInd; % ending in the last period

Nt = M0; % time steps
Ns = 10^4; % steps for iteration
m = 0;
pathPI = zeros(Nt,length(TI)); % store all the paths for each corresponding ending time
% loop through initial positions
for k = 1:length(IC)
    % loop through end positions
    for i = 1:length(TI)
        
        % set up time
        t = linspace(tOrbit(IC(k)), tOrbit(end), Nt);
        % set up iteration
        s = linspace(0,20,Ns);
        
        % plot initial condition, setting up as step function
        X0 = E_s2(1)*(t>tm) + E_s1(1)*(t<=tm);
        plot(t,X0);
        
        dt = t(2) - t(1);
        
        % solving the gradient flow using pdepe
        sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,t,s);
        % Extract the first solution component as u.
        X = sol(:,:,1);
        optimal_path = X(end,:);
        
        % compute the energy
        Nx = length(X(end,:)); % dimension in x-direction
        % fourth order apporximation of xprime
        % need higher order because of sharp transition to improve energy computation accuracy
        dM = 1/12/dt*(diag(ones(Nx-2,1),-2) + diag(-8*ones(Nx-1,1),-1) - diag(ones(Nx-2,1),2) + diag(8*ones(Nx-1,1),1));
        dM(1,1:5) = [-25, 48, -36, 16, -3]/12/dt;
        dM(2,1:5) = [-3, -10, 18, -6, 1]/12/dt;
        dM(end, end-4:end) = [3,-16, 36, -48, 25]/12/dt;
        dM(end-1, end-4:end) = [-1, 6, -18, 10, 3]/12/dt;
        xprime = dM*X(end,:)';
        energyPI = trapz(t,(xprime - f(X(end,:)',t')).^2)+sigma^2*trapz(t,fx(X(end,:)',t'));
        pathPI(:,i) = X(end,:)';
        
    end
end


% ====== uncomment if saving the data
% save([dataname '.mat']);
% saveas(h,figurename, 'fig');
