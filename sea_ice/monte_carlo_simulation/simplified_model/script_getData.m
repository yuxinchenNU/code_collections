clearvars;
close all;

% model parameters
alpha = 0.05;
beta = 0;
eps = 0.01;
A = 0.3;
% A = 2*eps;
sigma = 0.1;

% find bifurcation point Ac
Ac1 = alpha + 2/9*sqrt(3+beta^2) + beta/27*(9+4*beta)*(beta+sqrt(3+beta^2));
Ac2 = 1/27*(27*alpha - 6*sqrt(3 + beta^2) + beta*(9 + 2*beta* (beta - sqrt(3 + beta^2))));
Ac = min(abs(Ac1),abs(Ac2));
fprintf('The thresholding A is %3.6g\n',Ac);

% First I need to get initial conditions for period orbits
years = 10;
% lower stable orbit
[t1,X1]=ODE45SolverStable(-1,A,beta,alpha,eps,years);
% upper stable orbit
[t2,X2]=ODE45SolverStable(1,A,beta,alpha,eps,years);
% unstable orbit
[t3,X3]=ODE45SolverUnstable(0,A,beta,alpha,eps,years);

E0_s1 = X1(end); % IC for lower stable orbit
E0_s2 = X2(end); % IC for upper stable orbit
E0_unstable = X3(end); % IC for unstable orbit
fprintf('Initial conditions are %3.12g, %3.12g, and %3.12g\n',X1(end),X2(end),X3(end));

M0 = 144;
randomseed = 1:4;
max_year = 200;
E0 = E0_s1; % initial condition for the stochastic run
% note abs(alpha) < 1
if alpha < 0.1
    path = ['../data/eps0p' num2str(eps*100) '/Alpha0p0' num2str(abs(alpha)*1000) '/sigma0p' num2str(sigma*1000)];
else
    path = ['../data/eps0p' num2str(eps*100) '/Alpha0p' num2str(abs(alpha)*1000) '/sigma0p' num2str(sigma*1000)];
end
mkdir(path);
c = abs(rand(length(randomseed),3));
for k = 1:length(randomseed)
    % obtain trajectories by sovling dx = f(x,t)dt + sigma*dW using
    % Euler-Maruyama method
	[t,E,E_unstable,E_s2,E_s1]=bg_additive_noise(E0,eps,sigma,A,beta,alpha,M0,max_year,...
                E0_unstable,E0_s1,E0_s2,randomseed(k));
            
	plot(t,E_s1,t,E_unstable,t,E_s2,t,E,'linewidth',2);
	figureName = [path '/eps' num2str(floor(eps)) 'p' num2str((eps-floor(eps))*100)...
                'sigma' num2str(floor(sigma)) 'p' num2str((sigma-floor(sigma))*100)...
                'A' num2str(floor(A)) 'p' num2str((A-floor(A))*100)...
                'rs' num2str(randomseed(k))];
    % save data        
	dataName = [figureName '.mat'];
	save(dataName);

end







