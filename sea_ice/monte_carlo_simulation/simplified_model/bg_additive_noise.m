function [t1,E,E_unstable,E_s2,E_s1]=bg_additive_noise(E0,epsilon,sigma,param,param2,param3,M0,max_year,E0_unstable,E0_stable1,E0_stable2,randomseed)
% implement Euler-Maruyama method to obtain tipping trajectories for
% dx = f(x,t)dt + sigma*dW

global eps A beta alpha

rng(randomseed);
% === Parameters ===
% setting model parameters
eps = epsilon;
A = param;
beta = param2; % quadratic coeff
alpha = param3;
% numerical parameters
m = 10;
M = m*M0; % total time steps
 
% ====== integration backward to get unstable steady state
t0 = linspace(0,-1,M+1); % run for a year
dt = abs(t0(2) - t0(1)); % time step size
E_unstable = zeros(M,1);
E_unstable(1) = E0_unstable;

% forward Euler
for i = 1:M
    dE_unstable = -dEdtStable(t0(i),E_unstable(i)).*dt;
    E_unstable(i+1) = E_unstable(i) + dE_unstable;
end
E_unstable = E_unstable(1:m:end-1);
E_unstable = repmat(E_unstable,max_year,1);
E_unstable = [E_unstable;E0_unstable];
E_unstable = E_unstable(end:-1:1);
t0 = -t0; % reverse the order of time variable
% integrate to get lower stable orbit
E_s1 = zeros(M+1,1);
E_s1(1) = E0_stable1;

% forward Euler
for i = 1:M
    dE_stable1 = dEdtStable(t0(i),E_s1(i)).*dt;
    E_s1(i+1) = E_s1(i) + dE_stable1;
end
E_s1 = E_s1(1:m:end-1);
E_s1 = repmat(E_s1,max_year,1);
E_s1 = [E_s1;E0_stable1];

% ====== integration for getting perennially ice free stable state
E_s2 = zeros(M+1,1);
E_s2(1) = E0_stable2;

% forward Euler
for i = 1:M
    dE_stable2 = dEdtStable(t0(i),E_s2(i)).*dt;
    E_s2(i+1) = E_s2(i) + dE_stable2;
end
E_s2 = E_s2(1:m:end-1);
E_s2 = repmat(E_s2,max_year,1);
E_s2 = [E_s2;E0_stable2];


% === integration for solution with noise ===

% numerical parameters for Euler-Maruyama
N = M0*max_year;
E = zeros(N,1);

t1 = linspace(0,max_year,N+1);
dt1 = t1(2)-t1(1);
E(1) = E0; % initial condition
% Euler Maruyama method
% integration, M-1 iterations in time
    for i = 1:N
        dW = normrnd(0,sqrt(dt1));
        dE = dEdtStable(t1(i),E(i)).*dt1 + sigma.*dW;    
        E(i+1) = E(i) + dE;
        
    end
    
end

% set up RHS 
function dE = dEdtStable(t,E)
global eps A beta alpha
    dE = 1/eps*(E-E^3+beta*E^2+alpha+A*cos(2*pi*t));
end

