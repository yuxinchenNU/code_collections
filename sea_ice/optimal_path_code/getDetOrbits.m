function [ t,E_unstable,E_s2,E_s1 ] = getDetOrbits(epsilon,Ain,betain,alphain,M0,max_year,E0_unstable,E0_stable1,E0_stable2 )
% The function to get deterministic orbits for one cycle
    % === Parameters ===
% setting model parameters

global eps A alpha beta

eps = epsilon;
A = Ain;
alpha = alphain; % constant
beta = betain; % quadratic coeff
% numerical parameters
% max_year = 100;
m = 10;
M = m*M0;

% ====== integration backward to get unstable steady state
t0 = linspace(0,-1,M+1); % run for a year
% t0 = linspace(-1,0,M+1); % run for a year
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
t0 = -t0;
% t0 = t0+max_year; % reverse the order of time variable
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
    
t = linspace(0,max_year,length(E_s1));
% t = linspace(0,1,M0+1);
end

% set up RHS
function dE = dEdtStable(t,E)
global eps A beta alpha
    dE = 1/eps*(E - E^3 - beta*E^2 + alpha + A*cos(2*pi*t));
%     dE = 2*pi/eps*(E-E^3-beta*E^2-alpha+A*cos(2*pi*t));
end
