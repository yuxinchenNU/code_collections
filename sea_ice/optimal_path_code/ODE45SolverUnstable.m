function [t,X]=ODE45SolverUnstable(X0,Ain,betain,alphain,epsilon,years)

% This function numerically integrates the simplified model (smooth) for
% one period
% parameters are derived from mathematica notebook simpler_model_Oct17
global eps A alpha beta

% setting model parameters
eps = epsilon; % smoothing param
A = Ain;
alpha = alphain; % constant
beta = betain; % quadratic coeff
% set numerical parameters
% max_dur = 2*pi*3/eps;
% max_dur = -2*pi*years;
max_dur = -years;
RelTol=1e-8; % for ODE solver
AbsTol=1e-9; % for ODE solver


% === integration ===
X=[]; t=[];
options=odeset('RelTol',RelTol,'AbsTol',AbsTol);
% while intyr<=max_dur
    [tt,XX]=ode15s(@dEdt,[0 max_dur],X0,options);
    X=[X; XX]; t=[t; tt];
%     X = X(end:-1:1);
% end
% plot(t,E);
t = t+years;
% t = -t;
end

% set up RHS
function dX = dEdt(t,X)
global eps A alpha beta
%     dX = 1/eps*(X-X^3-beta*X^2-alpha+A*cos(t));
    dX = 1/eps*(X-X^3+beta*X^2+alpha+A*cos(2*pi*t));
end