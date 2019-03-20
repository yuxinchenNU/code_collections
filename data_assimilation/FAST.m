function [ Xdot ] = FAST( ~,X,theta)
% FAST model
% X: state variables
% theta contains all the model parameters
% see: A fast intensity simulator for tropical cyclone risk analysis,

    % parameters, same size with t
    Vp = theta(1,:); 
    S = theta(2,:); 
    Gamma = theta(3,:);
    hm = theta(4,:); % ocean mixed layer depth
    uT = theta(5,:);
    Cd = theta(6,:);
    h = theta(7,:);
    epsilon = theta(8,:);
    kappa = theta(9,:);
    beta = 1 - epsilon - kappa;
    
    % variables
    V = X(1,:);
    m = X(2,:);
    
    % some notation used in Emanuel 2017
    z = 0.01*Gamma.^(-0.4).*hm.*uT.*Vp./V;
    alpha = 1 - 0.87*exp(-z);
    gamma = epsilon + alpha.*kappa;
    
    % model equations
    Vdot = 1/2*Cd/h*(alpha.*beta.*Vp.^2.*m.^3 - (1 - gamma.*m.^3).*V.^2);
    mdot = 1/2*Cd/h*((1 - m).*V - 2.2*S.*m);
    
    Xdot = [Vdot; mdot];
end

