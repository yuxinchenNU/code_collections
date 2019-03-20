function [x] = EulerM_hist_theta_time(f,tvals,dt,xIC,sigma,theta)
% solving model using Euler Maruyama
tini = tvals(1);
M = (tvals(end)-tvals(1))/dt; %total number of time steps
N = (tvals(2)-tvals(1))/dt; %number of time steps between obs points
dimx = length(xIC);
x=zeros(dimx,M+1);
x(:,1) = xIC;

for n=1:M  
    k1 = dt * feval(f,tini+(n-1)*dt, x(:,n), theta(:,n));
    x(:,n+1) = x(:,n) + k1 + sigma*sqrt(dt).*randn(dimx,1);
end

x = x(:,1:N:end); %