function x = EulerM_theta(f,tvals,dt,x,sigma,theta)
tini = tvals(1);
M = ceil((tvals(end)-tvals(1))/dt); 
[x1,x2] = size(x);
for n=1:M  
    k1 = dt * feval(f,tini+(n-1)*dt, x, theta);
    x = x + k1 + sigma*sqrt(dt).*randn(x1,x2);
end