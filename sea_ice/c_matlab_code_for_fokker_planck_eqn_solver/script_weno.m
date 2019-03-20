% solves fokker-planck equation for random walk with a time-varying
% external force, using WENO (a conservative PDE scheme), introduced in
% "High-order Finite Difference and Finite Volume WENO Schemes and
% Discontinuous Galerkin Methods for CFD"
% needs a conservative scheme to make sure prob = 1

clearvars;
close all;

% model params
A = 0.36;
epsilon = 1;
sigma = 0.8/sqrt(epsilon);
% sigma = sigma/sqrt(epsilon);
ic = 0.6;%-0.562123043682555; % initial condition
D = sigma^2/2; % diffusion coeff

% numerical params
% dt = 0.000025;
T = 10*2*pi; % terminal time
L = 4; % x\in [-L, L]
N0 = 401; % grid points in space
M = 1300000; 
savedM = M/100;
t = linspace(0,T,M); % time vector
dt = t(2) - t(1);
dx = 2*L/(N0-1); % grid size
const = D*dt/2/dx^2; % simplified notation used in CN step
N = N0 + 2; % add two ghost points for non-flux boundary condition

% set up grid 
x = linspace(-L,L,N0); % regular grid
x = [-L-dx,x,L+dx]; % including ghost point x0 and x_N
xh = x + 1/2*dx; % reconstructed grid, x{1/2}, ... x{N-1/2}, x{N+1/2}

% set up stencils
I = zeros(N, 2);
I(:,1) = xh(1:N);
I(:,2) = I(:,1) + dx;


% set initial condition
s = 0.1; % to approximate delta function
prob = 1/sqrt(2*pi)/s*exp(-(x-ic).^2/2/s^2); % use Gaussian
Prob = zeros(N,savedM);

drift =@(x,t) 1/epsilon*(x-x.^3 + A*cos(t)); % deterministic term
driftx = @(x) 1-3*x.^2; % x derivative of deterministic

% allocate variables
b = zeros(N,1); % vector to store rhs of the discretized equations

u = prob'; % store actual solution at grid points, N by 1, initialized to be gaussian
vm = zeros(length(xh),1); % fhat_minus at x_{i+1/2}
vp = zeros(length(xh),1); % fhat_plus at x_{i+1/2}
fhat = zeros(length(xh),1);

a = zeros(length(xh),1); % roe speed
%%%%%%%%%%%%%%%%%%%%%
% WENO reconstruction
% input: fu for finite difference; ubar_i for finite volume
% output: fhat for finite difference; u_{i+1/2} for finite volume
% define matrix for implicit solver, Crank Nicolson
vd = (1+2*const)*ones(1,N); % vector contains diagonal entries of M
vu = -const*ones(1,N-1); % vector contains upper diagonal entries of M
vl = -const*ones(1,N-1); % vector contains lower diagonal entries of M
Mat = diag(vd) + diag(vu,1) + diag(vl,-1);% store the diffusion matrix
% set up boundary conditions
Mat(1,1) = D/2/dx;
Mat(1,3) = -D/2/dx;
Mat(N,N) = -D/2/dx;
Mat(N,N-2) = D/2/dx;
%%% some vectors for WENO reconstruction
fuAppr = zeros(length(xh),3); % 3 approx in three substencils
% parameters used in the approximation from different neighboring
% substencils
c = [1/3, 5/6, -1/6; -1/6, 5/6, 1/3; 1/3, -7/6, 11/6];
    
fuApprt = zeros(length(xh),3); % 3 approx in three substencils
ct = [11/6, -7/6, 1/3; 1/3, 5/6, -1/6; -1/6, 5/6, 1/3];
for n = 2:M
    %%%%%% advection term update
    
    fu = drift(x,t(n)).*u'; % vector to store f(u_i), used to compute h^bar_i
    %%% step 1, compute the approximations from three stencils
    for i = 3:N-2
        % v(r)'s, eqn 2.51
        fuAppr(i,3) = c(1,:)*[fu(i-2); fu(i-1); fu(i)];
        fuAppr(i,2) = c(2,:)*[fu(i-1); fu(i); fu(i+1)];
        fuAppr(i,1) = c(3,:)*[fu(i); fu(i+1); fu(i+2)];
    end
    
    % step 2: obtain the k reconstructed values v_{i-1/2}^(r) using (2.10) based on
    % stencils (2.50)
    for i = 3:N-2
        % v(r)_{i-1/2}
        % fuApprt only has information at xh(2:N-3)
        fuApprt(i-1,1) = ct(1,:)*[fu(i); fu(i+1); fu(i+2)];
        fuApprt(i-1,2) = ct(2,:)*[fu(i-1); fu(i); fu(i+1)];
        fuApprt(i-1,3) = ct(3,:)*[fu(i-2); fu(i-1); fu(i)];
    end
    
    % step 3: find nonlinear weights (use the weights given in Jiang & Shu, JCP 1996)
    % k = 3
    d = [3/10, 3/5, 1/10]; % consider k = 3
    dT = d(end:-1:1);
    beta = zeros(1,3); % depending on x positions
    alpha = zeros(size(beta));
    w = zeros(size(beta));
    wt = zeros(size(beta));
    alphat = zeros(size(beta));
    eps = 10^-6; % used in determining nonlinear weights
    for i = 3:N-2
        beta(3) = 13/12*(fu(i-2) - 2*fu(i-1) + fu(i)).^2 + 1/4*(fu(i-2) - 4*fu(i-1) + 3*fu(i)).^2; % 2.63
        beta(2) = 13/12*(fu(i-1) - 2*fu(i) + fu(i+1)).^2 + 1/4*(fu(i-1) - fu(i+1)).^2;
        beta(1) = 13/12*(fu(i) - 2*fu(i+1) + fu(i+2)).^2 + 1/4*(3*fu(i) - 4*fu(i+1) + fu(i+2)).^2;
        alpha(1) = d(1)/(eps + beta(1))^2; % eqn 2.59
        alpha(2) = d(2)/(eps + beta(2))^2;
        alpha(3) = d(3)/(eps + beta(3))^2;
        w(1) = alpha(1)/(alpha(1)+alpha(2)+alpha(3)); % eqn 2.58
        w(2) = alpha(2)/(alpha(1)+alpha(2)+alpha(3));
        w(3) = alpha(3)/(alpha(1)+alpha(2)+alpha(3));
        % step 4 in procedure 2.2
        alphat(1) = dT(1)/(eps + beta(1))^2;
        alphat(2) = dT(2)/(eps + beta(2))^2;
        alphat(3) = dT(3)/(eps + beta(3))^2;
        wt(1) = alphat(1)/(alphat(1)+alphat(2)+alphat(3));
        wt(2) = alphat(2)/(alphat(1)+alphat(2)+alphat(3));
        wt(3) = alphat(3)/(alphat(1)+alphat(2)+alphat(3));
        % compute vminus{i+1/2} step 5
        % only has information at xh(3:N-2), need 2 stencils to compute
        % vm(2)
        vm(i) = w(1).*fuAppr(i,1) + w(2).*fuAppr(i,2) + w(3).*fuAppr(i,3);
        % compute vplus{i-1/2} step 5
        % only has information at xh(2:N-3), need 2 stencils to compute
        % vp(N-2)
        vp(i-1) = wt(1).*fuApprt(i-1,1) + wt(2).*fuApprt(i-1,2) + wt(3).*fuApprt(i-1,3);
    end
    
    %%% compute vm(2)
    % compute v(r)_{i+1/2}, i = 2
    c0 = [3/2, -1/2; 1/2, 1/2; -1/2, 3/2];
    v1 = c0(2,1)*fu(2) + c0(2,2)*fu(3); % r = 1;
    v2 = c0(3,1)*fu(1) + c0(3,2)*fu(2); % r = 2;
    % compute vm(2)
    d1 = 2/3;
    d2 = 1/3;
    beta1 = (fu(3)-fu(2))^2;
    beta2 = (fu(2)-fu(1))^2;
    alpha1 = d1/(eps + beta1)^2;
    alpha2 = d2/(eps + beta2)^2;
    w1 = alpha1/(alpha1+alpha2);
    w2 = alpha2/(alpha1+alpha2);
    vm(2) = w1*v1 + w2*v2;
    
    %%% compute vp(N-2)
    v1t = c0(1,1)*fu(N-1) + c0(1,2)*fu(N);
    v2t = c0(2,1)*fu(N-2) + c0(2,2)*fu(N-1);
    % compute vp(N-2)
    dt1 = 1/3;
    dt2 = 2/3;
    alphat1 = dt1/(eps+beta1)^2;
    alphat2 = dt2/(eps+beta2)^2;
    wt1 = alphat1/(alphat1+alphat2);
    wt2 = alphat2/(alphat1+alphat2);
    vp(N-2) = wt1*v1t + wt2*v2t;
    
     %%% compute vm(1) and vp(1) at xh(1), ie x(1+1/2) where x(1) is the
     %%% ghost point
    vm(1) = fu(1); % w1 = d1 = 1, c11 = 1
    vp(1) = fu(2); 
     
    %%% compute vm(N-1) and vp(N-1) at xh(N-1), ie x(N-1+1/2) where x(N) is the
    %%% ghost point
    vm(N-1) = fu(N-1); % w1 = d1 = 1, c11 = 1
    vp(N-1) = fu(N); 
    
    % Compute Roe speed
    for i = 1:N-1
        a(i) = (fu(i+1)-fu(i))/(u(i+1)-u(i));
    end
    
    % loop through the faces starting from 1 to N-1
    for i = 1:N-1
       if a(i) >= 0
           fhat(i) = vm(i);
       else
           fhat(i) = vp(i);
       end
    end
    
    
    %%%%%% diffusion term update
    % left boundary
    Mat(1,2) = drift(x(2),t(n));
    % right boundary
    Mat(N,N-1) = drift(x(N-1),t(n));
    
    B = Mat\eye(N);
    % update rhs of disretized system
    for i = 2:N-1
       b(i) = u(i) - dt/dx*(fhat(i) - fhat(i-1)) + const*(u(i+1) - 2*u(i) + u(i-1)); 
    end
    u = B*b;
    
    if mod(n,50)==0
        plot(x,u,'.');
        tit = sprintf('Total prob %e, time = %f',trapz(x,u), t(n));
        title(tit,'fontsize',16);
        drawnow;
        if sum(isnan(u)) ~= 0
            break;
        end
    end
    
    if mod(n, 100) == 0
        Prob(:,round(n/100)) = u;
    end
end


