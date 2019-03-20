function [t1,E,E_unstable,E_s2,E_s1]=EW09_additive_noise(sigma,M0,max_year,dF,E0_unstable,E0_stable1,E0_stable2,randomseed)
% solve dE = f(E,t)dt + sigma*dW
% where details about f(E,t) can be found in Eisenman & Wettlaufer's 2009 paper
% called in SCRIPT_get_data.m to obtain tipping trajectories
global Fbot ai ao ki Li Hml cml F0input Fswinput FTinput tinput v0 tanha Tlin Tsmelt intmeth

rng(randomseed);
% === Parameters ===
% Measuring time in years for d/dt, using W/m^2 for fluxes and then
% multiplying d/dt by number of seconds in year.
yr=3.16e7; % seconds in year
Fbot=2;
ai=0.68;
ao=0.2;
ki=2;
Li=3*10^8/yr;
Hml=50;
cml=4*10^6/yr;
Tsmelt=0; % sfc temperature for onset of melt
tanha=0.5*Li; % if not equal to zero, gives width of albedo dependence
Tlin=0; % linearity of ice surface temperature: 0 for sea ice model, 1 for as mixed layer (linear)
v0=0.1; % ice export

% = For computing Ftop =
% atmmod: 0 for surface fluxes specified, 1 for active atmosphere with KLW
% from cloud fraction, 2 for specified F0(tinput) and FT(tinput).
atmmod=1;
tinput=(0.5:11.5)/12; % times that forcing is input at
% Stefan-Boltzmann linearization
sigma0=316;
sigmaT=3.9;
% surface fluxes from Maykut & Untersteiner (1971)
Fswinput=15.9*[0 0 1.9 9.9 17.7 19.2 13.6 9.0 3.7 0.4 0 0];
SH=15.9*[1.18 0.76 0.72 0.29 -0.45 -0.39 -0.30 -0.40 -0.17 0.10 0.56 0.79];
LH=15.9*[0 -0.02 -0.03 -0.09 -0.46 -0.70 -0.64 -0.66 -0.39 -0.19 -0.01 -0.01];
LWd=15.9*[10.4 10.3 10.3 11.6 15.1 18.0 19.1 18.7 16.5 13.9 11.2 10.9];
% atmospheric model parameters
% Cloud fraction from Maykut & Church (1973)
Cl=[4.9 4.7 4.2 5.9 8.1 8.0 8.3 8.8 8.9 8.4 7.1 5.2]/10;
% 0-70N mean temp based on NCEP-NCAR reanalysis
Tsouth=18*(1-1/3*cos((tinput-20/365)*2*pi));
kD=2.9;
tau0=0.5; % optical thickness independent of cloud fraction
tauc=3.6; % optical thickness dependence on cloud fraction

intmeth='linear'; % method for interpolating input data in time

% = atmospheric fluxes =
if atmmod==0 % specified surface fluxes
    F0input=-LH-SH-LWd+sigma0;
    FTinput=sigmaT+zeros(1,12);
elseif atmmod==1 % active equilibrium atmospheric model with KLW from cloud fraction, v0>0
    KLW=1./(tau0+tauc*Cl);
    F0input=KLW*sigma0-dF-kD*Tsouth/2;
    FTinput=KLW*sigmaT+kD/2;
elseif atmmod==2 % F0,FT specified as input
    F0input=F0;
    FTinput=FT;
end

% numerical parameters
m = 10;
M = m*M0;
t0 = linspace(0,1,M+1); % run for a year

dt = t0(2) - t0(1); % time step size

% ====== integration for getting unstable steady state
E_unstable = zeros(M,1);
E_unstable(1) = E0_unstable;

% forward Euler
for i = 1:M
    dE_unstable = dEdt(t0(i),E_unstable(i)).*dt;
    E_unstable(i+1) = E_unstable(i) + dE_unstable;
end
E_unstable = E_unstable(1:m:end-1);
E_unstable = repmat(E_unstable,max_year,1);
E_unstable = [E_unstable;E0_unstable];


% ====== integration for getting perennially ice stable state
E_s1 = zeros(M+1,1);
E_s1(1) = E0_stable1;

% forward Euler
for i = 1:M
    dE_stable1 = dEdt(t0(i),E_s1(i)).*dt;
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
    dE_stable2 = dEdt(t0(i),E_s2(i)).*dt;
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
E0 = E0_stable1; % initial value of E=-Li*hi+cml*Hml*Tml
E(1) = E0; % initial condition
% Euler Maruyama method
% integration, M-1 iterations in time
    for i = 1:N
        dW = normrnd(0,sqrt(dt1));
        dE = dEdt(t1(i),E(i)).*dt1 + sigma.*dW;    
        E(i+1) = E(i) + dE;
    
    end
    
end


% === model equations ===
function [F Tsrf Ftop F0 FT Fsw]=dEdt(t,E)
% see Eisenman & Wettlaufer's 2009 paper for model details
global Fbot cml Hml ki Li F0input Fswinput FTinput tinput ai ao intmeth v0 tanha Tlin Tsmelt intmeth
F0=interp1([tinput(end)-1 tinput tinput(1)+1],F0input([end 1:end 1]),mod(t,1),intmeth);
Fsw=interp1([tinput(end)-1 tinput tinput(1)+1],Fswinput([end 1:end 1]),mod(t,1),intmeth);
FT=interp1([tinput(end)-1 tinput tinput(1)+1],FTinput([end 1:end 1]),mod(t,1),intmeth);
Tsrf0=-(E.*F0-E*(1-ai).*Fsw)./(E.*FT-ki*Li); % ice surface temperature when no surface melt
Tsrfi=Tsrf0.*(Tsrf0<Tsmelt);
Tsrf=Tsrfi.*(E<0)+E/(cml*Hml).*(E>=0);

if Tlin==1 % for sea ice as ocean mixed layer with only albedo changing when Ts<0
    Tsrf=E/(cml*Hml).*(E<0)+E/(cml*Hml).*(E>=0);
end

if tanha>0 % using gradual albedo transition
    a=ao+(ai-ao)/2*(1-tanh(E/tanha));
else
    a=ai*(E<0)+ao*(E>=0);
end



F=-F0+(1-a).*Fsw-FT.*Tsrf+Fbot-v0*E.*(E<0);
Ftop=F0+a.*Fsw+FT.*Tsrf; % output for diagnostic plots
end
