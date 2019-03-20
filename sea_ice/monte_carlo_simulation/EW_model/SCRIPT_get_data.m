% script to run EW model + noise to obtain tipping trajectories
% 
clearvars;
close all;

% load initial conditions
load('IC.mat'); % IC.mat contains initial condition for dF = [18, 19, 20, 21, 22]
E0_unstable = IC(:,2); % IC for unstable orbit
E0_stable1 = IC(:,1); % IC for lower stable orbit
E0_stable2 = IC(:,3); % IC for upper stable orbit
sigma = 4; % noise intensity
max_year = 10;
M0 = 144;% number of time steps in each year
randomSeeds = 1:3;
dir = 'Data';
mkdir(dir);
for j = 1:length(dF)
      for i = 1:length(randomSeeds)
         % solve dx = f(x,t)dt + sigma*dW using euler-maruyama method
         [t,E,E_unstable,E_s2,E_s1]=EW09_additive_noise(sigma,M0,max_year,dF(j),E0_unstable(j),E0_stable1(j),E0_stable2(j),randomSeeds(i));
         t = t';
         
         %%%%%% saving data
         dataName = [dir '/DF_' num2str(floor(dF(j))) 'p' num2str((dF(j)-floor(dF(j)))*100) '_randomseed_' num2str(randomSeeds(i)) '_sigma_' num2str(floor(sigma)) 'p' num2str((sigma-floor(sigma))*100) '.mat'];
         save(dataName,'t','E_unstable', 'E_s1', 'E_s2', 'E');
         
         %%%% plotting
         figure;
         figureName = [path '/DF_' num2str(floor(dF(j))) 'p' num2str((dF(j)-floor(dF(j)))*100) '_randomseed_' num2str(randomSeeds(i)) '_sigma_' num2str(floor(sigma)) 'p' num2str((sigma-floor(sigma))*100)];        
         plot(t,E,t,E_s2,t,E_unstable,t,E_s1);
         figureTitle = sprintf('dF = %g, Randomseed = %g, sigma = %g', dF(j),randomSeeds(i),sigma);
         title(figureTitle,'fontsize',14);
         set(gca,'fontsize',14);
         hold off;
      end
       
end

