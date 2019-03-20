% this is the script to get the right transition from lower stable orbit to
% the upper one. continue retrieving the data until one period after the
% period when tipping happens. Plot at least 3 periods.
clearvars;
close all;

alpha = 0.25;
beta = 0;
eps = 0.1;
A = 0.1;
sigma = 0.4;

M0 = 256;
% randomseed = 1:1000;
% randomseed = 1:8*10^4;
randomseed = 1:500;
% num_samp = 35441;
num_samp = 16000;

max_year = 20;
% note abs(alpha) < 1
if alpha < 0.1
    dir = ['../data/eps0p' num2str(eps*100) '/Alpha0p0' num2str(abs(alpha)*1000) '/sigma0p' num2str(sigma*1000)];
else
    dir = ['../data/eps0p' num2str(eps*100) '/Alpha0p' num2str(abs(alpha)*1000) '/sigma0p' num2str(sigma*1000)];
end
wid_frac = 12; % width fraction
throw_frac = 1/10;
y =1; % years to plot the transition trajectories
for j = 1:length(sigma)
    tot_tipOne = [];
    tot_tip = [];
    tippedE1 = [];%zeros(1,M0*y); % will store y periods of data where tipping occurs
    tippedE2 = [];
    tippedE3 = [];
    tippedE = [];
    time_escape1 = []; % store escaping time in first period
    time_escape2 = []; % store escaping time in second period
    time_escape3 = []; % store escaping time in third period
    meanEOne = zeros(1,M0*y);
    meanE = zeros(1,M0*y);
    fig1 = figure;
    for k = 1:length(randomseed)
        
        tPassS2 = []; % time when the path leaves the neighborhood of the initial orbit
        tPassS1 = []; % time when the path arrives the neighborhood of the final orbit
        TIP = [];
        [indS2,indS1,tip,t,E_s1,E_s2,E_unstable,E] = upTransitionBG( sigma(j),A,beta,alpha,randomseed(k),eps,dir,wid_frac,throw_frac );
        % indS2 is the index when the path leaves the neighborhood of the
        % initial orbit
        if ~isempty(tip)
            TIP = [TIP;tip'];
            tPassS2 = [tPassS2; t(indS2)];
            c = abs(rand(length(TIP),3));
            for i = 1:length(tPassS2)
                % get the starting index. want to start at the period before
                % tipping event happens
                startInd = floor(tPassS2(i))*M0+M0/2+1;
                if startInd > M0
                    if tPassS2(i) >= 0 && tPassS2(i) < 19 && startInd+M0-1 <= length(E) && E(startInd)<E_unstable(startInd) 
                        tippedE = [tippedE; E(startInd:startInd+M0-1)'];
                    end
                    if tPassS2(i) >= 0 && tPassS2(i) < 19 && startInd+M0-1 <= length(E) && E(startInd)>=E_unstable(startInd) 
                        tippedE = [tippedE; E(startInd - M0 + 1:startInd)'];
                    end
                    
                else
                    
                    if tPassS2(i) >= 0 && tPassS2(i) < 19 && startInd+M0-1 <= length(E)
                        tippedE = [tippedE; E(startInd:startInd+M0-1)'];
                    end
                end
                
            end
            
        end
    end
    
    time = t; 
    
    % pick the right number of bins in the histogram based on Freedman?Diaconis rule
    Nt = size(tippedE, 2); % number of time steps
    numBinSize = zeros(1,Nt);
    for i = 1:Nt
        param = iqr(tippedE(1:num_samp,i));
        n = size(tippedE(1:num_samp,:),1);
        numBinSize(i) = 2*param/n^(1/3);
    end
    Nb = floor((max(max(tippedE(1:num_samp,:))) - min(min(tippedE(1:num_samp,:))))/mean(numBinSize)); % number of bars in x
    % plot time-evolved histogram
    figure;
    counts = zeros(Nb, Nt);
    modeE = [];
    for i = 1:Nt
        [counts(:,i),centers] = hist(tippedE(1:num_samp,i),linspace(-1.5,1.5,Nb));
        counts(:,i) = counts(:,i)/sum(counts(:,i));
        [~, maxInd] = max(counts(:,i));
        modeE = [modeE, centers(maxInd)];
    end
    imagesc(time(1:M0*y),centers,counts(:,1:end)); colormap parula;
    set(gca,'YDir','normal')
    hold on;
    stdE = std(tippedE(1:num_samp,:));
    hold on;
    plot(time(1:M0*y),E_s1([M0/2+1:M0,1:M0/2]),time(1:M0*y),E_unstable([M0/2+1:M0,1:M0/2]),'--',time(1:M0*y),E_s2([M0/2+1:M0,1:M0/2]),'color','k','linewidth',4);
    
    ylim([-1.5,1.5]);
    colorbar;
        caxis([0 0.214])
    set(gca,'fontsize', 50);
    xlabel('$$t$$','interpreter','latex','fontsize',70)
    ylabel('$$x(t)$$','interpreter','latex','fontsize',70)
    set(gcf,'color','w');
    load('opEps0p1alpha0p25A0p1sigma0p4.mat');
    hold on;plot(t(1:M0),op(M0/2+1+M0*2:M0*3+M0/2),'r','linewidth',4);
    hold off
    % ==== plot on a cylinder
        p = max_year;
        figure;
        % plot histogram in the first period
        imagesc(time(round(M0/16):M0*1),centers,counts(:,round(M0/16):M0*1)); colormap parula;
        set(gca,'YDir','normal')
        hold on;
        for period = 2:p
            imagesc(time(M0*(period-1):M0*period),centers,counts(:,M0*(period-1):M0*period)); colormap parula;
            set(gca,'YDir','normal')
        end
        ylim([-2,2]);
        colorbar;
        caxis([0 0.45]);
    
        % ==== plot tipped trajectories on a cylinder
        p = max_year;
        nperiod = M0;
        figure;
        % plot histogram in the first period
        imagesc(time(round(M0/16):M0*1),centers,counts(:,round(M0/16):M0*1)); colormap parula;
        set(gca,'YDir','normal')
        hold on;
        for period = 2:p
        imagesc(time(M0*(period-1):M0*period),centers,counts(:,M0*(period-1):M0*period)); colormap parula;
        set(gca,'YDir','normal')
        end
        ylim([-2,2]);
        colorbar;
        caxis([0 0.45]);
    
    % =========== Convergence Test for Monte Carlo simulations
    numSamples = size(tippedE(1:num_samp,:), 1); % number of samples
    Nt = size(tippedE, 2); % number of time steps
    HN = zeros(Nb, Nt);
    for i = 1:Nt
        [HN(:,i),centers] = hist(tippedE(1:round(numSamples/2),i),linspace(-2,2,Nb));
        HN(:,i) = HN(:,i)/sum(HN(:,i));
    end
    
    % compute H_2N
    H2N = zeros(Nb, Nt);
    for i = 1:Nt
        [H2N(:,i),centers] = hist(tippedE(1:num_samp,i),linspace(-2,2,Nb));
        H2N(:,i) = H2N(:,i)/sum(H2N(:,i));
    end
    Error = 1/Nb/Nt*sum(sum((HN - H2N).^2));
end


