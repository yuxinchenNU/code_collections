function [ indS,indE,tip,t,E_s1,E_s2,E_unstable,E ] = upTransitionBG( sigma,A,beta,alpha,k,epsilon,path,wid_frac,throw_frac)
% This function returns tipping transitions that go upward from computed
% data
% load data
% set data name
figureName = [path '/eps' num2str(floor(epsilon)) 'p' num2str((epsilon-floor(epsilon))*100)...
                'sigma' num2str(floor(sigma)) 'p' num2str((sigma-floor(sigma))*100)...
                'A' num2str(floor(A)) 'p' num2str((A-floor(A))*100)...
                'rs' num2str(k)];

dataname = [figureName '.mat'];

load(dataname);


% get DT
E = E(1:end-1); % throw out the last entry
t = t(1:end-1);

width1 = (E_unstable - E_s1)/wid_frac;
width2 = (E_s2 - E_unstable)/wid_frac;


index = []; % record when trajectory crosses unstable orbit
indStart = []; % record index when trajectories cross a neighborhood of the unstable orbit

% obtain index when trajectories cross unstable orbit
for l = 1:length(E)-1
    if (E(l) > E_unstable(l) && E(l+1) < E_unstable(l+1)) || (E(l)<E_unstable(l)&&E(l+1)>E_unstable(l+1))
        index = [index;l];
    end
    
end

% obtain index when trajectories cross a neighborhood of the unstable orbit
for l = 1:length(index)
    m = index(l);
    while ((E(m) <= E_s2(m) - width2(m)&& E(m)>E_unstable(m)) ||...
            (E(m) >= E_s1(m) + width1(m)&& E(m) < E_unstable(m)))&& m>1
        m = m-1;
    end
    indStart = [indStart;m];
end

% obtain index when trajectories leave a neighborhood of the unstable orbit
indEnd = zeros(size(indStart));

for l = 1:length(indStart)
    m = indStart(l);
    while m< length(E)&& ((E(m+1)< E_unstable(m+1)&&E(m+1)>=E_s1(m+1)+width1(m+1)) ||...
            (E(m+1)>E_unstable(m+1)&&E(m+1)<=E_s2(m+1)-width2(m+1)))
        m = m+1;
        
    end
    indEnd(l) = m;
end


% we want the last time when trajectories enter and leave the neighborhood
% of the unstable orbit. non-uniqueness comes from the fact that flow near
% the unstable orbit is weak
indS = unique(indStart);
indE = unique(indEnd);

rmInd = [];
% remove the non-tipping cases, sometimes trajectories cross the unstable
% orbit but tip back
for l = 1:length(indS)
    if (E(indS(l))-E_unstable(indS(l)))*(E(indE(l))-E_unstable(indE(l)))>=0
        rmInd = [rmInd;l];
    end
end

indS(rmInd) = [];
indE(rmInd) = [];

% remove down transtion
rmInd = [];
for l = 1:length(indS)
    if E(indS(l)) > E(indE(l))
        rmInd = [rmInd;l];
    end
end

indS(rmInd) = [];
indE(rmInd) = [];

zone = zeros(3,1);
dtInd1 = zeros(1,length(indS)); % keep track of the beginning and ending indices in zone(2), which is near unstable orbit
dtInd2 = zeros(2,length(indS));
for l = 1:length(indS)
    zone2ind = []; % keep track of indices in zone(2)
    zone1ind = [];
    zone3ind = [];
    for m = indS(l):indE(l)
        % obtain when entering neighborhood of the UPPER stable orbit
        if E(m) <= E_s2(m) - width2(m) && E(m) > E_unstable(m) + width2(m)
            zone(1) = zone(1) + 1;
            zone1ind = [zone1ind;m];
        else if E(m) <= E_unstable(m) + width2(m) && E(m) >= E_unstable(m) - width1(m)
                zone(2) = zone(2) + 1;
                zone2ind = [zone2ind;m];
                % obtain when entering neighborhood of the LOWER stable orbit
            else if E(m) < E_unstable(m) - width1(m) && E(m) >= E_s1(m) + width1(m)
                    zone(3) = zone(3) + 1;
                    zone3ind = [zone3ind;m];
                end
            end
        end
    end
    %     zone2ind
    dtInd1(l) = min(zone2ind);
    dtInd2(l) = max(zone2ind);
end

% calculate the time spent closed to unstable orbit for each tipping event
DT = zeros(size(indS));

for l = 1:length(indS)
    DT(l) = (t(dtInd2(l)) - t(dtInd1(l)))/2/pi;
    % throw out the transitions that hang around the unstable orbit for
    % too long (over 1/throw_frac of the period)
    if DT(l)>=1/throw_frac
        dtInd1(l) = 0;
        dtInd2(l) = 0;
        indS(l) = 0;
        indE(l) = 0;
    end
end

dtInd1(dtInd1==0) =[];
dtInd2(dtInd2==0) = [];
indS(indS==0) = [];
indE(indE==0) = [];

% Now get the time crosses E_unstable

if ~isempty(dtInd1)
    tip = zeros(1,length(dtInd1)); % record the time crosses E_unstable
    for l = 1:length(dtInd1)
        for m = dtInd1(l):dtInd2(l)
            if E(m) < E_unstable(m) && E(m+1) > E_unstable(m+1)
                tip(l) = t(m);
            end
        end
        
    end
else
    tip = [];
end

end

