% driver for GF_pdepe
% taking different initial conditions and solving solutions to the gradient
% flow
% the stationary solutions are maximizers of the probability formulated via
% the path integral by Feyman

clearvars;
close all;
global E_u E_s1 E_s2 tOrbit
% here we take step functions transitioning at different times as our
% initial guesses
tmVec = linspace(0.5, 4.6, 20);
op = cell(length(tmVec),1);
hold on;
for i = 1:length(tmVec)
    [t, path] = GF_tm(tmVec(i));
    op{i} = path;
    plot(t, path, 'color', 'r', 'linewidth', 4);
    drawnow;
end

set(gca, 'fontsize', 28);

ylim([-1.5 1.5]);
hold on;
plot(tOrbit,E_s1,tOrbit,E_u,'--',tOrbit,E_s2,'color','k','linewidth',4);
ax = gca;
ax.XTick = ([0 1 2 3 4 5]);
ax.XTickLabel = ({'0','1','2','3','4','5'});
hold off;