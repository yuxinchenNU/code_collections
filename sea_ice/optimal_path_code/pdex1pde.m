function [c,g,S] = pdex1pde(t,s,X,DuDx)
    global ft f fx fxx sigma
    c = 1;
    g = 2*DuDx;
    S = -2*ft(X,t) - 2*f(X,t).*fx(X,t) - sigma^2*fxx(X,t);
end