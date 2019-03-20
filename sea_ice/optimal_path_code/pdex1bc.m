function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
    global E0_stable2 E0_stable1
    pl = ul - E0_stable1;
    ql = 0;
    pr = ur - E0_stable2;
    qr = 0;
end