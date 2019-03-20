function X0 = pdex1ic(t)
global E_s2 E_s1 X1 X2 X3 X4 tm

X0 = E_s2(1)*(t>tm) + E_s1(1)*(t<=tm);
% X0 = E_s2(1);

% 

% X0(t<=0.5) = X1(t<=t(Nt/2));
% X0(t>t(Nt/2)&t<=t(Nt/2+tI)) = X2(t>t(Nt/2)&t<=t(Nt/2+tI));
% X0(t>t(Nt/2+tI)&t<=t(end-tI)) = X3(t>t(Nt/2+tI)&t<=t(end-tI));
% X0(t>t(end-tI)) = X4(t>t(end-tI));
% %                 X = (E_s2(TI(i)) - E0_stable2)/(tOrbit(TI(i)) - tOrbit(IC(k)))^6*(t-tOrbit(IC(k))).^6 + E0_stable2;
% X0 = X0';


end