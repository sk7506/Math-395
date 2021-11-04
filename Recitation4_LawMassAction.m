% Law of mass action - enzyme kinetics
% Problem 1 - simulate kinetics of the enzyme
S0 = 10; % 1
E0 = 10; % 2
ES0 = 0; % 3
P0 = 0;  % 4
y0 =[S0;E0;ES0;P0];
kf = 1;
kr = 1;
kcat = 1;
odefun = @(t,y) EnzymeKinetics(t,y,kf,kr,kcat);
[t,y] = ode45(odefun, [0 2], y0);
% Forward Euler method
% Tf=10;
% dt=1e-3;
% nSt = Tf/dt;
% ys(1,:)=y0;
% for iT=1:nSt
%     ys(iT+1,:)=ys(iT,:)+dt*EnzymeKinetics(t,ys(iT,:),kf,kr,kcat)';
% end
% ts=(0:nSt)*dt;
% Michaelis-Menten approximation
slope=kcat*S0/(S0+kr/kf)*E0;
perror = ((y(end,4)-slope*t(end))/(y(end,4)))*100



function dydt = EnzymeKinetics(t,y,kf,kr,kcat)
    S = y(1);
    E = y(2);
    ES = y(3);
    P = y(4);
    dydt = [-kf*E*S + kr*ES; ...
        -kf*E*S + kr*ES + kcat*ES; ...
        kf*E*S - kr*ES - kcat*ES; ...
        kcat*ES];
end