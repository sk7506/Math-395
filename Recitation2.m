% add if statement: if object below the ground add this force
% add if statement: if utan is 0 the frictional force is 0
%force of gravity mg(0, -1) vector

%force as a function of x with velocity u equals force of gravity + spring
%force as a function of x and u + friction force as a function of x and u

% F(x, u) = Fg + (Fsp(x, u) + Ffric(x, u))delta for 
% delta = 1 if H(x) less than or equal to 0 and delta = 0 for H(x) greater than 0
%%
% % from last time:
% x2(1) = 1; % goal: produce a big array of the values of x2 and u2
% u2(1) = 1; % array indices in MATLAB begin with 1
% r = 1; % u is velocity and x is position
% s = 1; % increasing the stiffness also increases the amount of times the spring bounces
% m = 1; % increasing the mass makes the spring swing higher out and with a longer period
% d = 0; % increasing this value further dampens the spring
% t = 0;
% finalindex = 10/delt;
% delt = 0.01;
% for n = 1:finalindex % create a for loop for a specified range
%     u2(n+1) = (delt/m)*((-s*(x2(n)-r) - d*u2(n))) + u2(n); %solves first equation for the new u2
%     x2(n+1) = (u2(n+1)*delt) + x2(n);
% end
% plot(0:delt:10, x2)
% hold on % prevents figure from closing and allows us to plot multiple lines on one figure
%%
% goal: plot energy over time for when the object is above the ground
m = 1;
g = 9.8;
mg = m*g;
mu = 0;
S = 1000;
Tf = 2;
mu = 0;
u0 = [0,0];
x0= [0,1];
delta = 0.01;
finalindex = Tf/delta;
%%
%sym x y
% set up function that gives H and gradient of H
% input exactly what you want to output (same names)

for n = 1:finalindex % create a for loop for a specified range
    %call calcForce
    Force = calcForce(x0,u0,mg,S,mu);

    u0=u0+delta/m*Force; % evolves velocity

    %x(iT+1,:)=x(iT,:)+dt*u(iT+1,:);
    x0 = x0 + delta*u0;
    plot(x0(1),x0(2),'o')
    hold on
    plot(xlim,[0 0],'-k')
    ylim([-1 1])
    drawnow
    hold off

end
plot(0:delta:10, x0(2))
hold on

function [h, grdh] = calcH(x)
    h = x(2);
    grdh = [0 1];

    end

function Force = calcForce(x0,u0,mg,S,mu)
%H finds how far away from the ground our object is; call calcH
[h, grdh] = calcH(x0);
% delta: 1 if H(x) less than or equal to 0, 0 if H(x) greater than 0
delta = (h <=0);

%write equation for utan
utan = u0 - (grdh/norm(grdh))*dot(u0,grdh/norm(grdh));
utanhat = utan/norm(utan);

if (norm(utan)==0)
    utanhat = [0 0];
end

Fsp = S*(-h/norm(grdh))*(grdh/norm(grdh));
Ffric = mu*S*(-h/norm(h))*-utanhat;

% is Fg a constant? Fg is mg
Force = -mg*[0 1] + (Fsp + Ffric)*delta;
end


