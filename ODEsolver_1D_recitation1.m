x2(1) = 1; % goal: produce a big array of the values of x2 and u2
u2(1) = 1; % array indices in MATLAB begin with 1
r = 1; % u is velocity and x is position
s = 1; % increasing the stiffness also increases the amount of times the spring bounces
m = 1; % increasing the mass makes the spring swing higher out and with a longer period
d = 0; % increasing this value further dampens the spring
t = 0;
finalindex = 10/delt;
delt = 0.01;
for n = 1:finalindex % create a for loop for a specified range
    u2(n+1) = (delt/m)*((-s*(x2(n)-r) - d*u2(n))) + u2(n); %solves first equation for the new u2
    x2(n+1) = (u2(n+1)*delt) + x2(n);
end
plot(0:delt:10, x2)
hold on % prevents figure from closing and allows us to plot multiple lines on one figure


