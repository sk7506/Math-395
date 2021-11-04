M = 5;
N = 100; 
R = 1;
g = 9.8;
Mg = M*g;
mu = 0;
S = 10000;
Tf = 10;

dthet = (2*pi)/N;
thet = (0: N-1)'*dthet;
delta = 1e-2;
finalIndex = Tf/delta; %10 divided by 0.01

omega0 = -[0 0 1];

X = [R*cos(thet), (R*sin(thet) + R), zeros(N,1)];
Xcm = [0 R 0];
Xtild = X - Xcm;

Ucm = [0 0 0];
I = ((M*(R)^2)/2)*[1 0 0; 0 1 0; 0 0 2];

L = I*omega0';

PointMasses = zeros(finalIndex+1, 3);
PointMasses(1,:) = Xcm;
omega = zeros(finalIndex, 1);

for step = 1:finalIndex
    omega0 = (I\L)';
    omega(step) = omega(3);
    X = Xtild + Xcm;
    U = Ucm + cross(omega0.*ones(N,3), Xtild);
    
    pointforce = zeros(N, 3);
    for points = 1:N
        pointforce(points, :) = calcForce(X(points, :), U(points, :), (M/N)*g, S, mu);
    end
    
    for points = 1:N
        Xtild(points, :) = rotate(Xtild(points, :), omega0, delta); % this is now a row vector
    end
    
    Ucm = Ucm + (delta/M)*(sum(pointforce));
    Xcm = Xcm + delta*Ucm;
    pointMasses(step+1,:) = Xcm;    
    L = L + delta*sum(cross(Xtild, pointforce))';
    
    plot( [X(:, 1); X(1,1)] , [X(:,2); X(1,2)] )
    hold on
    plot( X(3,1), X(3,2), 'ro')
    axis equal
    drawnow
    hold off
    
end
    
%%  (for step = 1:finalIndex)
% 1) compute:
%    omega(t) = I(t)*L(t)
     %omega = (I\L)';
%    Xk(t) = Xtildk(t) + Xcm(t)
     %X = Xtild + Xcm
%    Uk(t) = Ucm(t) + cross(omega, Xtildk)
     %U = Ucm + cross(omega0, Xtild)
%    Fk(t) = F(Xk(t), Uk(t)) compute force at each time step from each
     %pointforce = zeros(N, 3);
     %for points = 1:N
        %pointforce(points,:) = calcForce(X(points,:), U(points,:), mu, S, (M/N)*g);
     %end
%    point position and actual velocity (not u relative to center of mass)
% 2) rotate X:
%    Xtildk(t+delta) = rotate(Xtildk, omega*delta)
     %Xtild = rotate(Xtild, omega*delta)
% 3) update Ucm:
%    Ucm(t+delta) = Ucm(t) + (delta/M)*(sum(Fk(t))
     %Ucm = Ucm + (delta/M)*(sum(pointforce));
% 4) update Xcm:
%    Xcm(t+delta) = Xcm(t) + delta*Ucm(t+delta)
     %Xcm = Xcm + delta*Ucm
     %pointMasses(step+1,:) = Xcm
% 5) update L:
%    L(t+delta) = L(t) + delta*(sum(cross(Xtildk(t+delta)))
     %L = L + delta*sum(cross(Xtild, pointforce))';
  
%for n = 1:finalindex % create a for loop for a specified range
%call calcForce
%Force = calcForce(x0,u0,mg,S,mu);

%u0=u0+delta/m*Force; % evolves velocity

%x(iT+1,:)=x(iT,:)+dt*u(iT+1,:);
%x0 = x0 + delta*u0;

function [h, grdh] = calcH(x) % this is the function that makes our ground
    h = x(2);
    grdh = [0 1 0];

    end

function force = calcForce(x0,u0,mg,S,mu) % changed for 3D
%H finds how far away from the ground our object is; call calcH
    [h, grdh] = calcH(x0);
    force = -mg*[0 1 0];
    % delta: 1 if H(x) less than or equal to 0, 0 if H(x) greater than 0
    if (h <=0)
        %write equation for utan
        utan = u0 - (grdh/norm(grdh)) * dot(u0, grdh/norm(grdh));
        if (norm(utan) > 0.000000000000001)
            utanhat = utan/norm(utan);
        else
            utanhat = [0 0 0];
        end
        % Fg is mg
        force = force + S*(-h/norm(grdh)*(grdh/norm(grdh)-mu*utanhat));
    end
end

function updatedx = rotate(x, omega0, delta) % resolved issue about arguments! Yields one row vector

    lomega = omega0*delta;
    omegahat = lomega/norm(lomega);
    updatedx = ((omegahat*x')*omegahat') + cos(norm(lomega))*(x'-(omegahat*x')*omegahat') + sin(norm(lomega))*(cross(omegahat, x))'; 
    updatedx = updatedx';
    % our first line gives us a column vector, and we would like to end up with a row vector
    % rotatedx = (number * column) + number * (column - (number * column)) + (number * column)
end    