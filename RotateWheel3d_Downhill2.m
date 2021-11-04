M = 5;
N = 100; 
R = 1;
g = 9.8;
Mg = M*g;
mu = 10;
S = 10000;
Tf = 10;
% acceleration = -(2/3)*g*sin(pi/2) % negative because it's rolling downhill

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
    
    pointforce = zeros(N, 3); % evolve the force of each point using the calcForce function
    for points = 1:N
        pointforce(points, :) = calcForce(X(points, :), U(points, :), (M/N)*g, S, mu);
    end
    
    for points = 1:N % evolve the location of the points 
        Xtild(points, :) = rotate(Xtild(points, :), omega0, delta); % this is now a row vector
    end
    
    Ucm = Ucm + (delta/M)*(sum(pointforce)); 
    Xcm = Xcm + delta*Ucm;
    PointMasses(step+1,:) = Xcm;    
    L = L + delta*sum(cross(Xtild, pointforce))';
    
    plot( [X(:, 1); X(1,1)] , [X(:,2); X(1,2)] )
    hold on
    plot( X(3,1), X(3,2), 'ro')
    plot([0,10],[10,10])
    xlim([-10, 10])
    axis equal
    drawnow
    hold off
end
   
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