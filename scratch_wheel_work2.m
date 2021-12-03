% initializing our variables and wheel components
M = 5;      % mass
N = 100;    % number of points on the wheel
R = 1;      % radius of the wheel
g = 9.8;    % downward force of gravity
Mg = M*g;   % mass times the downward force of gravity
mu = 10;    % friction constant (must be a positive number)
S = 100; % stiffness of the ground (must be positive)
Tf = 10000;   % value we will use to create timesteps
slope=-1;   % this is the slope of our hill

dthet = (2*pi)/N;       % unit circle equally divided by our number of points
thet = (0: N-1)'*dthet; % angle theta equally divided by our unit circle
delta = 8*1e-4;           % some small (adjustable) value to change timestep by
finalIndex = Tf/delta;  % 10000 divided by delta for our total timesteps
forcesSaved = zeros(finalIndex/1000, 3); % array that will save our forces for debugging purposes

omega0 = -[0 0 1]; % initial angular velocity 
% Position ball at (0,0) sitting on ramp
[h0, grdh0] = calcH([0 0 0], slope);
grdH0=grdh0/norm(grdh0);
Xcm = R*grdH0; % center of mass of wheel adjusted to put the wheel exactly at the sloped ground
thefig=figure;
movieframes=getframe(thefig);
                                                     
X = [R*cos(thet), R*sin(thet)+R, zeros(N,1)]+Xcm; % initial wheel constrained by the radius
Xtild = X - Xcm; % displacement between the wheel and the wheel's center of mass

Ucm = [0 0 0];                           % initial velocity
I = ((M*(R)^2)/2)*[1 0 0; 0 1 0; 0 0 2]; % initial moment of inertia matrix

L = I*omega0'; % angular momentum of the system about its center of mass

PointMasses = zeros(finalIndex+1, 3); % initial array of each points' mass
PointMasses(1,:) = Xcm;               % each mass updated by the system's center of mass
omega = zeros(finalIndex/1000, 1);    % overall angular velocity of the system at each point

step=1; % we begin with the first step, and while our timestep is less than the length..
        % ..of the hill, we keep the ball rolling
while (Xcm(1) < 200/sqrt(2)) % for each timestep, we want to compute a few things:
                             % 1) omega 2) X
                             % 3) U     4) new forces F given X and U
    omega0 = (I\L)';
    omega(step) = omega(3);
    X = Xtild + Xcm;         % force depends on position and translational velocity of each point
    U = Ucm + cross(omega0.*ones(N,3), Xtild);
    
    pointforce = zeros(N, 3); % evolve the force of each point using the calcForce function
    for points = 1:N
        pointforce(points, :) = calcForce(X(points, :), U(points, :), (M/N)*g, S, mu,slope);
        
        if (mod(step, 1000) == 0) % this helps save the force values for debugging purposes
            forcesSaved(step/1000, :) = forcesSaved(step/1000, :) + pointforce(points, :);
        end
    end
    
    for points = 1:N % evolve the location of the points 
        Xtild(points, :) = rotate(Xtild(points, :), omega0, delta); % this is now a row vector
    end
    
    Ucm = Ucm + (delta/M)*(sum(pointforce));      % update velocity of the center of mass
    Xcm = Xcm + delta*Ucm;                        % update the position of the center of mass
    PointMasses(step+1, :) = Xcm;            
    L = L + delta*sum(cross(Xtild, pointforce))'; % update angular momentum according to torque
    
    % fix the dimensions to make this work
    kineticEnergy = (1/2).*I.*omega.^2 + (1/2)*M.*(sum(Ucm.^2));
    potentialEnergy = Mg*(R); % change the R if you change the height with which you drop the wheel
    groundEnergy = calcForce(2);
    
    totalEnergy = kineticEnergy + potentialEnergy + groundEnergy;
    
    if (mod(step, 500) == 0)  % and now we plot! 
                            % to save time memory, we plot every 1000 timesteps
                            % this also helps the animation play more smoothly
        plot( [X(:, 1); X(1,1)] , [X(:,2); X(1,2)] )
        hold on
        plot( X(3, 1), X(3, 2), '.') % small red circle helps us see the wheel is spinning
        xlim([0 100])
        ylim([-50 10])
        plot([0 50],slope*[0 50])
        % axis equal
        pbaspect([2 1 1]) % this means the x-axis is twice as large as the other axes
        movieframes(length(movieframes)+1)=getframe(thefig); % append new frame to movieframes
        hold off
        plot(totalEnergy)
    end
    
    step=step+1;
end
%plot(forcesSaved)
movieframes(1)=[]; % this line deletes the first frame (blank frame)
v = VideoWriter('3DWheel_NoBounce.mp4', 'MPEG-4');
open(v);
writeVideo(v, movieframes);
close(v);

% plot(energy)

function [h, grdh] = calcH(x,slope) % this is the function that makes our ground
    % real life parameters make c = 1
    h = x(2)-slope*x(1);          % change the constant coefficient c to change the slope
    grdh = [-slope 1 0];          % this is the gradient of h
 
end

function [force, energy] = calcForce(x0,u0,mg,S,mu,slope) % changed for 3D
%H finds how far away from the ground our object is; call calcH
    [h, grdh] = calcH(x0,slope);
    
    force = -mg*[0 1 0];
    energy = 0;
     % make sure to zero out energy of the ground each time so it
     % doesn't accumulate
        
    % delta: 1 if H(x) less than or equal to 0, 0 if H(x) greater than 0
    if (h <= 0) % if the object is below the ground we introduce a restorative force
        
        % write equation for utan
        utan = u0 - (grdh/norm(grdh)) * dot(u0, grdh/norm(grdh));
        if (norm(utan) > 0.000000000000001)
            utanhat = utan/norm(utan);
        else
            utanhat = [0 0 0];
        end
        force = force + S*(-h/norm(grdh)*(grdh/norm(grdh)-mu*utanhat));
        energy = energy + (1/2)*S*(h/norm(grdh))^2; % one half times stiffness times distance squared
       
    end
end

function updatedx = rotate(x, omega0, delta) % Yields one row vector
    % Rodrigues Rotation Formula (figure 6 in the writeup)
    lomega = omega0*delta;
    omegahat = lomega/norm(lomega);
    updatedx = ((omegahat*x')*omegahat') + cos(norm(lomega))*(x'-(omegahat*x')*omegahat') + sin(norm(lomega))*(cross(omegahat, x))'; 
    updatedx = updatedx';
    % our first line gives us a column vector, and we would like to end up with a row vector
    % rotatedx = (number * column) + number * (column - (number * column)) + (number * column)
end