% Simulating a wheel rolling on the ground
% Calculate moment of inertia tensor
M = 5;
N = 100;
r = 1;
g = 9.8;
S = 10000;
mu = 0;
Tf = 10;

thet = (0:N-1)'*2*pi/N;
Omega = -[0 0 1];

x = [r*sin(thet) r*cos(thet) zeros(N,1)]+[0 r 0];
x_cm = mean(x);
x_tilde = x-x_cm;

u_cm = [0 0 0];
I = zeros(3);

for k=1:N
    I = I+M/N*(norm(x_tilde(k,:))^2*eye(3)-x_tilde(k,:)'*x_tilde(k,:));
end

% I does not change for a circle
L = I*Omega';
dt = 1e-2;
nSteps=Tf/dt;

xcms = zeros(nSteps+1,3);
xcms(1,:)=x_cm;

Omegas = zeros(nSteps,1);

for iStep=1:nSteps
% Compute I and angular velocity from L
% I = zeros(3);
% for k=1:N
%     I = I+M/N*(norm(x_tilde(k,:))^2*eye(3)-x_tilde(k,:)'*x_tilde(k,:));
% end 
% I
% Compute forces on every pt (need real xs for this)
% Numerical method
Omega = (I \ L)';
Omegas(iStep)=Omega(3);
x = x_tilde+x_cm;
u = cross(Omega.*ones(N,3),x_tilde)+u_cm;
F = zeros(N,3);
for iPt=1:N
    F(iPt,:) = GroundForce3D(x(iPt,:),u(iPt,:),mu,S,M/N*g);
end
TotalF = sum(F);
% Rotate the Xks
for iPt=1:N
    x_tilde(iPt,:) = rotate(x_tilde(iPt,:),Omega*dt);
end
% Update center of mass
u_cm = dt/M*TotalF+u_cm;
L = L+dt*sum(cross(x_tilde,F))';
x_cm = x_cm+dt*u_cm;
xcms(iStep+1,:)=x_cm;
plot([x(:,1);x(1,1)],[x(:,2);x(1,2)])
hold on
plot(x(3,1),x(3,2),'ro')
plot(xlim,[0 0],'-k')
axis equal
drawnow
hold off
end

function force = GroundForce3D(x,u,mu,S,mg)
    [h,gradH] = Hvals3D(x);
    force = -mg*[0 1 0];
    if (h <= 0)
        % Add the ground forces
        n = gradH/norm(gradH); 
        Utan = u - dot(u,n)*n;
        if (norm(Utan) > 1e-10)
            UtanHat = Utan/norm(Utan);
        else
            UtanHat = [0 0 0];
        end
        force = force + S*(-h/norm(gradH)*(n-mu*UtanHat));
    end
end
    

function [h,gradH] = Hvals3D(x) % specifies the ground
    h = x(2);
    gradH = [0 1 0];
end