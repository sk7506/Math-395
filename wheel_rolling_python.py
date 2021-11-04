#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 18:24:34 2021

@author: skilpatrick
"""

"""
dt = 0.01
lst_x = []
lst_y = []
t = 0
while t < 10: #for instance
    t += dt
    a = get_acceleration(function, x)
    x += v * dt + 0.5 * a * dt * dt
    v += a * dt
    y = get_position(fuction, x)
    lst_x.append(x)
    lst_y.append(y)
    
"""
"""
import numpy as np
import matplotlib.pyplot as plt

theta = np.arange(0, np.pi * 2, (0.01 * np.pi))
x = np.arange(-50, 1, 1)
y = x + 10
plt.figure()
for t in np.arange(0, 4, 0.1):
    plt.plot(x, y)
    xc = ((-9.81 * t**2 * np.sin(np.pi / 2)) / 3) + (5 * np.cos(theta))
    yc = ((-9.81 * t**2 * np.sin(np.pi / 2)) / 3) + (5 * np.sin(theta))
    plt.plot(xc, yc, 'r')
    xp = ((-9.81 * t**2 * np.sin(np.pi / 2)) / 3) + (5 * np.cos(np.pi * t))
    yp = ((-9.81 * t**2 * np.sin(np.pi / 2)) / 3) + (5 * np.sin(np.pi * t))
    plt.plot(xp, yp, 'bo')
    plt.pause(0.01)
    plt.cla()
plt.show()
"""

"""
# Simulating a wheel rolling on the ground
# Calculate moment of inertia tensor
M = 5;
N = 100;
r = 1;
g = 9.8;
S = 10000;
mu = 0;
Tf = 10;

thet = np.arrange(N-1)
np.thet.transpose
thet*2*pi/N

Omega = -1*np.array([0, 0, 1])

x = np.array([r*sin(thet), r*cos(thet), zeros(N,1)]) + np.array([0 r 0])

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

"""
import matplotlib as plt
import numpy as np

xs = np.array([0,1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20])
ys = np.array([])

for x in range(0,10):
    y = 10 - x
    ys.np.append(y)
    
for x in range(10,20):
    y = 0
    ys.np.append(y)
    
plt.plot(xs,ys)
plt.show()
    