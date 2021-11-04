function [kmax,lmax,X,jj,kk,S,D,Rzero,M] = ...
wheel(r,a,n,M_rim,M_axle,S_rim,D_rim,S_spoke,D_spoke,S_axle,D_axle)

%This function constructs a wheel with n points on the rim,
%and 2 axle points, one on each side of the plane of the rim.
%Each rim point is linked to its two neighbors on the rim,
%and it is also linked to each axle point by a spoke.
%Finally, the two axle points are linked to each other.

%input parameters:

% r = radius of wheel
% a = length of axle
% n = number of points on the rim
% M_rim = total mass of the rim
% M_axle = total mass of the axle
% S_rim = stiffness of each rim link
% D_rim = damping constant of each rim link
% S_spoke = stiffnes of each spoke
% D_spoke = damping constant of each spoke
% S_axle = stiffness of the axle
% D_axle = damping constant of the axle

%outputs:
% X(k,:) = coordinates of point k
% jj(l),kk(l) = indices of points connected by link l
% S(l) = stiffness of link l
% D(l) = damping constant of link l
% Rzero(l) = rest length of link l

% There are n points on the rim and 2 axle points:
kmax = n+2;

%There are n rim links, 2*n spokes, and 1 axle link:
lmax = 3*n+1

R_spoke = sqrt(r^2 + (a/2)^2)  %length of each spoke
R_rim = 2*r*sin((2*pi)/(2*n))  %length of each rim link (chord length)

X=zeros(n+2,3);

%points on rim:
for k=1:n
  theta=2*pi*k/n;
  X(k,:) = r*[cos(theta),sin(theta),0];
end

%points at each end of axle:
X(n+1,:) = (a/2)*[0,0,-1];
X(n+2,:) = (a/2)*[0,0, 1];

M=[(M_rim/n)*ones(n,1);(M_axle/2)*ones(2,1)];  %mass of each point

jj=zeros(lmax,1);kk=zeros(lmax,1);
for k=1:n
  jj(      k) = k; kk(      k)= k+1-(k==n)*n; %link to next point on rim    
  jj(  n + k) = k; kk(  n + k)= n+1;          %link to one axle point
  jj(2*n + k) = k; kk(2*n + k)= n+2;          %link to other axle point
end
jj(3*n+1) = n+1; kk(3*n+1) = n+2;           %link between axle points

%stiffness, damping constant, and rest length of each link:
S=    [S_rim*ones(n,1);S_spoke*ones(2*n,1);S_axle];
D=    [D_rim*ones(n,1);D_spoke*ones(2*n,1);D_axle];
Rzero=[R_rim*ones(n,1);R_spoke*ones(2*n,1);a     ];



