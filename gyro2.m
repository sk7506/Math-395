%gyroscope simulation
clear all
close all

g = 9.8      %gravitational acceleration (m/s^2)
omega = 10   %initial angular velocity (radians/s)

r = 1          %radius of wheel (m)
a = 0.5        %length of axle (m)
n = 20         %number of points on the rim
M_rim = 1      %total mass of the rim (Kg)
M_axle = 1     %total mass of the axle (Kg)
S_rim = 2000*(M_rim/n)*omega^2  %stiffness of each rim link (Kg/s^2)
D_rim = 0*(M_rim/n)*omega   %damping constant of each rim link (Kg/s)
S_spoke = S_rim %stiffnes of each spoke (Kg/s^2)
D_spoke = D_rim %damping constant of each spoke (Kg/s)
S_axle =  S_rim %stiffness of the axle (Kg/s^2)
D_axle =  D_rim %damping constant of the axle (Kg/s)

[kmax,lmax,X,jj,kk,S,D,Rzero,M] = ...
wheel(r,a,n,M_rim,M_axle,S_rim,D_rim,S_spoke,D_spoke,S_axle,D_axle)

%setup for animation
figure(1)
lr=1:n; %rim links
ls=n+1:3*n; %spoke links
la=3*n+1; %axle link
nskip=5

 hr=plot3([X(jj(lr),1),X(kk(lr),1)]',[X(jj(lr),2),X(kk(lr),2)]',...
          [X(jj(lr),3),X(kk(lr),3)]','r','linewidth',3)
 axis([-1.2*r,1.2*r,-1.2*r,1.2*r,-1.2*a,1.2*a])
 hold on
 hs=plot3([X(jj(ls),1),X(kk(ls),1)]',[X(jj(ls),2),X(kk(ls),2)]',...
          [X(jj(ls),3),X(kk(ls),3)]','b')
 ha=plot3([X(jj(la),1),X(kk(la),1)]',[X(jj(la),2),X(kk(la),2)]',...
          [X(jj(la),3),X(kk(la),3)]','k','linewidth',3)
 hold off
    axis equal
    axis manual
    axis(1.05*[-r,r,-r,r,-a,a])
    drawnow


fixed = n+1 %index of point that will be held fixed
forced = n+2 %index of point to which external force will be applied

tmax = 20   % duration of simulation (s)
clockmax = 10000 %number of time steps
dt = tmax/clockmax %(s)
X_save=zeros(clockmax,3);
t_save=zeros(clockmax,1)

t_ext_start = 2 %time that external force starts (s)
t_ext_stop =  3 %time that external force stops (s)
F_ext = 0.5*[1,0,0]*M_rim*g %external force (Kg.m/s^2)

%initial velocity
U = omega*[-X(:,2),X(:,1),zeros(kmax,1)];

for clock=1:clockmax
  t = clock*dt;
  DX = X(jj,:) - X(kk,:); %link vectors
  DU = U(jj,:) - U(kk,:); %link velocity difference vectors
  R = sqrt(sum(DX.^2,2)); %link lengths
  T = S.*(R-Rzero) + (D./R).*sum(DX.*DU,2); %link tensions
  TR=T./R; %link tensions divided by link lengths
  FF=[TR,TR,TR].*DX; %link force vectors

  F=zeros(kmax,3); %initialize force array for mass points 
  F(:,3) = - M*g;    %apply force of gravity to each link point

  % For each link, add its force with correct sign
  % to the point at each end of the link:
  for link=1:lmax
    F(kk(link),:)=F(kk(link),:)+FF(link,:);
    F(jj(link),:)=F(jj(link),:)-FF(link,:);
  end

  %apply external force during specified time interval:
  if((t_ext_start < t) && (t < t_ext_stop))
    F(forced,:) = F(forced,:) + F_ext;
  end

  U = U + dt*F./[M,M,M]; %update velocities of all points,
  %but if the index of the point is on the list "fixed",
  %set its velocity equal to zero:
  U(fixed,:)=0; 

  X = X + dt*U; %update positions of all points

  %store some results for future plotting
  X_save(clock,:)=X(n+2,:);
  t_save(clock) = t;
  
  %animation:
  if(mod(clock,nskip)==0)
    c=0;
    for l=lr
        c=c+1;
        hr(c).XData=[X(jj(l),1),X(kk(l),1)];
        hr(c).YData=[X(jj(l),2),X(kk(l),2)];
        hr(c).ZData=[X(jj(l),3),X(kk(l),3)];
    end
    c=0;
    for l=ls
        c=c+1;
        hs(c).XData=[X(jj(l),1),X(kk(l),1)];
        hs(c).YData=[X(jj(l),2),X(kk(l),2)];
        hs(c).ZData=[X(jj(l),3),X(kk(l),3)];
    end
    c=0;
    for l=la
        c=c+1;
        ha(c).XData=[X(jj(l),1),X(kk(l),1)];
        ha(c).YData=[X(jj(l),2),X(kk(l),2)];
        ha(c).ZData=[X(jj(l),3),X(kk(l),3)];
    end
    drawnow
    
  end
  
end
figure(2)
plot(t_save',X_save')
