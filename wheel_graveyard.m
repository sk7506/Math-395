% drawing the ground
% plot([0, 10], [10, 10])
% hold on
% plot([10, 20], 20-[10, 20])
% plot([20, 30], [0, 0])
% drawing the ground

%R = 1;
%plot([-R+11 R+10 10])
%xlim([-10, 20])
%ylim([-10, 20])
%hold off

%for i = 1:N
%    Xcm(i,:) = Xcm(i,2) + 10;
%end

%         phi = pi/4; % 45 degree angle down
%         y = @(x) R/sin(phi) - x*cos(phi)/sin(phi);
%         fplot(y, [-5, 30]); 

%R = 1;
    %phi = pi/4; % 45 degree angle down MATLAB only does things in radians

    % the @ symbol specifies which variable the function is respect to
    %h = R/sin(phi) - x*cos(phi)/sin(phi); % gives a value, not a function handle  
    %grdh = [cos(phi) sin(phi) 0];
    
%h = x(2);
    %grdh = [0 1 0];
    
    % if the ground is an incline, the professor showed me that
    % we can write the ground as x*cos(phi)+y*sin(phi) = R
    
    % The single point that touches the ground and the circle is 
    % equal to R(-cos(phi), -sin(phi), 0)
    
    % The normal force along this incline is n = (cos(phi), sin(phi), 0)
    
    % h(x,y) piecewise: for x(0,10)  h = y-10
%                   for x(10,20) h = y+x
%                   for x(20,30) h = y
%function [h, grdh] = calcH(t) % this is the function that makes our ground
%    if (t(2) <= 10) % does t(2) have to equal 'step' instead?
%        h = t(2) - 10;
%        grdh = [0 1 0]; 
%        
%    elseif (10 <= t(2) && t(2) <= 20)
%        h = t(2) + t(1);
%        grdh = [1 1 0]/sqrt(2);
%        
%    else
%        h = t(2);
%        grdh = [0 1 0];
%        
%    end
%end