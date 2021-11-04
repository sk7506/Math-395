R = 1;
phi = 45; % 45 degree angle down

% the @ symbol specifies which variable the function is respect to
y = @(x) R/sin(phi) - x*cos(phi)/sin(phi);

% fplot(f,xinterval) plots over the specified interval. 
% Specify the interval as a two-element vector of the form [xmin xmax].
fplot(y, [-10, 15]);
%axis equal

%for x = 1:10
%    plot(x, y);
%    xlim([-10,20])
%    axis equal
%end
