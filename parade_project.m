
% X = [0,1,1,0];
% Y = [0,0,2,2];
% C = [1,0,0; 0,1,0; 0,0,1; 1,1,0];
% patch('faces', X,'vertices',Y,'FaceVertexCData',C,'EdgeColor','flat','FaceColor','none','LineWidth',2)

X = [0,1,1,0];
Y = [0,0,2,2];
C = [1,0,0;0 0 0;0 0 0;0 0 0];
patch('faces',X,'vertices',Y,'FaceVertexCData',C,'EdgeColor','flat','FaceColor','none','LineWidth',2)

%set(gca,'visible', 'off');
%grid on;
%set(gca,'visible', 'off');
% rectangle('Position',[0 0 1 1]);
% frontbumper = 
% axis([-10 10 -10 10]);

% % network of blocks (one-way streets) and intersections
% % modeled for a triangle of streets all going in the same direction
% % and ten cars in the loop
% 
% % b = block index
% % c = car index
% % iX(b) (X changes with the number of intersections)
% % ... = indices of intersections connected by block b
% % ni = # of intersections
% % nbin(i) = # of blocks entering intersection i
% % bin(i,j) = index of jth block entering intersection i for j = 1...nbin(i)
% 
% % nbout(i) = # of blocks leaving intersection i
% % bout(i, j) = index of jth block leaving intersection i
% 
% % three intersections:
% nc = 10;
% ni = 5; % number of intersections
% nimax = 5; % final intersection is the last one
% 
% % fixing number of roads going INTO an intersection
% nbin(1) = 0;  % tuned to one-way direction, the number of roads going into each 
% for i = 2:nimax
%     nbin(i) = 1; 
% end 
% % fixing number of roads going OUT of an intersection
% nbout(5) = 0; 
% for i = 1:4
%     nbout(i) = 1;
% end
% 
% 
% 
% for i = 1:ni
%     nbin(i) = sum(i2 == i);
%     nbout(i) = sum(i1 == i);
% end
% 
% nbinmax = max(nbin);
% nboutmax = max(nbout);
% bin = zeros(ni, nbinmax);
% bout = zeros(ni, nboutmax);
% for i = 1:ni
%     bin(i, 1:nbin(i)) = find(i2 == 1);
%     bout(i, 1:nbout(i)) = find(i1 == 1);
% end
% 
% % as a check it should be the case that sum(nbin) = sum(nbout) = nb
% 
% % let jgreen(i) be an integer designating which block has the green light
% % 1 <= jgreen(i) <= nbin(i)
% 
% % let s(b) be the state of the light at the end of block b
% % where s = 0 denotes red light & s = 1 denotes green
% 
% s = zeros(1, nb);
% for i = 1:ni
%     b = bin(i, jgreen(i));
%     s(b) = 1;
% end
% 
% % now let's expand a notion of distance
% % x(i) and y(i) are our coordinates of intersection i
% % L(b) = length of block b
% 
% %(ux(b), uy(b)) = unit vector along block b in the traffic flow direction
% 
% % given x(i) and y(i) we can find L, ux, and uy
% 
% ux = xi(i2) - xi(i1);
% uy = yi(i2) - yi(i1);
% 
% L = sqrt(ux.^2 + uy.^2);
% 
% ux = ux./L;
% uy = uy./L;
% 
% % now, let p(c) be the position of car c on whatever block it's on
% % if car c is on block b, then 0 <= p(c) < L(b)
% % coordinates of car c are given by 
% x(c) = xi(i1(b)) + p(c)*ux(b);
% y(c) = yi(i1(b)) + P(c)*uy(b);
% 
% % to access all cars on a block (decreasing p order) we use a linked-list structure:
% % firstcar(b) = index of first car on block b
% % nextcar(c) = index of car immediately behind car c on the SAME block
% % lastcar(b) = index of last car on block b
% 
% % if the block is empty then firstcar and lastcar = 0
% 
% % cars enter the roadway at random times and locations
% % let R be the rate of entry (units 1/(time*length))
% 
% % choose a timestep small enough that R*Lmax*dt is MUCH less than 1
% 
% if (rand < dt*R*L(b))
%     nc = nc + 1;
%     p(nc) = rand* L(b);
% end
% 
% % let bd(c) be the block on which the destination lies
% % let pd(c) be the position on block bd(c), expressed as a distance from the start of the block
% 
% bd(c) = 1 + floor(rand*nb);
% pd(c) = rand*L(bd(c));
% while (pd(c) >= L(bd(c)))
%     bd(c) = 1 + floor(rand*nb);
%     pd(c) = rand*Lmax;
% end

