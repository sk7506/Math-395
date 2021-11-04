%time = -0.1:0.001:0.1;
%frequency = 5;
%phase = 5;
%phase_in_rad = degtorad(phase);
%y = sin(2 * pi * frequency * time + phase_in_rad);
%plot(time, y), xlabel('Time'), ylabel('Sine wave')

%h = plot(y)
%direction = [1, 0];
%rotate(y, direction, 45)

for n = 0:10
    y = 10 - x;
    plot(x, y);
    hold on
end
drawnow
for n = 10:20
    y = 0;
    plot(x, y);
end
drawnow
