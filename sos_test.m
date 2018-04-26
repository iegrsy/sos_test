clear;clc; 

ax = linspace(0,1,50);
ay = sin(ax*pi*2);

p = SOS(100,28,ax,ay);
y = polyval(p,ax);

subplot(2,1,1)
plot(ax,ay);
subplot(2,1,2)
plot(ax,y);