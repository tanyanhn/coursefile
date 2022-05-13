t = [0 1 3 4 6 7 9 10];
y = [8 6 5 2 1.5 1.3 1.1 1];
x = 0:0.2:10;
yy=interp1(t,y,x,'pchip');
plot(t,y,'o',x,yy,'linewidth',2)
