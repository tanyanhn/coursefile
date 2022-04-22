% data points for interpolation
y = (0:1:8)';
t = y.^2;
u = linspace(0,64)'; % for plotting
% Perform cubic spline interpolation
v = spline(t,y,u);
plot(t,y,'o',u,v,'-')
axis([0 64 -01 8.01])
hold on
plot(u, sqrt(u), 'b')
hold off
legend('data points', 'cubic spline interpolation', 'sqrt(t)', 'Location', 'NorthWest')