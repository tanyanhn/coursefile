% data points for interpolation
y = (0:1:8)';
t = y.^2;
u = linspace(0,1)'; % for plotting
% Perform Lagrange polynomial interpolation
n = length(t)
v = zeros(size(u));
for k = 1:n
    w = ones(size(u));
    for j=[1:k-1 k+1:n]
        w = (u-t(j))./(t(k)-t(j)).*w;
    end
    v = v + w*y(k);
end
plot(u,v,'-')
axis([0 1 -.01 1.01])
hold on
% cubic spline interpolation
vv = spline(t,y,u);
plot(u,vv,'--r')
axis([0 1 -.01 1.01])
plot(u, sqrt(u), 'b')
hold off
legend('interpolating polynomial', 'cubic spline interpolation', 'sqrt(t)', 'Location', 'NorthWest')