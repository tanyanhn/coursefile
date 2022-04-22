% data points for interpolation
y = (0:1:8)';
t = y.^2;
u = linspace(0,64)'; % for plotting
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
plot(t,y,'o',u,v,'-')
axis([0 64 -5 100])
hold on
plot(u, sqrt(u), 'b')
hold off
legend('data points', 'interpolating polynomial', 'sqrt(t)', 'Location', 'NorthWest')