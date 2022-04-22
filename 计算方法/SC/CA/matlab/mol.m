% Solve the heat equation using 
% method-of-lines
global m h
a = 0; b = 1;
m = 5; h = (b-a)/(m+1);
x = linspace(a,b,m+2)';
t0 = 0; tfinal = 0.1;
y0 = sin(pi*x(2:end-1));
options = odeset('RelTol',1e-9,'AbsTol',1e-12);
[t,u] = ode23s('func_mol',[t0,tfinal],y0,options);
plot(x,[0 u(1,:) 0]);
hold on
grid
for i=2:length(t)
    plot(x, [0 u(i,:) 0]);
end
uTrue = exp(-pi^2*tfinal)*sin(pi*x);
plot(x, uTrue, '*');
hold off
err = norm([0 u(end,:) 0]'-uTrue, inf)