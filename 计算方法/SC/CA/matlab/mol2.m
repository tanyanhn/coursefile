global m h
a = 0; b = 1;
m = 320; h = (b-a)/(m+1);
x = linspace(a,b,m+2)';
t0 = 0; tfinal = 0.1;
y0 = cos(pi*x);
options = odeset('RelTol',1e-9,'AbsTol',1e-12);
[t,u] = ode23s('func_mol2',[t0,tfinal],y0,options);
plot(x,u(1,:));
hold on
grid
for i=2:length(t)
    plot(x, u(i,:));
end
uTrue = exp(-pi^2*tfinal)*cos(pi*x);
plot(x, uTrue, '*');
hold off
err = norm(u(end,:)'-uTrue, inf)