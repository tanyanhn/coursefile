function up = func_mol2(t, u)
global m h
up = zeros(size(u));
up(1) = (-2*u(1)+2*u(2))/h^2;
for i=2:m+1
    up(i) = (u(i-1)-2*u(i)+u(i+1))/h^2;
end
up(m+2) = (2*u(m+1)-2*u(m+2))/h^2;
end