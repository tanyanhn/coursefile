a = 0; b = 1; n = 1280;
func_u = @(x,y) (x.^2-x.^4).*(y.^4-y.^2);
func_rhs = @(x,y) 2*((1-6*x.^2).*y.^2.*(1-y.^2)+(1-6*y.^2).*x.^2.*(1-x.^2));
h = abs(b-a)/n; t = a + (0:n)*h;
neqn = (n+1)*(n+1); u_old = zeros(neqn,1);
[xx, yy] = meshgrid(t);
X = reshape(xx',neqn,1);
Y = reshape(yy',neqn,1);
[A, idx_bnd, idx_inner] = spLaplacian(a,b,n);
uu = func_u(X,Y); b = func_rhs(X,Y);
b(idx_bnd) = func_u(X(idx_bnd), Y(idx_bnd));
uh = A\b;
mesh(xx', yy', reshape(abs(uh-uu), n+1, n+1));
norm(uh-uu,inf)