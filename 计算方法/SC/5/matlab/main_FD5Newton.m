a = 0; b = 1; n = 80;
func_u = @(X, Y) sin(pi*X).*sin(pi*Y);
func_rhs = @(u) 2*pi*pi*u + u.*u.*u;
Newton_F = @(u) u.*u.*u; Newton_DF = @(u) 3*u.*u;
h = abs(b-a)/n; t = a + (0:n)*h; [xx, yy] = meshgrid(t); % build mesh
neqn = (n+1)*(n+1); u_old = zeros(neqn,1);
X = reshape(xx',neqn,1); Y = reshape(yy',neqn,1);
[A, idx_bnd, idx_inner] = spLaplacian(a,b,n);
NMaxNewtonIter = 100; TolNewton = 1e-8; iter = 1; error = 100.0;
while iter <= NMaxNewtonIter && error >= TolNewton
    f = Newton_F(u_old); df = Newton_DF(u_old);
    Mat = A + sparse(idx_inner, idx_inner, df(idx_inner), neqn, neqn);
    b = A*u_old + f - func_rhs(func_u(X,Y));
    b(idx_bnd) = func_u(X(idx_bnd), Y(idx_bnd)); % Dirichlet boundary
    u_new = u_old - Mat\b;
    error = norm(u_new - u_old, 2);
    % fprintf('%d: %20.15f\n', iter, error);
    u_old = u_new; iter = iter + 1;
end
mesh(xx', yy', reshape(u_new, n+1, n+1)); axis tight;