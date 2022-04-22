n = 1000;

% set up the right-hand side
rhs = 1.0/(n*n*n*n)*ones(n, 1);

% backward substitution
rhs(n-1) = rhs(n-1) + 2*rhs(n);
for k = n-2:-1:1
    rhs(k) = rhs(k) + 2*rhs(k+1) - rhs(k+2);
end
rhs(1) = rhs(1)/2;

% forward substitution
rhs(1) = rhs(1)/2;
rhs(2) = rhs(2) + 2*rhs(1);
for k = 3:n
    rhs(k) = rhs(k) + 2*rhs(k-1) - rhs(k-2);
end

% set up the coefficient matrix
a = ones(n, 1);
b = -4*ones(n, 1);
b(end-1) = -2;
c = 6*ones(n, 1);
c(1) = 9; c(end-1) = 5; c(end) = 1;
d = -4*ones(n, 1);
d(end) = -2;
A = spdiags([a b c d a], -2:2, n, n);

% set up the right-hand side
rhs2 = 1.0/(n*n*n*n) * ones(n, 1);
x = A\rhs2;
norm(A*x-rhs2, Inf)
norm(A*rhs-rhs2, Inf)

disp(['inf-norm of the difference of the two solutions is ', num2str(norm(x-rhs, Inf))])