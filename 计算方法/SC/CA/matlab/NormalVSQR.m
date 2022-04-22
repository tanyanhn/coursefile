m = 21; % number of data points
n = 11; % degree of the interpolating polynomial
epsilon = 1e-10; % perturbation parameter
t = linspace(0, 1, m)'; % equally spaced data points on [0,1]
y = zeros(size(t)); % polynomial values at t_i's
for k=0:n
    y = y + t.^k;
end
u = rand(m,1); % random numbers distributed in (0,1)
y = y + (2*u-1)*epsilon;
A = zeros(m, n+1);
for i=1:n+1
    A(:,i) = t.^(i-1);
end

% Normal equations approach
R = chol(A'*A); % Cholesky factorization
rhs = A'*y;
z = R'\rhs; % forward substitution
x = R\z; % backward substitution
fprintf('Normal equations approach: the inf-norm of the error is: \n%d\n', norm(x-ones(n+1,1), inf));

% QR factorization approach
[Q, R] = qr(A); % QR factorization
rhs = Q'*y;
x = R(1:n+1,:)\rhs(1:n+1); % backward substitution
fprintf('QR factorization approach: the inf-norm of the error is: \n%d\n', norm(x-ones(n+1,1), inf));