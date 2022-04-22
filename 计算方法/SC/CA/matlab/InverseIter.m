A = [6 2 1;
    2 3 1;
    1 1 1];
mu = 2;
M = 5; % number of iterations
n = 3; % order of the matrix
v = zeros(n,1); v(1) = 1;
for k=1:M
    w = (A-mu*eye(n,n))\v; % apply (A-mu*I)^-1
    v = w/norm(w); % normalize
    lambda = v'*A*v; % Rayleight quotient
end
fprintf('The eigenvalue nearest to 2 is:\n  %.4d\n', lambda);
fprintf('a corresponding eigenvector is:\n');
v