% This routine solves the linear system 
% in Programming Assignment (a) 
% using a sparse storage for the coefficient matrix

% Dimension of the matrix
n = 1000;

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
rhs = 1.0/(n*n*n*n) * ones(n, 1);
tic
x = A\rhs;
toc

