% This routine solves the linear system 
% in Programming Assignment (a) 
% using a dense storage for the coefficient matrix

% Dimension of the matrix
n = 100;

% set up the coefficient matrix
A = zeros(n, n);
for i=1:n-2
    A(i,i)=6;
    A(i,i+1)=-4; A(i+1,i)=-4;
    A(i, i+2)=1; A(i+2,i)=1;
end
A(1,1)=9; A(n-1,n-1)=5; A(n,n)=1;
A(n,n-1)=-2; A(n-1,n)=-2;

% set up the right-hand side
rhs = 1.0/(n*n*n*n)*ones(n, 1);
tic
x = A\rhs;
toc

disp(['The condition number of A is ', num2str(cond(A))])
