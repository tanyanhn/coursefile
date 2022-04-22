k = 10;
I = speye(k,k);
E = sparse(2:k,1:k-1,1,k,k);
D = E+E'-2*I;
% spy(D)
A = kron(D,I)+kron(I,D);
spy(A)
B = kron(A,I)+kron(I,A);
% spy(B)
