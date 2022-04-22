function [A, idx_bnd, idx_inner] = spLaplacian(a,b,n)
h = abs(b - a)/n;  h2 = 1/(h*h);
nequation = (n+1)*(n+1);  nunknown = (n-1)*(n-1);
%% find out the boundary
idx_bnd = zeros(4*n,1);  idx_bnd(1:n) = 1:n;
idx_bnd(n+1:2*n) = n+1:n+1:n*(n+1);
idx_bnd(2*n+1 : 3*n) = (n+1)*(n+1):-1:n*(n+1)+2;
idx_bnd(3*n+1 : 4*n) = n*(n+1)+1: -(n+1) : n+2;
% and the inner boundary
idx_inner = setdiff((1:(n+1)*(n+1))', idx_bnd);
%% This is for the non-zero entries of the sparse matrix
ii = zeros(nunknown, 5);  jj = ii;  nnz = ii;  
ii(:,1) = idx_inner;  jj(:,1) = idx_inner;  nnz(:,1) = 4*h2; % diagnal one
ii(:,2) = idx_inner;  jj(:,2) = idx_inner + 1;  nnz(:,2) = -h2; % east 
ii(:,3) = idx_inner;  jj(:,3) = idx_inner - 1;  nnz(:,3) = -h2; % west 
ii(:,4) = idx_inner;  jj(:,4) = idx_inner - (n+1);  nnz(:,4) = -h2; % south 
ii(:,5) = idx_inner;  jj(:,5) = idx_inner + (n+1);  nnz(:,5) = -h2; % north 
A = sparse(ii,jj,nnz, nequation, nequation) + ...
        sparse(idx_bnd, idx_bnd, 1, nequation, nequation);
end