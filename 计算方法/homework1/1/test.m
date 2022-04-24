k = 20;
m = rand(k, 1);
for n = 5:2:k;

A = tril(-ones(n));
A = A + 2 * eye(n);
A(:, end) = 1;
x = m(1:n, 1);
r = A * x;
b = GaussEliminateSolver(A, r);
disp([ num2str(n), ' & ' , num2str(log(norm(b - x, 2)), 10), ' & ', num2str(log(norm(A * b - r, 2)), 10), ' & ', num2str(cond(A), 10), ' \\']);

end