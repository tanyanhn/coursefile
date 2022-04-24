function [x] = GaussEliminateSolver(A, r) 
n = size(A, 1);
L = eye(n, n);
for k = 1: n-1
    [~, p] = max(abs(A(k : n, k)));
    p = p + k - 1;
    if A(p, k) == 0 
        continue;
    end
    if p ~= k
        swapM(A, k, p, 1);
        swapM(L, k, p, 1);
        swapM(r, k, p, 1);
    end
    for i = k + 1: n
        L(i, k) = A(i, k) / A(k, k);
        A(i, :) = A(i, :) - L(i, k) * A(k, :); 
    end
end

x = zeros(n, 1);
y = zeros(n, 1);
for i = 1:n
    y(i) = (r(i) - L(i, :) * y);
end
for i = n:-1:1
    x(i) = (y(i) - A(i, :) * x) / A(i, i);
end

end

function [A] = swapM(A, i, j, dim) 
    if dim == 1
        tmp = A(i, :);
        A(i, :) = A(j, :);
        A(j, :) = tmp;
    elseif dim == 2
        tmp = A(:, i);
        A(:, i) = A(:, j);
        A(:, j) = tmp;
    end
end