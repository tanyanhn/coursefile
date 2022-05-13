function yy=myfun(beta,x) %自定义拟合函数
yy=beta(1)+beta(2)*x(:, 1)+beta(3)*x(:, 2)+beta(4)*x(:, 3)+beta(5)*(x(:, 1).^2)+beta(6)*(x(:, 2).^2)+beta(7)*(x(:, 3).^2);
% yy=beta(1)+beta(2)*x(:, 1)+ beta(3)*(x(:, 1).^2) + beta(4)*(x(:, 1).^3);
end