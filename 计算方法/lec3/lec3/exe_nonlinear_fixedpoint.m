f = @(x1,x2) [(cos(x1) + sin(x2))/3;(sin(x1) + cos(x2))/4];
h = @(x1,x2) [3*x1-cos(x1) - sin(x2);4*x2-sin(x1) - cos(x2)];
x0 = [1.5;1]; err = 1; temp = x0; i = 1;
while err>1e-5
    x(:,i) = f(temp(1),temp(2));
    err = max(abs(h(x(1,i),x(2,i))));
    temp = x(:,i);     i = i+ 1;
end
