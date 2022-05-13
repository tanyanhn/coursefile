J = @(x1,x2) [3+sin(x1),-cos(x2);-cos(x1),4+sin(x2)];
h = @(x1,x2)[3*x1-cos(x1) - sin(x2);4*x2-sin(x1) - cos(x2)];
x0 = [1.5;1];  temp = x0; err = 1; i = 1;
while err > 1e-5 && i <= 10
    s = J(temp(1),temp(2))\(-h(temp(1),temp(2)));
    x(:,i) = temp + s;
    err = max(abs(h(x(1,i),x(2,i))));
    temp = x(:,i);  i = i+1;
end
