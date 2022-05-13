J = @(x1,x2) [3+sin(x1),-cos(x2);-cos(x1),4+sin(x2)];
h = @(x1,x2)[3*x1-cos(x1) - sin(x2);4*x2-sin(x1) - cos(x2)];
x0 = [1.5;1]; J0 = J(x0(1),x0(2));
tempx1 = x0; tempJ = J0; err = 1; i = 0;
while err > 1e-5 && i <= 10
    i = i + 1;  x(:,i) = tempx1;
    s = tempJ\(-h(tempx1(1),tempx1(2)));
    tempx2 = tempx1 + s;
    yk = h(tempx2(1),tempx2(2)) - h(tempx1(1),tempx1(2));
    tempJ = tempJ + ((yk - tempJ*s)*s')/(s'*s); % Update J
    err = max(abs(h(x(1,i),x(2,i))));
    tempx1 = tempx2;       
 %   printf("Step %d: x = %f, f(x) = %f.\n",i, tempx1, yk);
end