function [r1, r2] = exRoot(a, b, c)

bound = max(abs([a, b, c]));
a = a / bound;
b = b / bound;
c = c / bound;

m = sqrt(b * b - 4 * a * c);
if real(m) * b > 0
    if a > 10^-300
        r1 = (-b + m) / (2 * a);
    else
        r1 = (2 * c) / (-b - m);
    end
    if imag(r1) ~= 0
        r2 = abs(imag(r1));
        r1 = real(r1);
disp(["image result. r1 :", num2str(r1), "r2 :", num2str(r2)]);
        return;
    end
    r2 = (2 * c) / (-b + m);
else
    if  a > 10^-300
        r1 = (-b - m) / (2 * a);
    else
        r1 = (2 * c) / (-b + m);
    end
    if imag(r1) ~= 0 
        r2 = abs(imag(r1));
        r1 = real(r1);
disp(["image result. r1 :", num2str(r1), "r2 :", num2str(r2)]);
        return;
    end
    r2 = (2 * c) / (-b - m);
end
disp(["Real result. r1 :", num2str(r1), "r2 :", num2str(r2)]);
end