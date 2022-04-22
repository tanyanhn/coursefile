% Evaluate the integral 
%   I_k = exp(-1)\int_0^1 x^k*exp(x)dx
k = 21;
I = zeros(k,1);
fprintf(' k I_k\n');
for i=1:k
    I(i) = quad(@(x) x.^(i-1).*exp(x-1), 0, 1);
    fprintf('%2d %.4f\n', i-1, I(i));
end