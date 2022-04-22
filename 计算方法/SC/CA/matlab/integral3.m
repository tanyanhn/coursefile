% Evaluate the integral 
%   I_k = exp(-1)\int_0^1 x^k*exp(x)dx
% using the recurrence (unstable)
%   I_{k-1} = (1-I_k)/k
n = 22; % try different n
k = 21;
I = zeros(n,1);
fprintf(' k I_k\n');
for i=n:-1:2
    I(i-1) = (1 - I(i))/i;
end
I = [1-I(1); I];
% output the result
for i=1:k
    fprintf('%2d %.4f\n', i-1, I(i));
end
