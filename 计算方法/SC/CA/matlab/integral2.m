% Evaluate the integral 
%   I_k = exp(-1)\int_0^1 x^k*exp(x)dx
% using the recurrence (unstable)
%   I_k = 1 - k*I_{k-1}
k = 21;
I = zeros(k,1);
fprintf(' k I_k\n');
I(1) = 1 - exp(-1);
for i=2:k
    I(i) = 1 - i*I(i-1);
    fprintf('%2d %.2e\n', i-1, I(i));
end