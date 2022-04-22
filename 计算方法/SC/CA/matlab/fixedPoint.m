E = 0; % initial guess
tol = 1e-6; % tolerance
while abs(E-1-sin(E)/2)>tol
    E = 1 + sin(E)/2; % fixed-point iteration formula
end
fprintf('E: %d\n', E)