E = 1; epsilon = 1e-6; M = 20;
for k=0:M
    u = E - sin(E)/2 - 1;
    if abs(u)<epsilon
        break;
    end
    E = E - u/(1 - cos(E)/2);
end
fprintf('E: %d\n', E)