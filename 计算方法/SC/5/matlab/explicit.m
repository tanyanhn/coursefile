N = 21; % Solve U_t = \nu U_{xx}
a = 0; b = 1; nu = 1.0; T = 0.5;
h = (b-a)/(N-1); x = linspace(a,b,N);
tau = 0.8*h*h/nu; mu = nu*tau/h/h;
NT = ceil(T/tau); uh = zeros(N,NT+1);
uh(:,1) = sin(pi*x); % u_0 = sin(\pi x);
for n = 1:NT
    for j = 2:N-1
        uh(j,n+1) = (1-2*mu)*uh(j,n) + ...
            mu*(uh(j-1,n) + uh(j+1,n));
    end
end
waterfall(uh'); xlabel('x'); ylabel('t');