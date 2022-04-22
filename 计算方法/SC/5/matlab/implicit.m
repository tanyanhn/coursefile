% Solve U_t = \nu U_{xx}
% Set up parameters
N = 21;% number of grid points
a = 0; b = 1; nu = 1.0; T = 0.5;
h = (b-a)/(N-1); x = linspace(a,b,N); % space discretization
tau = 0.8*h*h/nu; % time step 
mu = nu*tau/h/h;
NT = ceil(T/tau); uh = zeros(N,NT+1);
uh(:,1) = sin(pi*x); % u_0 = sin(\pi x);
main = (1+2*mu)*sparse(ones(N-2,1));
off = -mu*sparse(ones(N-3,1));
A = diag(main) + diag(off,1) + diag(off,-1);
for n = 1:NT
    uh(2:N-1,n+1) = A\uh(2:N-1,n);
end
waterfall(uh'); xlabel('x'); ylabel('t');