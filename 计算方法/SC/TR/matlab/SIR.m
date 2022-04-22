N = 2000; % total population
beta = 0.08; gamma = 0.04; 
SIRfunc = @(t, y) [ -beta*y(2)/N*y(1);
    beta*y(2)/N*y(1)-gamma*y(2);
    gamma*y(2)];
t0 = 0; tfinal = 400;
% initial conditions
I0 = 20; S0 = N-I0; R0 = 0;
y0 = [S0; I0; R0];
[t, y] = ode45(SIRfunc,[t0,tfinal],y0);
plot(t,y(:,1),'-',t,y(:,2),'r-',t,y(:,3),'g-','LineWidth',3);
legend('Susceptible','Infected','Recovered')
title('Numerical result of SIR model')