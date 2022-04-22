%1-D non linear convection in conservative form: Simulating a travelling
...shock wave(a Heaviside function)
... d(u)/dt+u*d(u)/dx=0 ->d(u)/dt=-0.5*d(u^2)/dx

%Specifying parameters
nx=400;             %Number of steps in space(x)
nt=200;              %Number of time steps 
dt=0.01;            %Width of each time step
dx=4/(nx-1);        %Width of space step
x=0:dx:4;           %Range of x (0,4) and specifying the grid points
u=zeros(1,nx);      %Preallocating u
un=zeros(1,nx);     %Preallocating un
Eexp=0.11;          %Artificial viscosity

%%
%Initial conditions
index=find(x>=1);
u(1:index(1))=1;    %The Heaviside step function
ustar=u;
b=zeros(nx,1);
I=speye(nx,nx);

%%
i=2:nx-1;
%Calculating the velocity profile at every time step
for it=1:nt
    un=u;
    h=plot(x,u,'k');       %plotting the velocity profile
    axis([0 3 -0.5 2]); grid on;
    title({'1-D Convection';['time(\itt) = ',num2str(dt*it)]})
    xlabel('Spatial co-ordinate (x) \rightarrow')
    ylabel('Transport property profile (u) \rightarrow')
    drawnow; pause
    refreshdata(h)
    %Explicit methods:
    
    %Lax-Friedrichs(explicit, first order)
    %u(i) = 0.5*(un(i+1) + un(i-1)) - 0.25*(dt/dx)*(un(i+1).^2 - un(i-1).^2);
    h1 = 0.5*(0.5*un(i).^2+0.5*un(i+1).^2-norm(un,inf)*(un(i+1)-un(i)));
    h2 = 0.5*(0.5*un(i-1).^2+0.5*un(i).^2-norm(un,inf)*(un(i)-un(i-1)));
    u(i+1) = un(i) - (dt/dx)*(h1-h2);
end
        