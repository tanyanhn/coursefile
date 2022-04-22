function cavity_simple
%%LID DRIVEN CAVITY FLOW EXAMPLE
%%The code is based on the SIMPLE algorithm

format longG

%% GRID SIZE AND OTHER PARAMETERS
%i runs along x-direction and j runs along y-direction 

imax=33;                        %grid size in x-direction 
jmax=33;                        %grid size in y-direction 
max_iteration=6000; 
maxRes = 1000;
iteration = 1;
mu = 0.01;                      %viscosity
rho = 1;                        %density
velocity=1;                     %velocity = lid velocity
dx=1/(imax-1);					%dx,dy cell sizes along x and y directions
dy=1/(jmax-1); 
x=dx/2:dx:1-dx/2; 
y=0:dy:1; 
alphaP = 0.1;                   %pressure under-relaxation
alphaU = 0.7;                   %velocity under-relaxation
tol = 1e-5;

%   u_star, v_star are Intermediate velocities
%   u and v = Final velocities

%Variable declaration
p   = zeros(imax,jmax);             %   p = Pressure
p_star   = zeros(imax,jmax);        
p_prime = zeros(imax,jmax);         %   pressure correction 
rhsp = zeros(imax,jmax);            %   Right hand side vector of pressure correction equation
divergence = zeros(imax,jmax); 

%Vertical velocity
v_star = zeros(imax,jmax+1);
vold   = zeros(imax,jmax+1);
vRes   = zeros(imax,jmax+1);
v      = zeros(imax,jmax+1);
d_v    = zeros(imax,jmax+1);    %velocity orrection coefficient

% Horizontal Velocity -----------
u_star = zeros(imax+1,jmax);
uold   = zeros(imax+1,jmax);
uRes   = zeros(imax+1,jmax);
u      = zeros(imax+1,jmax);
d_u    = zeros(imax+1,jmax);  %velocity orrection coefficient

%Boundary condition 
%Lid velocity (Top wall is moving with 1m/s)
u_star(1:imax+1,jmax)=velocity;
u(1:imax+1,jmax)=velocity;

%% ---------- iterations -------------------
while ( (iteration <= max_iteration) && (maxRes > tol) ) 
    iteration = iteration + 1;
    [u_star,d_u] = u_momentum(imax,jmax,dx,dy,rho,mu,u,v,p_star,velocity,alphaU);       %%Solve u-momentum equation for intermediate velocity u_star 
    [v_star,d_v] = v_momentum(imax,jmax,dx,dy,rho,mu,u,v,p_star,alphaU);                 %%Solve v-momentum equation for intermediate velocity v_star
    uold = u;
    vold = v; 
    
    rhsp=zeros((jmax-2)*(imax-2),1);    %vector of RHS for solving pressure corrections
    stride = jmax;

   % RHS is the same for all nodes except the p_prime(1,1)
   % because p(1,1) is set to be zero, it has no pressure correction
    for j=1:jmax
        for i=1:imax
            position = i + (j-1)*stride; 
            rhsp(position) = rho * (u_star(i,j)*dy - u_star(i+1,j)*dy + v_star(i,j)*dx - v_star(i,j+1)*dx); 
        end
    end

    % modify for p_prime(1,1)
    rhsp(1,1) = 0;    
    %%Calculate rhs vector of the Pressure Poisson matrix 
    N = imax*jmax;
    stride = jmax;
    Ap = zeros(N,N);
    %P-prime for the boundary nodes is set to zero
    %the interior nodes imax-2*jmax-2 solved implicitly for P-prime
    for j=1:jmax
        for i=1:imax
            position = i + (j-1)*stride; 
            aE = 0;
            aW = 0;
            aN = 0;
            aS = 0;
            %set BSc for four corners
            if(i == 1 && j == 1)
                Ap(position,position) = 1;                 
                %pressure correction at the first node is zero
                continue;
            end
            if (i == imax && j == 1)
                Ap(position,position-1) = -rho*d_u(i,j)*dy;                    
                aW = -Ap(position,position-1);
                Ap(position,position+stride) = -rho*d_v(i,j+1)*dx;       
                aN = -Ap(position,position+stride);
                aP = aE + aN + aW + aS;
                Ap(position,position) = aP; 
                continue;
            end
            if (i == 1 && j == jmax)
                Ap(position,position+1) = -rho*d_u(i+1,j)*dy;                    
                aE = -Ap(position,position+1);
                Ap(position,position-stride) = -rho*d_v(i,j)*dx;  
                aS = -Ap(position,position-stride);
                aP = aE + aN + aW + aS;
                Ap(position,position) = aP;                   
                continue;
            end
            if (i == imax && j == jmax)
                Ap(position,position-1) = -rho*d_u(i,j)*dy;                    
                aW = -Ap(position,position-1);
                Ap(position,position-stride) = -rho*d_v(i,j)*dx;       
                aS = -Ap(position,position-stride);
                aP = aE + aN + aW + aS;
                Ap(position,position) = aP;
                continue;
            end
            %set four boundaries
            if (i == 1)
                Ap(position,position+1) = -rho*d_u(i+1,j)*dy;                    
                aE = -Ap(position,position+1);
                Ap(position,position+stride) = -rho*d_v(i,j+1)*dx;       
                aN = -Ap(position,position+stride);
                Ap(position,position-stride) = -rho*d_v(i,j)*dx;       
                aS = -Ap(position,position-stride);
                aP = aE + aN + aW + aS;
                Ap(position,position) = aP;
                continue;
            end
            if (j == 1)
                Ap(position,position+1) = -rho*d_u(i+1,j)*dy;                    
                aE = -Ap(position,position+1);
                Ap(position,position+stride) = -rho*d_v(i,j+1)*dx;       
                aN = -Ap(position,position+stride);
                Ap(position,position-1) = -rho*d_u(i,j)*dy;                    
                aW = -Ap(position,position-1);
                aP = aE + aN + aW + aS;
                Ap(position,position) = aP;
                continue;
            end
            if (i == imax)
                Ap(position,position+stride) = -rho*d_v(i,j+1)*dx;       
                aN = -Ap(position,position+stride);
                Ap(position,position-stride) = -rho*d_v(i,j)*dx;       
                aS = -Ap(position,position-stride);
                Ap(position,position-1) = -rho*d_u(i,j)*dy;                    
                aW = -Ap(position,position-1);
                aP = aE + aN + aW + aS;
                Ap(position,position) = aP;
                continue;
            end
            if (j == jmax)
                Ap(position,position+1) = -rho*d_u(i+1,j)*dy;                    
                aE = -Ap(position,position+1);
                Ap(position,position-stride) = -rho*d_v(i,j)*dx;       
                aS = -Ap(position,position-stride);
                Ap(position,position-1) = -rho*d_u(i,j)*dy;                    
                aW = -Ap(position,position-1);
                aP = aE + aN + aW + aS;
                Ap(position,position) = aP;
                continue;
            end
            % interior nodes
            Ap(position,position-1) = -rho*d_u(i,j)*dy;                          %sub diagonal
            aW = -Ap(position,position-1);

            Ap(position,position+1) = -rho*d_u(i+1,j)*dy;                        %%upper diagonal
            aE = -Ap(position,position+1);

            Ap(position,position-stride) = -rho*d_v(i,j)*dx;                     %%sub sub diagonal
            aS = -Ap(position,position-stride);

            Ap(position,position+stride) = -rho*d_v(i,j+1)*dx;                   %%upper upper diagonal
            aN = -Ap(position,position+stride);
        
            aP = aE + aN + aW + aS;
            Ap(position,position) = aP;                                          %%main diagonal
        end
    end
    %%Form the Pressure Poisson coefficient matrix
    pressure = p_star;                                                               %   p = Pressure
    p_prime = zeros(imax,jmax);                                                 %   pressure correction 

    p_prime_interior = rhsp\Ap;%

    %convert pressure correction in to a matrix
    %update preesure values
    z=1; 
    for j=1:jmax
        for i=1:imax
            p_prime(i,j)=p_prime_interior(z); 
            z=z+1;
            pressure(i,j) = p_star(i,j) + alphaP*p_prime(i,j);
        end
    end
    pressure(1,1) = 0;    
    %%Solve pressure correction implicitly and update pressure
    v = zeros(imax,jmax+1);
    u = zeros(imax+1,jmax);

    %update interior nodes of u and v
    for i=2:imax
        for j=2:jmax-1
        
            u(i,j) = u_star(i,j) + d_u(i,j)*(p_prime(i-1,j)-p_prime(i,j));
        
        end
    end

    for i=2:imax-1
        for j=2:jmax
        
            v(i,j) = v_star(i,j) + d_v(i,j)*(p_prime(i,j-1)-p_prime(i,j));
        
        end
    end

    %update BCs
    v(1,1:jmax+1) = 0.0; %left wall
    v(imax,1:jmax+1) = 0.0; %right wall
    v(1:imax, 1) = -v(1:imax, 2); %bottom wall
    v(1:imax, jmax+1) = -v(1:imax, jmax); %top wall 

    u(1,1:jmax) = -u(2,1:jmax); %left wall
    u(imax+1,1:jmax) = -u(imax,1:jmax); %right wall
    u(1:imax+1, 1) = 0.0; %bottom wall
    u(1:imax+1, jmax) = velocity; %top wall             %%Update velocity based on pressure correction                  
    p_star = p;                                                                          %%use p as p_star for the next iteration
    
    %find maximum residual in the domain
    vRes = abs(v - vold);
    uRes = abs(u - uold);
    maxRes_u = max(max(uRes));
    maxRes_v = max(max(vRes));
    maxRes = max(maxRes_u, maxRes_v);
                                                                            %%Check for convergence 
    disp(['It = ',int2str(iteration),'; Res = ',num2str(maxRes)])
    if (maxRes > 2)
        disp('not going to converge!');
        break;
    end
end
%% plot


disp(['Total Iterations = ',int2str(iteration)])

figure 
contourf(x,y,u(2:imax,:)',50, 'edgecolor','none');colormap jet
colorbar;
axis([0 1 0 1]); 
title('steady Ux'); 

end


function [u_star,d_u] = u_momentum(imax,jmax,dx,dy,rho,mu,u,v,p,velocity,alpha)

u_star=zeros(imax+1,jmax);
d_u=zeros(imax+1,jmax);

De  = mu*dy / dx;  %convective coefficients
Dw  = mu*dy / dx; 
Dn  = mu*dx / dy; 
Ds  = mu*dx / dy; 

A = @(F,D)( max(0, (1-0.1 * abs(F/D))^5 ) );

%%compute u_star
for i = 2:imax
    for j = 2:jmax-1
        Fe  = .5*rho*dy*(u(i+1,j)+u(i,j));                                     
        Fw  = .5*rho*dy*(u(i-1,j)+u(i,j)); 
        Fn  = .5*rho*dx*(v(i,j+1)+v(i-1,j+1)); 
        Fs  = .5*rho*dx*(v(i,j)+v(i-1,j));
        
        aE = De * A(Fe,De) + max(-Fe,0);
        aW = Dw * A(Fw,Dw) + max(Fw,0);
        aN = Dn * A(Fn,Dn) + max(-Fn,0);
        aS = Ds * A(Fs,Ds) + max(Fs,0);
        aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
        
        pressure_term = (p(i-1,j)-p(i,j)) * dy;
        
        u_star(i,j) = alpha/aP * ( (aE*u(i+1,j)+aW*u(i-1,j)+aN*u(i,j+1)+aS*u(i,j-1)) + pressure_term ) + (1-alpha)*u(i,j);
        
        d_u(i,j) = alpha * dy / aP;   %refer to Versteeg CFD book
        
    end
end

%%set d_u for top and bottom BCs
%%they will be later used by the pressure correction equation 
%%they should not be zero, or BCs of pressure correction will get messed up
j = 1; %bottom
for i=2:imax
    Fe  = .5*rho*dy*(u(i+1,j)+u(i,j));                                     
    Fw  = .5*rho*dy*(u(i-1,j)+u(i,j)); 
    Fn  = .5*rho*dx*(v(i,j+1)+v(i-1,j+1)); 
    Fs  = 0;
        
    aE = De * A(Fe,De) + max(-Fe,0);
    aW = Dw * A(Fw,Dw) + max(Fw,0);
    aN = Dn * A(Fn,Dn) + max(-Fn,0);
    aS = 0;
    aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
    d_u(i,j) = alpha * dy / aP;
end

j = jmax; %top
for i=2:imax
    Fe  = .5*rho*dy*(u(i+1,j)+u(i,j));                                     
    Fw  = .5*rho*dy*(u(i-1,j)+u(i,j)); 
    Fn  = 0; 
    Fs  = .5*rho*dx*(v(i,j)+v(i-1,j));
        
    aE = De * A(Fe,De) + max(-Fe,0);
    aW = Dw * A(Fw,Dw) + max(Fw,0);
    aN = 0;
    aS = Ds * A(Fs,Ds) + max(Fs,0);
    aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
    d_u(i,j) = alpha * dy / aP;
end


%%Apply BCs
u_star(1,1:jmax) = -u_star(2,1:jmax); %left wall
u_star(imax+1,1:jmax) = -u_star(imax,1:jmax); %right wall
u_star(1:imax+1, 1) = 0.0; %bottom wall
u_star(1:imax+1, jmax) = velocity; %top wall 


return 
end

function [v_star,d_v] = v_momentum(imax,jmax,dx,dy,rho,mu,u,v,p,alpha)

v_star=zeros(imax,jmax+1);
d_v=zeros(imax,jmax+1);

De  = mu*dy / dx;  %convective coefficients
Dw  = mu*dy / dx; 
Dn  = mu*dx / dy; 
Ds  = mu*dx / dy; 

A = @(F,D)( max(0, (1-0.1 * abs(F/D))^5 ) );

%%compute u_star
for i = 2:imax-1
    for j = 2:jmax
        Fe  = .5*rho*dy*(u(i+1,j)+u(i+1,j-1));                                     
        Fw  = .5*rho*dy*(u(i,j)+u(i,j-1)); 
        Fn  = .5*rho*dx*(v(i,j)+v(i,j+1)); 
        Fs  = .5*rho*dx*(v(i,j-1)+v(i,j));
        
        aE = De * A(Fe,De) + max(-Fe,0);
        aW = Dw * A(Fw,Dw) + max(Fw,0);
        aN = Dn * A(Fn,Dn) + max(-Fn,0);
        aS = Ds * A(Fs,Ds) + max(Fs,0);
        aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
        
        pressure_term = (p(i,j-1)-p(i,j)) * dx;
        
        v_star(i,j) = alpha/aP * ( (aE*v(i+1,j)+aW*v(i-1,j)+aN*v(i,j+1)+aS*v(i,j-1)) + pressure_term ) + (1-alpha)*v(i,j);
        
        d_v(i,j) = alpha * dx / aP;   %refer to Versteeg CFD book
        
    end
end

%%set d_v for left and right BCs
%%they will be later used by the pressure correction equation 
%%they should not be zero, or BCs of pressure correction will get messed up
%%Apply BCs
i = 1;  %left BC
for j=2:jmax
    Fe  = .5*rho*dy*(u(i+1,j)+u(i+1,j-1));                                     
    Fw  = 0; 
    Fn  = .5*rho*dx*(v(i,j)+v(i,j+1)); 
    Fs  = .5*rho*dx*(v(i,j-1)+v(i,j));                                                       
        
    aE = De * A(Fe,De) + max(-Fe,0);
    aW = 0;
    aN = Dn * A(Fn,Dn) + max(-Fn,0);
    aS = Ds * A(Fs,Ds) + max(Fs,0);
    aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
    d_v(i,j) = alpha * dx / aP;
end

i = imax;  %right BC
for j=2:jmax
    Fe  = 0;                                     
    Fw  = .5*rho*dy*(u(i,j)+u(i,j-1)); 
    Fn  = .5*rho*dx*(v(i,j)+v(i,j+1)); 
    Fs  = .5*rho*dx*(v(i,j-1)+v(i,j));
        
    aE = 0;
    aW = Dw * A(Fw,Dw) + max(Fw,0);
    aN = Dn * A(Fn,Dn) + max(-Fn,0);
    aS = Ds * A(Fs,Ds) + max(Fs,0);
    aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
    d_v(i,j) = alpha * dx / aP;
end

%%apply BCs
v_star(1,1:jmax+1) = 0.0; %left wall
v_star(imax,1:jmax+1) = 0.0; %right wall
v_star(1:imax, 1) = -v_star(1:imax, 2); %bottom wall
v_star(1:imax, jmax+1) = -v_star(1:imax, jmax); %top wall 


return 
end