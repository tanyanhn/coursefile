%% Jake Jordan (2020). 2D Shallow Water Equation, 
%  MATLAB Central File Exchange. Retrieved May 14, 2020. 
%
%% 1. set up finite-difference mesh or grid:
nx = 51; ny = 51;  Lx = 100; Ly = 100;  dx = Lx/(nx-1); dy = Ly/(ny-1);
[x,y] = meshgrid(linspace(0,Lx,nx),linspace(0,Ly,ny));

%% 2. precomput index for reflective boundaries  
eta=zeros(nx,ny);  up=zeros(nx,ny);  vp=zeros(nx,ny);
etam=zeros(nx,ny);  u=zeros(nx,ny);   v=zeros(nx,ny);   
etap=zeros(nx,ny); um=zeros(nx,ny);  vm=zeros(nx,ny);
eta(:,1) = eta(:,2);      up(:,1) = up(:,2);       vp(:,1) = -vp(:,2);
eta(:,ny) = eta(:,ny-1);  up(:,ny) = up(:,ny-1);   vp(:,ny) = -vp(:,ny-1);
eta(1,:) = eta(2,:);      up(1,:) = -up(2,:);      vp(1,:) = vp(2,:);
eta(nx,:) = eta(nx-1,:);  up(nx,:) = -up(nx-1,:);  vp(nx,:) = vp(nx-1,:);

%% 3. setup the initial condition. etam at t=0
xc = 78; yc = 15; k = 5; 
eta = 10*exp(-1/k/k*((x - xc).^2 + (y - yc).^2));

%% 4. main integration
dt = .5; Nsteps = 500;  c = (37:100)'; cyan=[0*c c c]/100; g = 9.81; 
for k = (1:Nsteps)
    t = k*dt;
    
%     for i = 2:nx-1    %% alternative, slower however meaningful version
%         for j = 2:ny-1
%             up(i,j) = um(i,j) - g*(dt/dx)*(eta(i+1,j) - eta(i-1,j));
%             vp(i,j) = vm(i,j) - g*(dt/dy)*(eta(i,j+1) - eta(i,j-1));
%             etap(i,j) = 2*eta(i,j) - etam(i,j) + (((2*dt^2/dx*dy))*...
%               (eta(i+1,j) + eta(i-1,j) + eta(i,j+1) + eta(i,j-1) - 4*eta(i,j)));
%         end
%     end    
    
    up(2:nx-1,2:ny-1) = um(2:nx-1,2:ny-1)...
                    -g*(dt/dx)*(eta(3:nx,2:ny-1) - eta(1:nx-2,2:ny-1));
    vp(2:nx-1,2:ny-1) = vm(2:nx-1,2:ny-1)...
                    -g*(dt/dy)*(eta(2:nx-1,3:ny) - eta(2:nx-1,1:ny-2));
    etap(2:nx-1,2:ny-1) = 2*eta(2:nx-1,2:ny-1)...
                    -etam(2:nx-1,2:ny-1)+(((2*dt^2/dx*dy))*...
                    (eta(3:nx,2:ny-1)+eta(1:nx-2,2:ny-1)...
                    +eta(2:nx-1,3:ny)+eta(2:nx-1,1:ny-2)...
                    -4*eta(2:nx-1,2:ny-1)));
    
    etam = eta; eta = etap;   % prepare for next iteration

    subplot(1,2,1); z=etap; surf(x,y,eta); zlim([-5 5]); 
          colormap(cyan); %shading interp;  zlabel('H +/- \eta'); 
    subplot(1,2,2); surf(x,y,eta); view(0,90);shading flat; drawnow;
end









