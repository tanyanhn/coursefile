m = 20; [x,y]= meshgrid(-1: (2/(m-1)): 1);
U = exp(-5*(x.^2+y.^2));
clf; shg; h = surf(x,y,U); axis off; ax = axis; 
c = (37:100)'; cyan=[0*c c c]/100; colormap(cyan);
omega = 0.8;
n = [2:m 1]; e = n; s = [m 1:m-1]; w = s;
T = 500;  
for t = 1:T
    
    U = (1-omega)*U + omega*(U(n,:)+U(:,e)+U(s,:)+U(:,w))*0.25; 
    
    set(h,'zdata',U); title(['t = ', num2str(t)]); 
    axis(ax); drawnow; pause(0.2);
end