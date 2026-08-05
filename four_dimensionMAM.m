%u and v are heave component; x and y are angles

bu = @(u,v,x,y) v;
bv = @(u,v,x,y) param_h/mass*(0.5*gamma*x.^2-u)-kone*v;
bx = @(u,v,x,y) y;
by = @(u,v,x,y) c/inertia*x.*(x.^2/xv^2-1)+param_h/inertia*gamma*x.*(u-0.5*gamma ...
    *x.^2)-ktwo*y;


% all parameters are positive
% h: heave stiffness; positive
param_h = 4;
% c: depends on the ocean vehicle
c = 2;
% m and i are corresponding mass or moment of inertia, could be rescaled so
% assume to be 1
mass = 1;
inertia = 1;
% xv: angle of vanishing stability
xv = pi/4;
% the symmetric static variation z is equal to 0.5*gamma*x^2;
gamma = 2;
% the ks are damping parameters
kone = 0.2;
ktwo = 0.2;

%derivatives of u,v,x and y
dubu = @(u,v,x,y) 0*u; dvbu = @(u,v,x,y) 0*u+1; 
dxbu = @(u,v,x,y) 0*u; dybu = @(u,v,x,y) 0*u;

dubv = @(u,v,x,y) -param_h/mass+0*v; dvbv = @(u,v,x,y) 0*v-kone;
dxbv = @(u,v,x,y) param_h*gamma*x/mass; dybv = @(u,v,x,y) 0*v;

dubx = @(u,v,x,y) 0*x; dvbx = @(u,v,x,y) 0*x;
dxbx = @(u,v,x,y) 0*x; dybx = @(u,v,x,y) 0*x+1;


duby = @(u,v,x,y) param_h*gamma*x/inertia; dvby = @(u,v,x,y) 0*y; dyby = @(u,v,x,y) 0*y-ktwo;
dxby = @(u,v,x,y) (c/inertia/xv^2-param_h*gamma^2/2/inertia)*3*x.^2+(param_h*gamma*u/inertia-c/inertia);

%checking stability and eigenvector at non-trivial point
x_fixed = xv; u_fixed = gamma*x_fixed^2/2;
diff_matrix = [0,1,0,0; -param_h/mass,-kone,param_h*gamma*x_fixed/mass,0; 0,0,0,1;...
    param_h*gamma*x_fixed/inertia,0,(c/inertia/xv^2-param_h*gamma^2/2/inertia)*3*x_fixed^2+(param_h*gamma*u_fixed/inertia-c/inertia),-ktwo];
eigenvalue = eig(diff_matrix);
[V,D] = eig(diff_matrix);

%initialising
dt = 0.08;
Nt = 300;
total_time = Nt*dt;
iterations = 1e7;
alpha = 1e-5;
everyPlot = 10000;
v = zeros(1,Nt);
x = linspace(0,xv,Nt);
u = linspace(0,gamma*x(end)^2/2,Nt);
y = zeros(1,Nt);

for iter = 1:iterations
  udotdot = [0, (u(3:end) - 2*u(2:end-1) + u(1:end-2))/dt^2, 0];
  udot = [0, (u(3:end)-u(1:end-2))/(2*dt), 0];
  vdotdot = [0, (v(3:end) - 2*v(2:end-1) + v(1:end-2))/dt^2, 0];
  vdot = [0, (v(3:end)-v(1:end-2))/(2*dt), 0];
  xdotdot = [0, (x(3:end) - 2*x(2:end-1) + x(1:end-2))/dt^2, 0];
  xdot = [0, (x(3:end)-x(1:end-2))/(2*dt), 0];
  ydotdot = [0, (y(3:end) - 2*y(2:end-1) + y(1:end-2))/dt^2, 0];
  ydot = [0, (y(3:end)-y(1:end-2))/(2*dt), 0];


  dSu = - udotdot - vdot.*(dubv(u,v,x,y)-dvbu(u,v,x,y)) - ydot.*(duby(u,v,x,y)) ...
      + bv(u,v,x,y).*dubv(u,v,x,y) + by(u,v,x,y).*duby(u,v,x,y);
  dSv = - vdotdot - udot.*(dvbu(u,v,x,y)-dubv(u,v,x,y)) - xdot.*(-dxbv(u,v,x,y)) ...
      + bu(u,v,x,y).*dvbu(u,v,x,y) + bv(u,v,x,y).*dvbv(u,v,x,y);
  dSx = - xdotdot - vdot.*(dxbv(u,v,x,y)) - ydot.*(dxby(u,v,x,y)-dybx(u,v,x,y)) ...
      + bv(u,v,x,y).*dxbv(u,v,x,y) + by(u,v,x,y).*dxby(u,v,x,y);
  dSy = - ydotdot - udot.*(-duby(u,v,x,y)) - xdot.*(dybx(u,v,x,y)-dxby(u,v,x,y)) ...
      + bx(u,v,x,y).*dybx(u,v,x,y) + by(u,v,x,y).*dyby(u,v,x,y);
  
  
  u(1)=0; u(end)=gamma*x(end)^2/2;
  v(1)=0; v(end)=0;
  x(1)=0; x(end)=xv;
  y(1)=0; y(end)=0;
    
  u = u - alpha*dSu;
  v = v - alpha*dSv;
  x = x - alpha*dSx;
  y = y - alpha*dSy;


%   if mod(iter,everyPlot)==0
%     figure(1)
% 	plot(x,u,'x-')
% 	drawnow
%   end
% 
%  
%   if mod(iter,everyPlot)==0
%     figure(2)
% 	plot(x, y,'x-')
% 	drawnow
%   end
% 
%   if mod(iter,everyPlot)==0
% 	 figure(3)
%      plot(u,v,'x-')
%  	 drawnow
%   end

end

%%
figure(1)
plot(x,u,'x-')
hold on
title('roll-heave')
xlabel('roll configuration x')
ylabel('heave configuration u')
hold off

figure(2)
plot(x, y,'x-')
hold on
title('phase diagram of roll')
xlabel('roll configuration x')
ylabel('roll momentum y')
hold off

figure(3)
hold on
plot(u,v,'x-')
title('phase diagram of heave')
xlabel('heave configuration u')
ylabel('heave momentum v')
hold off


action = 0;
for i=1:Nt-1
    actionu = (u(i+1)-u(i))/dt - bu(u(i),v(i),x(i),y(i));
    actionv = (v(i+1)-v(i))/dt - bv(u(i),v(i),x(i),y(i));
    actionx = (x(i+1)-x(i))/dt - bx(u(i),v(i),x(i),y(i));
    actiony = (y(i+1)-y(i))/dt - by(u(i),v(i),x(i),y(i));
    action = action+actionu^2+actionv^2+actionx^2+actiony^2;
end
action = action*dt/2
