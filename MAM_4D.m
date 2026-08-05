%test whether the action has converged
action_Stest = [];

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
eigenvalue = eig(diff_matrix)
[V,D] = eig(diff_matrix);

%initialising
%the tested optimal dt value is 0.08, the greater the Nt, the smaller
%action, but action doesn't change much (e^-4 order) after Nt=600
dt = 0.08;
Nt = 300;
total_time = Nt*dt;
iterations = 2e5;
alpha = 1e-3;
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

  %display the newest numerical value of action
  lagrange = 0.5*((xdot-bx(u,v,x,y)).^2+(ydot-by(u,v,x,y)).^2+ ...
    (udot-bu(u,v,x,y)).^2+(vdot-bv(u,v,x,y)).^2);
  action_Stest = [action_Stest,sum(lagrange*dt)];

end

figure(1)
plot(u,v,'x-','DisplayName','MAM')
drawnow

figure(2)
plot(x,y,'x-','DisplayName','MAM')
drawnow

figure(4)
plot(x,u,'x-','DisplayName','MAM')
drawnow

umam = u;
vmam = v;
xmam = x;
ymam = y;