%u and v are heave component; x and y are angles

bu = @(u,v,x,y) v;
bv = @(u,v,x,y) param_h/mass*(0.5*gamma*x.^2-u)-kone*v;
bx = @(u,v,x,y) y;
by = @(u,v,x,y) c/inertia*x.*(x.^2/xv^2-1)+param_h/inertia*gamma*x.*(u-0.5*gamma ...
    *x.^2)-ktwo*y;


% all parameters are positive
% h: heave stiffness; positive
param_h = 2;
% c: depends on the ocean vehicle
c = 4;
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


steps = 2400;
N = 4000;
dt = 0.01;
total_time = steps*dt;
eps = 0.005;
everyPlot = 400;

u = zeros(N,1);
v = zeros(N,1);
x = zeros(N,1);
y = zeros(N,1);

us = zeros(N,steps);
vs = zeros(N,steps);
xs = zeros(N,steps);
ys = zeros(N,steps);

for step = 1:steps
  u_ = u + dt*bu(u,v,x,y) + sqrt(dt*eps)*randn(size(x));
  v_ = v + dt*bv(u,v,x,y) + sqrt(dt*eps)*randn(size(x));
  x_ = x + dt*bx(u,v,x,y) + sqrt(dt*eps)*randn(size(x));
  y_ = y + dt*by(u,v,x,y) + sqrt(dt*eps)*randn(size(x));
  x = x_; v = v_; u = u_; y = y_;

  x(find(abs(x)>2))=sign(x(find(abs(x)>2)))*2;
  y(find(abs(y)>2))=sign(y(find(abs(y)>2)))*2;

  us(:,step) = u;
  vs(:,step) = v;
  xs(:,step) = x;
  ys(:,step) = y;
  
  %if mod(step,everyPlot)==0
  %  figure(1)
  %  plot(u,v,'x')
  %  xlim([-2,6])
  %  ylim([-2,2])
  %  drawnow
  %end

  % if mod(step,everyPlot)==0
  %  figure(2)
  %  plot(x,y,'x')
  %  xlim([-5,5])
  %  ylim([-5,5])
  %  drawnow
  %end

  %if mod(step,everyPlot)==0
  %  figure(3)
  %  plot(xs',ys')
  %  drawnow
  %end

  %if mod(step,everyPlot)==0
  %  figure(4)
  %  plot(us',vs')
  %  drawnow
  %end
  
end

uplot = us(find(u>1.5),:)';
vplot = vs(find(u>1.5),:)';
xplot = xs(find(u>1.5),:)';
yplot = ys(find(u>1.5),:)';


%%
figure(5)
plot(uplot,vplot)
drawnow
hold on
title('phase diagram of heave')
xlabel('heave configuration')
ylabel('heave momentum')
hold off

figure(6)
plot(xplot,yplot)
drawnow
hold on
title('phase diagram of roll')
xlabel('roll configuration')
ylabel('roll momentum')
hold off

figure(7)
plot(xplot,uplot)
drawnow
hold on
title('heave-roll')
xlabel('roll configuration \theta')
ylabel('heave configuration z')
hold off