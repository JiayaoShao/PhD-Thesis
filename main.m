%u and v are heave component; x and y are angles
bu = @(u,v,x,y) v;
bv = @(u,v,x,y) param_h/mass*(0.5*gamma*x.^2-u)-kone*v;
bx = @(u,v,x,y) y;
by = @(u,v,x,y) c/inertia*x.*(x.^2/xv^2-1)+param_h/inertia*gamma*x.*(u-0.5*gamma ...
    *x.^2)-ktwo*y;

% IMPORTANT good approximation figure at c=h=gamma=2
% IMPORTANT action value should be around 0.0470
% all parameters are positive
% h: heave stiffness; positive
param_h = 2;
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

SDE_trajectory
figure(1)
hold on
figure(2)
hold on
figure(4)
hold on

MAM_4D

figure(1)
title('phase diagram of heave')
xlabel('heave configuration u')
ylabel('heave momentum v')
legend('Simulation','MAM')
hold off

figure(2)
title('phase diagram of roll')
xlabel('roll configuration x')
ylabel('roll momentum y')
legend('Simulation','MAM')
hold off

figure(4)
title('phase diagram of roll-heave coupling')
xlabel('roll configuration x')
ylabel('heave configuration u')
legend('Simulation','MAM')
hold off

lagrange = 0.5*((xdot-bx(u,v,x,y)).^2+(ydot-by(u,v,x,y)).^2+ ...
    (udot-bu(u,v,x,y)).^2+(vdot-bv(u,v,x,y)).^2);
action_S = sum(lagrange*dt);

epsilon = linspace(0.0005,0.05,200);

probability = exp(-action_S./epsilon);

figure(3)
semilogy(1./epsilon,probability,'k--')
drawnow


probability_escape = zeros(1,size(epsilon,2));
for i=1:size(epsilon,2)
    eps = epsilon(i);
    SDE_epsilon
    probability_escape(i) = escape/N;
end

figure(3)
hold on
title('Comparison between MAM and stochastic simulation')
xlabel('\epsilon^{-1}')
ylabel('capsize probability')
semilogy(1./epsilon,probability_escape,'-')
legend('MAM','Simulation')



