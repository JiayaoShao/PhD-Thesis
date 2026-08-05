Nt = 200;
dt = 0.08;
iteration = 1e7;

% initial adaptive epsilon
epsilon = 1e-7;
epsilon_max = 1e-5;
epsilon_min = 1e-14;

% line search parameters
maxLineSearch = 25;
shrinkFactor = 0.5;
growFactor = 1.1;

everyEL = 100000;

% c: depends on the ocean vehicle
c = 4;
h = 2;

% xv: angle of vanishing stability
xv = pi/4;
gamma = 2;
kone = 0.2;
ktwo = 0.2;

x_ = linspace(0,xv,Nt);
z_ = linspace(0,0.5*gamma*xv^2,Nt);
x = x_';
z = z_';

[x,z] = impose_bc(x,z,xv,gamma);

Gx = @(x,xdot,xdotdot,z,zdot,zdotdot) ...
      (3*h*gamma^2*x.^2-6*c*x.^2/xv^2+2*c-2*h*gamma*z-ktwo^2).*xdotdot ...
      -2*h*gamma*x.*zdotdot ...
      - (kone+ktwo)*h*gamma*x.*zdot ...
      - 2*h*gamma*zdot.*xdot ...
      + (h*gamma^2/2-c/xv^2)*6*x.*xdot.^2 ...
      - gamma*h^2*x.*(z-gamma*x.^2/2) ...
      + ((h*gamma^2/2-c/xv^2)*3*x.^2+c-h*gamma*z) ...
        .* ((h*gamma^2/2-c/xv^2)*x.^3+(c-h*gamma*z).*x);

Gz = @(x,xdot,xdotdot,z,zdot,zdotdot) ...
      (2*h-kone^2)*zdotdot ...
      - 2*h*gamma*x.*xdotdot ...
      + x.*xdot*h*gamma*(kone-ktwo) ...
      - gamma*h*xdot.^2 ...
      + h^2*(z-gamma*x.^2/2).*(1+gamma^2*x.^2) ...
      + c*h*gamma*x.^2.*(x.^2/xv^2-1);

% construct the fourth order operator matrix
a1 = [3,-14,26,-24,11,-2,zeros(1,Nt-6)];
a2 = [2,-9,16,-14,6,-1,zeros(1,Nt-6)];
abend = [zeros(1,Nt-6),-1,6,-14,16,-9,2];
aend = [zeros(1,Nt-6),-2,11,-24,26,-14,3];

D = [a1;a2];

for ii = 1:Nt-4
    a = zeros(1,Nt);
    a(ii:ii+4) = [1,-4,6,-4,1];
    D = [D;a];
end

D = [D;abend;aend];

I = eye(Nt);

% initial cost
cost_current = compute_cost(x,z,dt,h,gamma,c,xv,kone,ktwo);

best_cost = cost_current;
best_x = x;
best_z = z;

fprintf('Initial cost = %.12e\n', cost_current);

for iter = 1:iteration

    % compute derivatives from current x,z
    [xdot,xdotdot,zdot,zdotdot] = compute_derivatives(x,z,dt);

    % compute gradient-like residuals using current x,z only
    Gx_old = Gx(x,xdot,xdotdot,z,zdot,zdotdot);
    Gz_old = Gz(x,xdot,xdotdot,z,zdot,zdotdot);

    % backtracking line search
    eps_try = epsilon;
    accepted = false;

    for ls = 1:maxLineSearch

        M = I + eps_try*D/dt^4;

        % synchronous update
        x_try = M \ (x - eps_try*Gx_old);
        z_try = M \ (z - eps_try*Gz_old);

        % impose boundary conditions after proposal
        [x_try,z_try] = impose_bc(x_try,z_try,xv,gamma);

        % reject immediately if NaN/Inf appears
        if any(~isfinite(x_try)) || any(~isfinite(z_try))
            eps_try = shrinkFactor * eps_try;
            continue;
        end

        cost_try = compute_cost(x_try,z_try,dt,h,gamma,c,xv,kone,ktwo);

        if isfinite(cost_try) && cost_try <= cost_current
            accepted = true;
            break;
        else
            eps_try = shrinkFactor * eps_try;
        end

        if eps_try < epsilon_min
            break;
        end
    end

    if accepted
        x = x_try;
        z = z_try;
        cost_current = cost_try;

        % mildly increase epsilon after successful step
        epsilon = min(growFactor*eps_try, epsilon_max);

        if cost_current < best_cost
            best_cost = cost_current;
            best_x = x;
            best_z = z;
        end

    else
        fprintf('Line search failed at iter = %d\n', iter);
        fprintf('Current cost = %.12e, best cost = %.12e\n', cost_current, best_cost);
        fprintf('Last epsilon tried = %.3e\n', eps_try);
        break;
    end

    if mod(iter,everyEL)==0
        [xdot,xdotdot,zdot,zdotdot] = compute_derivatives(x,z,dt);

        fprintf('iter = %d, cost = %.12e, best = %.12e, epsilon = %.3e\n', ...
                iter, cost_current, best_cost, epsilon);

        figure(1)
        plot(x, xdot, 'x-')
        title('phase diagram of roll')
        xlabel('configuration of roll')
        ylabel('momentum of roll')
        drawnow

        figure(2)
        plot(z, zdot, 'x-')
        title('phase diagram of heave')
        xlabel('configuration of heave')
        ylabel('momentum of heave')
        drawnow

        figure(3)
        plot(x,z,'x-')
        title('roll-heave coupling')
        xlabel('configuration of roll')
        ylabel('configuration of heave')
        drawnow
    end
end

% use the best trajectory found
x = best_x;
z = best_z;

[xdot,xdotdot,zdot,zdotdot] = compute_derivatives(x,z,dt);

eta = zdotdot + h*(z-0.5*gamma*x.^2) + kone*zdot;

xi = xdotdot ...
     - c*x.*(x.^2/xv^2-1) ...
     - h*gamma*x.*(z-0.5*gamma*x.^2) ...
     + ktwo*xdot;

cost = 0.5*dt*sum(eta.^2+xi.^2);

fprintf('\nFinal reported cost = %.12e\n', cost);
fprintf('Best line-search cost = %.12e\n', best_cost);

figure(1)
hold on
plot(x,xdot,'x-')
title('phase diagram of roll')
xlabel('configuration of roll')
ylabel('momentum of roll')
hold off

figure(2)
hold on
plot(z,zdot,'x-')
title('phase diagram of heave')
xlabel('configuration of heave')
ylabel('momentum of heave')
hold off

figure(3)
hold on
title('roll-heave coupling')
xlabel('configuration of roll')
ylabel('configuration of heave')
plot(x,z,'x-')
drawnow
hold off

figure(4)
local_action = 0.5*(eta.^2 + xi.^2);
plot((0:Nt-1)'*dt, local_action, 'x-')
xlabel('t')
ylabel('local action density')
title('local action density')

figure(5)
cum_action = dt*cumsum(local_action);
plot((0:Nt-1)'*dt, cum_action, 'x-')
xlabel('t')
ylabel('cumulative action')
title('cumulative action')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,z] = impose_bc(x,z,xv,gamma)

    x(1) = 0;
    x(end) = xv;

    x(2) = 3*x(1)/4 + x(3)/4;
    x(end-1) = 3*x(end)/4 + x(end-2)/4;

    z(1) = 0;
    z(end) = gamma*xv^2/2;

    z(2) = 3*z(1)/4 + z(3)/4;
    z(end-1) = 3*z(end)/4 + z(end-2)/4;

end


function [xdot,xdotdot,zdot,zdotdot] = compute_derivatives(x,z,dt)

    xdot = [0; ...
            (x(3:end)-x(1:end-2))/(2*dt); ...
            0];

    xdotdot = [(-x(4)+4*x(3)-5*x(2)+2*x(1))/dt^2; ...
               (x(3:end)-2*x(2:end-1)+x(1:end-2))/dt^2; ...
               (2*x(end)-5*x(end-1)+4*x(end-2)-x(end-3))/dt^2];

    zdot = [0; ...
            (z(3:end)-z(1:end-2))/(2*dt); ...
            0];

    zdotdot = [(-z(4)+4*z(3)-5*z(2)+2*z(1))/dt^2; ...
               (z(3:end)-2*z(2:end-1)+z(1:end-2))/dt^2; ...
               (2*z(end)-5*z(end-1)+4*z(end-2)-z(end-3))/dt^2];

end


function cost = compute_cost(x,z,dt,h,gamma,c,xv,kone,ktwo)

    [xdot,xdotdot,zdot,zdotdot] = compute_derivatives(x,z,dt);

    eta = zdotdot + h*(z-0.5*gamma*x.^2) + kone*zdot;

    xi = xdotdot ...
         - c*x.*(x.^2/xv^2-1) ...
         - h*gamma*x.*(z-0.5*gamma*x.^2) ...
         + ktwo*xdot;

    cost = 0.5*dt*sum(eta.^2 + xi.^2);

end