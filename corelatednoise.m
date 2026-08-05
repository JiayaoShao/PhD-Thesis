frequency_ratio = 2;
c = 2;
param_h = frequency_ratio*c;
gamma = 2;
xv = pi/4;
kone = 0.2;
ktwo = 0.2;
roll_period = 2*pi/sqrt(c);


%construction of covariance matrix a
noisev = 1;
noisey = 1;
rho = 0.6;
a = [noisev^2,rho*noisev*noisey;rho*noisev*noisey,noisey^2];
params.a = a;

%to sketch energy level sets
V = @(x,u) 0.5*c*xv^2*(x.^2/xv^2-0.5*x.^4/xv^4)+0.5*param_h*(u-0.5*gamma*x.^2).^2;
V_phi = @(x) 0.5*c*x.^2 - 0.25*c*x.^4/xv^2;
Hamilton = @(x,y) V_phi(x) + 0.5*y.^2;

bu = @(u,v,x,y) v;
bv = @(u,v,x,y) param_h*(0.5*gamma*x.^2-u)-kone*v;
bx = @(u,v,x,y) y;
by = @(u,v,x,y) c*x.*(x.^2/xv^2-1)+param_h*gamma*x.*(u-0.5*gamma*x.^2)-ktwo*y;

%derivatives of u,v,x and y
dubu = @(u,v,x,y) 0*u; dvbu = @(u,v,x,y) 0*u+1; 
dxbu = @(u,v,x,y) 0*u; dybu = @(u,v,x,y) 0*u;

dubv = @(u,v,x,y) -param_h+0*v; dvbv = @(u,v,x,y) 0*v-kone;
dxbv = @(u,v,x,y) param_h*gamma*x; dybv = @(u,v,x,y) 0*v;

dubx = @(u,v,x,y) 0*x; dvbx = @(u,v,x,y) 0*x;
dxbx = @(u,v,x,y) 0*x; dybx = @(u,v,x,y) 0*x+1;

duby = @(u,v,x,y) param_h*gamma*x; dvby = @(u,v,x,y) 0*y; dyby = @(u,v,x,y) 0*y-ktwo;
dxby = @(u,v,x,y) (c/xv^2-param_h*gamma^2/2)*3*x.^2+(param_h*gamma*u-c);

H = @(pu,pv,px,py,thu,thv,thx,thy) bu(pu,pv,px,py).*thu + bv(pu,pv,px,py).*thv + ...
    bx(pu,pv,px,py).*thx + by(pu,pv,px,py).*thy + 0.5*(a(1,1)*thv.^2 + 2*a(1,2)*thv.*thy + a(2,2)*thy.^2);

dH_dphiu = @(pu,pv,px,py,thu,thv,thx,thy) dubu(pu,pv,px,py).*thu+dubv(pu,pv,px,py).*thv+dubx(pu,pv,px,py).*thx+duby(pu,pv,px,py).*thy;
dH_dphiv = @(pu,pv,px,py,thu,thv,thx,thy) dvbu(pu,pv,px,py).*thu+dvbv(pu,pv,px,py).*thv+dvbx(pu,pv,px,py).*thx+dvby(pu,pv,px,py).*thy;
dH_dphix = @(pu,pv,px,py,thu,thv,thx,thy) dxbu(pu,pv,px,py).*thu+dxbv(pu,pv,px,py).*thv+dxbx(pu,pv,px,py).*thx+dxby(pu,pv,px,py).*thy;
dH_dphiy = @(pu,pv,px,py,thu,thv,thx,thy) dybu(pu,pv,px,py).*thu+dybv(pu,pv,px,py).*thv+dybx(pu,pv,px,py).*thx+dyby(pu,pv,px,py).*thy;

dH_dthetau = @(pu,pv,px,py) bu(pu,pv,px,py);
dH_dthetav = @(pu,pv,px,py,thv,thy) bv(pu,pv,px,py) + a(1,1)*thv + a(1,2)*thy;
dH_dthetax = @(pu,pv,px,py) bx(pu,pv,px,py);
dH_dthetay = @(pu,pv,px,py,thv,thy) by(pu,pv,px,py) + a(2,2)*thy + a(1,2)*thv;

params.dH_dphiu = dH_dphiu;
params.dH_dphiv = dH_dphiv;
params.dH_dphix = dH_dphix;
params.dH_dphiy = dH_dphiy;
params.dH_dthetau = dH_dthetau;
params.dH_dthetav = dH_dthetav;
params.dH_dthetax = dH_dthetax;
params.dH_dthetay = dH_dthetay;

[x_,u_] = meshgrid(linspace(-1,1,100), linspace(-1.5,1.5,100));
[x_,y_] = meshgrid(linspace(-1.5,1.5,100), linspace(-2,2,100));


Nt = 1600;
T = 16;
t = linspace(0,T,Nt);
dt = t(2)-t(1);

params.Nt = Nt;
params.dt = dt;

phi0u = 0;
phi0v = 0;
phi0x = 0;
phi0y = 0;
phiTu = 0.5*gamma*xv^2;
phiTv = 0;
phiTx = xv;
phiTy = 0;

params.phi0u = phi0u;
params.phi0v = phi0v;
params.phi0x = phi0x;
params.phi0y = phi0y;
params.phiTu = phiTu;
params.phiTv = phiTv;
params.phiTx = phiTx;
params.phiTy = phiTy;

pu = linspace(phi0u,phiTu,Nt);
pv = linspace(phi0v,phiTv,Nt);
px = linspace(phi0x,phiTx,Nt);
py = linspace(phi0y,phiTy,Nt);
muu = zeros(1,params.Nt-1);
muv = zeros(1,params.Nt-1);
mux = zeros(1,params.Nt-1);
muy = zeros(1,params.Nt-1);
thu = zeros(1,params.Nt-1);
thv = zeros(1,params.Nt-1);
thx = zeros(1,params.Nt-1);
thy = zeros(1,params.Nt-1);

params.lambda = 0.5;
params.betau = 0;
params.betav = 0;
params.betax = 0;
params.betay = 0;
beta_multiplier = zeros(2,100);


for outer=1:100
  [gradvold,gradyold] = gradient(thv,thy,params);
  dtau = linesearch(thv,thy,-gradvold,-gradyold,params);
  thv = thv - dtau*gradvold;
  thy = thy - dtau*gradyold;
  congradv = gradvold;
  congrady = gradyold;
  for iter=1:1e4
	%if mod(iter,50)==0
	%  % plotting
	%  for s=2:Nt
    %   pu(s) = pu(s-1) + params.dt*params.dH_dthetau(pu(s-1),pv(s-1),px(s-1),py(s-1));
    %   pv(s) = pv(s-1) + params.dt*params.dH_dthetav(pu(s-1),pv(s-1),px(s-1),py(s-1),thv(s-1),thy(s-1));
	%   px(s) = px(s-1) + params.dt*params.dH_dthetax(pu(s-1),pv(s-1),px(s-1),py(s-1));
	%   py(s) = py(s-1) + params.dt*params.dH_dthetay(pu(s-1),pv(s-1),px(s-1),py(s-1),thv(s-1),thy(s-1));
    %  end
    %  figure(1)
	%  contourf(x_,u_,V(x_,u_), 13); hold on;
    %  plot(pu,pv,'r-x'); 
    %  plot(px,pu,'g-x'); hold off

    % figure(2)
    %  contourf(x_,y_,Hamilton(x_,y_), 13); hold on;
    %  plot(px, py, 'k-x'); hold off
 
	%  drawnow
	%  action = 0.5*sum((params.a(1,1)*thv.^2 + 2*params.a(1,2)*thv.*thy + params.a(2,2)*thy.^2)*dt);
	%  display(sprintf('action=%g, normv=%g, normy=%g, dtau=%g', action, norm(gradv),norm(grady),dtau))
    %end

    [gradv,grady] = gradient(thv,thy,params);
    cgbeta = (sum(gradv.*gradv)+sum(grady.*grady))/(sum(gradvold.*gradvold)+sum(gradyold.*gradyold));

    congradv = -gradv + cgbeta*congradv;
    congrady = -grady + cgbeta*congrady;

    dtau = linesearch(thv,thy,congradv,congrady,params);  % get best step size
	thv = thv + dtau*congradv;
    thy = thy + dtau*congrady;              % gradient descent

    gradvold = gradv;
    gradyold = grady;

    threshold = sqrt(norm(gradv)^2+norm(grady)^2);
    if threshold < 5e-5
       break
    end
  end

  params.betau = params.betau + params.lambda*(pu(end)-phiTu);
  params.betav = params.betav + params.lambda*(pv(end)-phiTv);
  params.betax = params.betax + params.lambda*(px(end)-phiTx);
  params.betay = params.betay + params.lambda*(py(end)-phiTy);
  params.lambda = params.lambda*1.1;
  display(sprintf('OUTER LOOP %d, NEXT STEP after %d steps: lambda=%g, betax=%g, betay=%g',...
				  outer, iter, params.lambda, params.betax, params.betay));
  beta_multiplier(:,outer)=[params.betau;params.betax];
  
  for s=2:Nt
      pu(s) = pu(s-1) + params.dt*params.dH_dthetau(pu(s-1),pv(s-1),px(s-1),py(s-1));
      pv(s) = pv(s-1) + params.dt*params.dH_dthetav(pu(s-1),pv(s-1),px(s-1),py(s-1),thv(s-1),thy(s-1));
      px(s) = px(s-1) + params.dt*params.dH_dthetax(pu(s-1),pv(s-1),px(s-1),py(s-1));
      py(s) = py(s-1) + params.dt*params.dH_dthetay(pu(s-1),pv(s-1),px(s-1),py(s-1),thv(s-1),thy(s-1));
  end

  %%
  figure(1)
  contourf(x_,u_,V(x_,u_), 13); hold on; 
  plot(px,pu,'g-x'); 
  title('heave roll')
  xlabel('roll configuration x')
  ylabel('heave configuration u')
  hold off

  figure(2)
  contourf(x_,y_,Hamilton(x_,y_), 13); hold on;
  plot(px, py, 'k-x')
  title('phase diagram of roll')
  xlabel('roll configuration x')
  ylabel('roll momentum y')
  hold off

  drawnow
  action = 0.5*sum((params.a(1,1)*thv.^2 + 2*params.a(1,2)*thv.*thy + params.a(2,2)*thy.^2)*dt);
  display(sprintf('action=%g, normv=%g, normy=%g, dtau=%g', action, norm(gradv),norm(grady),dtau))

  theta = [thv;thy];
  force = a*theta;
  forcev = force(1,:);
  forcey = force(2,:);
  figure(3)
  plot(forcev,'r'); hold on;
  title('external force')
  plot(forcey,'g'); 
  legend('noise on heave','noise on roll')
  hold off;
end


function ret=cost(thv,thy,params)
  pu = zeros(1,params.Nt);
  pv = zeros(1,params.Nt);
  px = zeros(1,params.Nt);
  py = zeros(1,params.Nt);
  pu(1) = params.phi0u;
  pv(1) = params.phi0v;
  px(1) = params.phi0x;
  py(1) = params.phi0y;

  %integrate the phi from time 0 to time T
  for s=2:params.Nt
    pu(s) = pu(s-1) + params.dt*params.dH_dthetau(pu(s-1),pv(s-1),px(s-1),py(s-1));
    pv(s) = pv(s-1) + params.dt*params.dH_dthetav(pu(s-1),pv(s-1),px(s-1),py(s-1),thv(s-1),thy(s-1));
	px(s) = px(s-1) + params.dt*params.dH_dthetax(pu(s-1),pv(s-1),px(s-1),py(s-1));
	py(s) = py(s-1) + params.dt*params.dH_dthetay(pu(s-1),pv(s-1),px(s-1),py(s-1),thv(s-1),thy(s-1));
  end
  
  %cost function (theta + boundary condition and panelty)
  ret = 0.5*sum((params.a(1,1)*thv.^2 + 2*params.a(1,2)*thv.*thy + params.a(2,2)*thy.^2)*params.dt) ...
      + params.lambda*((pu(params.Nt)-params.phiTu).^2+(pv(params.Nt)-params.phiTv).^2+...
      (px(params.Nt)-params.phiTx).^2+(py(params.Nt)-params.phiTy).^2) ...
      + params.betau*(pu(params.Nt)-params.phiTu) + params.betav*(pv(params.Nt)-params.phiTv) ...
      + params.betax*(px(params.Nt)-params.phiTx) + params.betay*(py(params.Nt)-params.phiTy);
end

function [ret1,ret2]=gradient(thv,thy,params)
  pu = zeros(1,params.Nt);
  pv = zeros(1,params.Nt);
  px = zeros(1,params.Nt);
  py = zeros(1,params.Nt);
  pu(1) = params.phi0u;
  pv(1) = params.phi0v;
  px(1) = params.phi0x;
  py(1) = params.phi0y;
  for s=2:params.Nt  % forward equation
    pu(s) = pu(s-1) + params.dt*params.dH_dthetau(pu(s-1),pv(s-1),px(s-1),py(s-1));
    pv(s) = pv(s-1) + params.dt*params.dH_dthetav(pu(s-1),pv(s-1),px(s-1),py(s-1),thv(s-1),thy(s-1));
	px(s) = px(s-1) + params.dt*params.dH_dthetax(pu(s-1),pv(s-1),px(s-1),py(s-1));
	py(s) = py(s-1) + params.dt*params.dH_dthetay(pu(s-1),pv(s-1),px(s-1),py(s-1),thv(s-1),thy(s-1));
  end
  muu=zeros(1,params.Nt-1);
  muv=zeros(1,params.Nt-1);
  mux=zeros(1,params.Nt-1);
  muy=zeros(1,params.Nt-1);
  % mu = partial boundary condition partial phi
  muu(params.Nt-1) = -(2*params.lambda*(pu(params.Nt)-params.phiTu)+params.betau);
  muv(params.Nt-1) = -(2*params.lambda*(pv(params.Nt)-params.phiTv)+params.betav);
  mux(params.Nt-1) = -(2*params.lambda*(px(params.Nt)-params.phiTx)+params.betax);
  muy(params.Nt-1) = -(2*params.lambda*(py(params.Nt)-params.phiTy)+params.betay);
  for s=params.Nt-1:-1:2  % backward/adjoint equation
    muu(s-1) = muu(s) + params.dt*params.dH_dphiu(pu(s),pv(s),px(s),py(s),muu(s),muv(s),mux(s),muy(s));
    muv(s-1) = muv(s) + params.dt*params.dH_dphiv(pu(s),pv(s),px(s),py(s),muu(s),muv(s),mux(s),muy(s));
	mux(s-1) = mux(s) + params.dt*params.dH_dphix(pu(s),pv(s),px(s),py(s),muu(s),muv(s),mux(s),muy(s));
	muy(s-1) = muy(s) + params.dt*params.dH_dphiy(pu(s),pv(s),px(s),py(s),muu(s),muv(s),mux(s),muy(s));
  end

  update = [thv-muv;thy-muy];
  grad = params.a*update;
  ret1 = grad(1,:)*params.dt;
  ret2 = grad(2,:)*params.dt;
end

%deciding for each iteration what is the optimal step size (alpla)
function ret=linesearch(thv,thy,p1,p2,params)
  dtau = 1;
  for kk=1:100
    thv_ = thv + dtau*p1;
   	thy_ = thy + dtau*p2;
   	if cost(thv_,thy_,params) < cost(thv,thy,params)
   	  break
    else
   	  dtau = dtau*0.5;
   	end
  end
  ret = dtau;
end
