c = 1;
xv = 1;
ktwo = 0.2;

V = @(x) 0.5*c*x.^2 - 0.25*c*x.^4/xv^2;
Hamilton = @(x,y) V(x) + 0.5*y.^2;

bx = @(x,y,z) y;
by = @(x,y,z) c*x.^3/xv^2 - c*x - ktwo*y + z;
bz = @(x,y,z) -z;
dbxdx = @(x,y,z) 0; dbxdy = @(x,y,z) 1; dbxdz = @(x,y,z) 0;
dbydx = @(x,y,z) 3*c*x.^2/xv^2 -c; dbydy = @(x,y,z) -ktwo; dbydz = @(x,y,z) 1;
dbzdx = @(x,y,z) 0; dbzdy = @(x,y,z) 0; dbzdz = @(x,y,z) -1;

% define linearisation around sink
A_s = [0,1,0;-c,-ktwo,1;0,0,-1];
sigma = [0;0;1];
covmat = sigma*sigma';
params.covmat = covmat;
Q = lyap(A_s,covmat);
params.Q = inv(Q);

% define linearisation and coordinate transformation around saddle
A_u = [0,1,0;2*c,-ktwo,1;0,0,-1];
[P,D,W] = eig(A_u);
r = P(:,1);
l = W(:,1);
r = r/(r'*l);
params.l = l;
U = [r(1),-l(2)/l(1),-l(3)/l(1);r(2),1,0;r(3),0,1];
params.U = U;
transigma = U\sigma;
saddle = [xv,0,0];
params.saddle = saddle;
q = transigma(1)^2/2/D(1,1);
params.q = q;

% define ellipsoid
Ahat = U\A_u*U;
Ahat(1,1) = -Ahat(1,1);
Ahat(2,1) = -transigma(2)*transigma(1)/q;
Ahat(3,1) = -transigma(3)*transigma(1)/q;
ellipsoid_matrix = lyap(Ahat',eye(3));
params.ellipsoid_matrix = ellipsoid_matrix;
eig(ellipsoid_matrix)
deviation = 0.05;
params.deviation = deviation^2;
% sketch energy level and ellipsoid
[x_,y_] = meshgrid(linspace(-1.25,1.25,100), linspace(-1,1,100));
ellipsoid_origin = inv(U)'*ellipsoid_matrix*inv(U);
[P_ellip,D_ellip] = eig(ellipsoid_origin);
xc = xv; yc = 0; zc = 0;
xr = deviation/sqrt(D_ellip(1,1));
yr = deviation/sqrt(D_ellip(2,2));
zr = deviation/sqrt(D_ellip(3,3));
n=9;
[X,Y,Z] = ellipsoid(0,0,0,xr,yr,zr,n);
Xrot = zeros(n,n);
Yrot = zeros(n,n);
Zrot = zeros(n,n);
for i = 1:n
    for j = 1:n
        rotated = P_ellip*[X(i,j),Y(i,j),Z(i,j)]'; 
        Xrot(i,j) = rotated(1);
        Yrot(i,j) = rotated(2);
        Zrot(i,j) = rotated(3);
    end
end
Xrot = Xrot+xv;


H = @(x,y,z,thx,thy,thz) bx(x,y,z).*thx+by(x,y,z).*thy+bz(x,y,z).*thz+0.5*(thz.^2);
dH_dx = @(x,y,z,thx,thy,thz)dbxdx(x,y,z).*thx+dbydx(x,y,z).*thy+dbzdx(x,y,z).*thz;
dH_dy = @(x,y,z,thx,thy,thz)dbxdy(x,y,z).*thx+dbydy(x,y,z).*thy+dbzdy(x,y,z).*thz;
dH_dz = @(x,y,z,thx,thy,thz)dbxdz(x,y,z).*thx+dbydz(x,y,z).*thy+dbzdz(x,y,z).*thz;
dH_dthx = @(x,y,z) bx(x,y,z);
dH_dthy = @(x,y,z) by(x,y,z);
dH_dthz = @(x,y,z,thz) bz(x,y,z)+thz;

params.bx = bx; params.by = by; params.bz = bz;
params.dH_dx = dH_dx; params.dH_dy = dH_dy; params.dH_dz = dH_dz;
params.dH_dthx = dH_dthx; params.dH_dthy = dH_dthy; params.dH_dthz = dH_dthz;

Nt = 601;
T = 30;
t = linspace(0,T,Nt);
dt = t(2)-t(1);
params.Nt = Nt;
params.dt = dt;

phinull = 0.01*ones(1,3);
thz = 0.1*randn([1,Nt-1]);

params.lambda = 1;
params.beta = 0;
beta_multiplier = zeros(1,501);
actionchange = [];

for outer=1:500
  [gradzold,gradphiold] = gradient(thz,phinull,params);
  dtau = linesearch(thz,phinull,-gradzold,-gradphiold,params);
  thz = thz - dtau*gradzold; 
  phinull = phinull - dtau*gradphiold;
  congradz = gradzold;
  congradphi = gradphiold;

  for iter=1:1e4
	[gradz,gradphi] = gradient(thz,phinull,params);         % compute gradient
    cgbeta = (sum(gradz.*gradz)+sum(gradphi.*gradphi))/(sum(gradzold.*gradzold)+sum(gradphiold.*gradphiold));
    congradz = -gradz + cgbeta*congradz;
    congradphi = -gradphi + cgbeta*congradphi;
	dtau = linesearch(thz,phinull,congradz,congradphi,params); % get best step size
	thz = thz + dtau*congradz;              % gradient descent
    phinull = phinull + dtau*congradphi;
    gradzold = gradz;
    gradphiold = gradphi;

    if sqrt(norm(gradz)^2*dt+norm(gradphi)^2)<5e-5
	  break
    end
  end

  phix = zeros(1,Nt);
  phiy = zeros(1,Nt);
  phiz = zeros(1,Nt);
  phix(1) = phinull(1);
  phiy(1) = phinull(2);
  phiz(1) = phinull(3);

  for s=2:Nt
      phix(s) = phix(s-1) + params.dt*params.dH_dthx(phix(s-1),phiy(s-1),phiz(s-1));
	  phiy(s) = phiy(s-1) + params.dt*params.dH_dthy(phix(s-1),phiy(s-1),phiz(s-1));
      phiz(s) = phiz(s-1) + params.dt*params.dH_dthz(phix(s-1),phiy(s-1),phiz(s-1),thz(s-1));
  end

  pend = [phix(end),phiy(end),phiz(end)];
  pend_ = params.U\(pend-params.saddle)';
  
  figure(1) 
  hold on
  contourf(x_,y_,Hamilton(x_,y_), 13)
  plot(phix, phiy, 'k-x','Color',"#D95319")
  xlabel('x'); ylabel('y')
  title('2D Trajectory')
  hold off
  figure(2)
  plot3(phix,phiy,phiz,'x-');hold on
  surf(Xrot,Yrot,Zrot)
  xlabel('x'); ylabel('y'); zlabel('z')
  title('3D Trajectory')
  hold off

  figure(3)
  plot(t,phiz,'x-')
  
  figure(4)
  plot(t(2:Nt),sigma(3)*thz, 'x-')
  title('Diffusion')
  xlabel('time')
  drawnow


  params.beta = params.beta + params.lambda*(pend_'*params.ellipsoid_matrix*pend_-params.deviation);
  params.lambda = params.lambda*1.1;

  action = 0.5*sum(covmat(3,3)*thz.^2*dt)+0.5*phinull*params.Q*phinull'+pend_(1)^2/2/params.q;
  fprintf('OUTER LOOP %d, NEXT STEP after %d steps: lambda=%g, beta=%g',...
				  outer, iter, params.lambda, params.beta);

  actionchange(end+1) = action; 
  
  beta_multiplier(outer+1) = params.beta;
  diff = beta_multiplier(outer+1)-beta_multiplier(outer);
  if abs(diff) < 1e-6
      break
  end
end


%%
multiplier = Q\phinull';
N=201;
sinkt = linspace(-20,0,N);
deltat = sinkt(2)-sinkt(1);
phit = zeros(3,N);
etasink = zeros(1,N);
for ii = 1:N
    phit(:,ii)=Q*expm(-A_s'*sinkt(ii))*multiplier;
    etasink(:,ii)=sigma'*expm(-A_s'*sinkt(ii))*multiplier;
end
phix_sink = phit(1,:);
phiy_sink = phit(2,:);
phiz_sink = phit(3,:);

figure(4)
hold on
plot(sinkt,etasink,'Color',"#D95319")
hold off
drawnow


%% saddle analytical trajectory
pend_ = U\(pend-saddle)';
multiplier_saddle = pend_(1)/q;
tspan = linspace(0,20,201);
y0 = [pend_(1),pend_(2),pend_(3)];
[saddlet,y] = ode78(@(saddlet,y) ode(saddlet,y,Ahat),tspan,y0);
normalpend = saddle' + U*y';
phix_saddle = normalpend(1,:);
phiy_saddle = normalpend(2,:);
phiz_saddle = normalpend(3,:);

etasaddle = zeros(1,201);
for ii=1:201
    etasaddle(:,ii)=-transigma(1)*exp(Ahat(1,1)*saddlet(ii))*pend_(1)/q;
end

figure(1)
hold on
plot(phix_sink,phiy_sink,'y','LineWidth', 1)
plot(phix_saddle,phiy_saddle,'y','LineWidth', 1)
hold off

figure(2)
hold on
plot3(phix_sink,phiy_sink,phiz_sink,'LineWidth', 1)
plot3(phix_saddle,phiy_saddle,phiz_saddle,'LineWidth', 1)
hold off

saddlet = saddlet+T;
figure(3)
hold on
plot(saddlet,phiz_saddle,'Color',"#D95319")
plot(sinkt,phiz_sink,'Color',"#D95319")
title('filtered noise')
xlabel('time')
ylabel('z')
hold off

figure(4)
hold on
plot(saddlet,etasaddle,'Color',"#D95319")
hold off


%%

function ret=cost(thz,phinull,params)
  phix = zeros(1,params.Nt);
  phiy = zeros(1,params.Nt);
  phiz = zeros(1,params.Nt);
  phix(1) = phinull(1);
  phiy(1) = phinull(2);
  phiz(1) = phinull(3);

  for s=2:params.Nt
      phix(s) = phix(s-1) + params.dt*params.dH_dthx(phix(s-1),phiy(s-1),phiz(s-1));
	  phiy(s) = phiy(s-1) + params.dt*params.dH_dthy(phix(s-1),phiy(s-1),phiz(s-1));
      phiz(s) = phiz(s-1) + params.dt*params.dH_dthz(phix(s-1),phiy(s-1),phiz(s-1),thz(s-1));
  end

  pend = [phix(end),phiy(end),phiz(end)];
  pend_ = params.U\(pend-params.saddle)';

 % cost function (theta + boundary condition and panelty)
  costbyphi = 0.5*sum(params.covmat(3,3)*thz.^2)*params.dt;
  costbyphinull = 0.5*phinull*params.Q*phinull';
  costbyphiend = pend_(1)^2/2/params.q + params.beta*(pend_'*params.ellipsoid_matrix*pend_-params.deviation) ...
      + params.lambda*(pend_'*params.ellipsoid_matrix*pend_-params.deviation)^2;
  ret = costbyphi + costbyphinull + costbyphiend;
end


function [ret1,ret2]=gradient(thz,phinull,params)
  phix = zeros(1,params.Nt);
  phiy = zeros(1,params.Nt);
  phiz = zeros(1,params.Nt);
  phix(1) = phinull(1);
  phiy(1) = phinull(2);
  phiz(1) = phinull(3);

  for s=2:params.Nt
      phix(s) = phix(s-1) + params.dt*params.dH_dthx(phix(s-1),phiy(s-1),phiz(s-1));
	  phiy(s) = phiy(s-1) + params.dt*params.dH_dthy(phix(s-1),phiy(s-1),phiz(s-1));
      phiz(s) = phiz(s-1) + params.dt*params.dH_dthz(phix(s-1),phiy(s-1),phiz(s-1),thz(s-1));
  end
  pend = [phix(end),phiy(end),phiz(end)];

  mux=zeros(1,params.Nt-1);
  muy=zeros(1,params.Nt-1);
  muz=zeros(1,params.Nt-1);

  pend_ = params.U\(pend-params.saddle)';
  pend_ = pend_';
  matrix = inv(params.U)'*params.ellipsoid_matrix*inv(params.U);
  muend = - 2*params.beta*(pend-params.saddle)*matrix - (pend-params.saddle)*params.l*params.l'/params.q ...
      - 2*params.lambda*2*(pend-params.saddle)*matrix*(pend_*params.ellipsoid_matrix*pend_'-params.deviation);

  mux(end) = muend(1);
  muy(end) = muend(2);
  muz(end) = muend(3);

  for s=(params.Nt-1):-1:2  % backward/adjoint equation
	mux(s-1) = mux(s) + params.dt*params.dH_dx(phix(s),phiy(s),phiz(s),mux(s),muy(s),muz(s));
	muy(s-1) = muy(s) + params.dt*params.dH_dy(phix(s),phiy(s),phiz(s),mux(s),muy(s),muz(s));
    muz(s-1) = muz(s) + params.dt*params.dH_dz(phix(s),phiy(s),phiz(s),mux(s),muy(s),muz(s));
  end

  epsilon = zeros(1,3);
  epsilon(1) = - mux(1) - params.dt*params.dH_dx(phix(1),phiy(1),phiz(1),mux(1),muy(1),muz(1));
  epsilon(2) = - muy(1) - params.dt*params.dH_dy(phix(1),phiy(1),phiz(1),mux(1),muy(1),muz(1));
  epsilon(3) = - muz(1) - params.dt*params.dH_dz(phix(1),phiy(1),phiz(1),mux(1),muy(1),muz(1));
  
  ret1 =params.covmat(3,3)*(thz-muz)*params.dt;
  ret2 = phinull*params.Q + epsilon;
end


function ret=linesearch(thz,phinull,p1,p2,params)
  dtau = 1;
  for kk=1:100
   	thz_ = thz + dtau*p1;
    phinull_ = phinull + dtau*p2;
   	if cost(thz_,phinull_,params) < cost(thz,phinull,params)
   	  break
   	else
   	  dtau = dtau*0.5;
   	end
  end
  ret = dtau;
end

function dpdt = ode(t,y,Ahat)
  dpdt = Ahat*y;
end

