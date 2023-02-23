%%%%% Numerical simulations for a density-dependant gene drive model       %%%%%%
%%%%% Basis: semi-implicit finite difference scheme                        %%%%%%
%%%%% Goal: export time snapshots of the densities in space                %%%%%%

clear
clf
set (gcf, "papersize", [5,4])
set (gcf, "paperposition", [0, 0, 5, 4])
set(0,'DefaultTextInterpreter','latex')

% Biological parameters
s = 0.7; % fitness cost
r = 0.5; % Malthusian growth rate

% Numerical scheme parameters
T = 50; % final time
L = 150; % spatial length
Mt = 1000; % number of time steps
Mx = 1000; % number of space steps
dt = T/Mt; % size of time step
dx = L/Mx; % size of space step
X = [-Mx:0]*dx; % discretized spatial domain
A = spdiags([ones(Mx+1,1) -2*ones(Mx+1,1) ones(Mx+1,1)], [-1, 0, 1], Mx+1, Mx+1); % 1D Laplacian
A(1,1) = A(1,1)+1; % Neumann left-boundary condition
A(end,end) = A(end,end)+1; % Neumann right-boundary condition
Aimp = eye(Mx+1)-(dt/dx^2)*A; % Implicit Laplacian

% Initialisation of the time loop
D = [zeros(Mx/2+1,1);.95*ones(Mx/2,1)]; % n_D at t=0
O = [ones(Mx/2+1,1);.05*ones(Mx/2,1)]; % n_O at t=0
plot(X,D,'-b','LineWidth',2,X,O,'-g','LineWidth',2) % plots at t=0
axis([-Mx*dx 0 -0.1 1.1])
xlabel('$x$','fontsize',10)
ylabel('$(n_D,n_O)$','fontsize',10)
grid on
print('evolution_to_KPP_0','-dpdflatex') % export in .pdf + .tex format

% Time loop
for i=[1:Mt]
  d = Aimp\(D + dt*D.*(-1+(r*(1-(D+O))+1).*(1-s).*(2-D./(D+O)))); 
  o = Aimp\(O + dt*O.*(-1+(r*(1-(D+O))+1).*(1-D./(D+O)))); 
  D = d;
  O = o;
  if mod(i,Mt/10)==0 % At time t=0.1*T,0.2*T,…
    % Plot the densities in space at time t
    plot(X,D,'-b','LineWidth',2,X,O,'-g','LineWidth',2)
    axis([-Mx*dx 0 -0.1 1.1])
    xlabel('$x$','fontsize',10)
    ylabel('$(n_D,n_O)$','fontsize',10)
    grid on
    drawnow;
    % Export in .pdf + .tex format
    print(['evolution_to_KPP_',num2str(i*10/Mt)],'-dpdflatex')
  endif
end
