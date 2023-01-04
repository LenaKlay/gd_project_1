%%%%% Numerical simulations for a density-dependant gene drive model       %%%%%%
%%%%% Basis: semi-implicit finite difference scheme                        %%%%%%
%%%%% Goal: export the heatmap associated with the parameters $r$ and $s$  %%%%%%

clear
clf
set (gcf, "papersize", [5,4])
set (gcf, "paperposition", [0, 0, 5, 4])
set(0,'DefaultTextInterpreter','latex')

% Numerical scheme parameters
T = 150; % final time
L = 900; % spatial length
Mt = 2000; % number of time steps
Mx = 3000; % number of space steps
dt = T/Mt; % size of time step
dx = L/Mx; % size of space step
X = [-Mx/2:Mx/2]*dx; % discretized spatial domain
A = spdiags([ones(Mx+1,1) -2*ones(Mx+1,1) ones(Mx+1,1)], [-1, 0, 1], Mx+1, Mx+1); % 1D Laplacian
A(1,1) = A(1,1)+1; % Neumann left-boundary condition
A(end,end) = A(end,end)+1; % Neumann right-boundary condition
Aimp = eye(Mx+1)-(dt/dx^2)*A; % Implicit Laplacian

% Parameters s and r vary in discretized intervals
N = 200;
smin = 0.4;
smax = 0.75;
rmin = 0;
rmax = smax/(1-smax);
ds = (smax-smin)/N;
dr = (rmax-rmin)/N;

% Initialization of the outside (s,r) loop
speed = [];
wb = waitbar(0,'0%…'); % Display a waitbar
for i=[1:(N+1)^2] % Parameter loop starts
is = floor((i-1)/(N+1))+1;
ir = mod((i-1),(N+1))+1;
s = smin+(is-1)*ds;
r = rmin+(ir-1)*dr;

% Initialization of the inside time loop
D = [zeros(Mx/2+1,1);.95*ones(Mx/2,1)]; % n_D at t=0
O = [ones(Mx/2+1,1);.05*ones(Mx/2,1)]; % n_O at t=0
position = find((O<1/2),1);

if mod(i,(N+1))==0
  waitbar(i/(N+1)^2,wb,sprintf('%d%%…',100*i/(N+1)^2)) % Display a waitbar
end

% Time loop starts
for j=[1:Mt]
  d = Aimp\(D + dt*D.*(-1+(r*(1-(D+O))+1).*(1-s).*(2-D./(D+O)))); 
  o = Aimp\(O + dt*O.*(-1+(r*(1-(D+O))+1).*(1-D./(D+O)))); 
  D = d;
  O = o;
  position = [position find((O<1/2),1)]; % Estimate of the position of the front
end % Time loop ends

c =  mean(diff(position)(4*Mt/5:end))*dx/dt;
speed(is,ir) = c;

end % Parameter loop ends

% Save the speed values in a file speed.mat for future use
save speed.mat speed

% Formatting the heatmap
BrownBlueBaseMap = [144 100 44;
             187 120  54;
             225 146  65;
             248 184 139;
             244 218 200;
             241 244 245;
             207 226 240;
             160 190 225;
             109 153 206;
              70  99 174;
              24  79 162]/255;
idx1 = linspace(0,1,11);
idx2 = linspace(0,1,256);
BrownBlueMap = interp1(idx1,BrownBlueBaseMap,idx2);
BlueBrownMap = flipud(BrownBlueMap)
colormap(BlueBrownMap);
imagesc([smin smax],[rmin rmax],speed');
axis xy;
colorbar;
M = max(abs([min(speed'(:)),max(speed'(:))]));
caxis([-M M])
xlabel('$s$','fontsize',10)
ylabel('$r$','fontsize',10)

% Add curves on the heatmap
hold on
contour(smin+[0:N]*ds,rmin+[0:N]*dr,speed',[0 0],'--k','LineWidth',1); % 0-level set of the speeds
plot(smin+[0:N]*ds,(smin+[0:N]*ds)./(1-smin-[0:N]*ds),'k','LineWidth',1); % Eradication curve
sign_curve = max((smin-(1/2)+[0:N]*ds)./(1-smin-[0:N]*ds),(smin+[0:N]*ds).^(3).*(3*(smin+[0:N]*ds)-2)./((2*(smin+[0:N]*ds)-1).^(3)-(smin+[0:N]*ds).^(3).*(3*(smin+[0:N]*ds)-2)));
plot(smin+[0:N]*ds,sign_curve,'k','LineWidth',1); % Boundaries of the unknown sign area
plot([0.5 0.5],[rmin rmax],'k','LineWidth',1); % Monostability-bistability threshold
plot([2/3 2/3],[rmin rmax],'k','LineWidth',1); % Asymptotic 0-level set
hold off

% Export in .pdf + .tex format
print('heatmap','-dpdflatex')

