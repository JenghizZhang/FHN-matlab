% FHN model in 1d
%     du/dt = (a-u) * (u-1) * u-v + D * Laplacian u;
%     dv/dt = eps * (beta * u - gamma * v - delta);
    
% set some plotting defaults
set(0,'defaultaxesfontsize',16);
set(0,'defaultlinelinewidth',2);

% parameter values
a=0.1; beta=0.5; gamma=1; delta=0; eps=0.01; % excitable
%a=-0.1; beta=0.5; gamma=1; delta=0; eps=0.01; % oscillatory

% numerical and stimulation parameters
dt = 0.5;   % time step size
dx = 1;  % spatial resolution
endtime=800; % simulation duration
nsteps=ceil(endtime/dt);    % calculate number of time steps
diff=0.1; % diffusion coefficient
nx=100;     % domain size

% initial values for state variables
u = zeros(nx,1);       % zeros everywhere
v = zeros(nx,1);    % 0.5 everywhere
% start a propagating wave
u(1,1)=0.8;

% arrays for time and space values (for plotting)
t = dt:dt:endtime;
xx=1:nx;
xx=xx*dx;

% laplacian matrix with Neumann boundary conditions
lapa = speye(nx,nx);
lapb = sparse(2:nx,1:nx-1,1,nx,nx);
A = lapb+lapb'-2*lapa;
A(1,2)=2;
A(nx,nx-1)=2;
ddt_o_dx2=dt*diff/(dx*dx);

%     du/dt = (a-u) * (u-1) * u-v + D * Laplacian u;
%     dv/dt = eps * (beta * u - gamma * v - delta);

% time loop
for ntime=1:nsteps

    % calculate derivatives
    % dots in ".*" allow vectorization; apply everywhere at once (no spatial loops)
    du = (a-u).*(u-1).*u-v; 
    dv = eps*(beta*u-gamma*v-delta);
    
    % calculate coupling term
    xlap=reshape(A*reshape(u,nx,1),nx,1);

    % update variables using forward Euler
    u = u + dt*du + ddt_o_dx2*xlap;
    v = v + dt*dv;

    % plot every 10 time steps
    if(mod(ntime,10)==0)
        figure(1)
        plot(xx,u),title(['t = ',num2str(ntime*dt)]),ylim([-0.4 1.1]),drawnow
    end
	
end

