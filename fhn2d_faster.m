% FHN model in 2d
% here we speed things up by using vectorization
% we don't need to use space loops at all!
%     du/dt = (a-u) * (u-1) * u-v + D*Laplacian u;
%     dv/dt = eps * (beta * u - gamma * v - delta);

% parameter values
a=0.1; beta=0.5; gamma=1; delta=0; eps=0.01; % excitable
%a=-0.1; beta=0.5; gamma=1; delta=0; eps=0.01; % oscillatory

% numerical and stimulation parameters
dt = 0.5;   % time step size
dx = 1;  % spatial resolution
endtime=1000; % simulation duration
nsteps=ceil(endtime/dt);    % calculate number of time steps
diff=0.1; % diffusion coefficient
nx=100;     % domain size
ny=nx;      % make it square

% initial values for state variables
u = zeros(nx,ny);       % zeros everywhere
v = zeros(nx,ny);    % 0.5 everywhere
% some variations to start a spiral wave
u(1:ceil(nx/2),1:ceil(ny/2)) = 0.8;   % excite part
v(1:ceil(ny/2),:) = 0.1;       % make part refractory

% arrays for time and space values (for plotting)
t = dt:dt:endtime;
xx=1:nx;
xx=xx*dx;

% laplacian matrix with Neumann boundary conditions
I=speye(nx,nx);
E=sparse(2:nx,1:nx-1,1,nx,nx);
D=E+E'-2*I;
D(1,2)=2;
D(nx,nx-1)=2;
A=kron(D,I)+kron(I,D);  % scary (but convenient) kronecker tensor product
ddt_o_dx2=diff*dt/(dx*dx);  % useful combination; avoid computing repeatedly

% time loop
for ntime=1:nsteps

    % calculate derivatives
    % dots in ".*" allow vectorization; apply everywhere at once (no spatial loops)
    du = (a-u).*(u-1).*u-v; 
    dv = eps*(beta*u-gamma*v-delta);
    
    % calculate coupling term
    xlap=reshape(A*reshape(u,nx*nx,1),nx,nx);

    % update variables using forward Euler
    u = u + dt*du + ddt_o_dx2*xlap;
    v = v + dt*dv;

    % plot every 20 time steps
    if(mod(ntime,20)==0)
        figure(1)
        pcolor(u),shading interp,caxis([-0.4 1.1]),colorbar,title(['t = ',num2str(ntime*dt)]),daspect([1 1 1]),set(gca,'xtick',[],'ytick',[]),drawnow
    end
    
end

