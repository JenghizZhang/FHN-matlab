% FHN model in 1d
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

% initial values for state variables
u = zeros(nx,1);       % zeros everywhere
v = zeros(nx,1);    % 0.5 everywhere
% start a propagating wave
u(1,1)=0.8;

% arrays for time and space values (for plotting)
t = dt:dt:endtime;
xx=1:nx;
xx=xx*dx;

ddt_o_dx2=diff*dt/(dx*dx);  % useful combination; avoid computing repeatedly
unew = zeros(nx,1);        % array to store updated u values

% time loop
for ntime=1:nsteps
    
    for i=1:nx
        % calculate derivatives at this point
        du = (a-u(i))*(u(i)-1)*u(i)-v(i);
        dv = eps*(beta*u(i)-gamma*v(i)-delta);
        
        % calculate coupling term
        if(i==1)
            xlap=2*(u(2)-u(1));
        elseif(i==nx)
            xlap=2*(u(nx-1)-u(nx));
        else
            xlap=u(i-1)-2*u(i)+u(i+1);
        end
        
        % update variables using forward Euler
        unew(i) = u(i) + dt*du + ddt_o_dx2*xlap;
        v(i) = v(i) + dt*dv;
        
    end
    
    % now that all locations have been updated, overwrite u
    u = unew;
    
    % plot every 10 time steps
    if(mod(ntime,10)==0)
        figure(1)
        plot(xx,u),title(['t = ',num2str(ntime*dt)]),ylim([-0.4 1.1]),drawnow
    end
    
end


