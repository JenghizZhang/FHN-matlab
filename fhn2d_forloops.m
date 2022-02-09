% FHN model in 2d
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

ddt_o_dx2=diff*dt/(dx*dx);  % useful combination; avoid computing repeatedly
unew = zeros(nx,ny);        % array to store updated u values

% time loop
for ntime=1:nsteps

    for j=1:ny
        for i=1:nx
            % calculate derivatives at this point
            du = (a-u(i,j))*(u(i,j)-1)*u(i,j)-v(i,j); 
            dv = eps*(beta*u(i,j)-gamma*v(i,j)-delta);
    
            % calculate coupling term
            if(i==1)
                xlap1=2*(u(2,j)-u(1,j));
            elseif(i==nx)
                xlap1=2*(u(nx-1,j)-u(nx,j));
            else
                xlap1=u(i-1,j)-2*u(i,j)+u(i+1,j);
            end
            if(j==1)
                xlap2=2*(u(i,2)-u(i,1));
            elseif(j==ny)
                xlap2=2*(u(i,ny-1)-u(i,ny));
            else
                xlap2=u(i,j-1)-2*u(i,j)+u(i,j+1);
            end
            xlap=xlap1+xlap2;

            % update variables using forward Euler
            unew(i,j) = u(i,j) + dt*du + ddt_o_dx2*xlap;
            v(i,j) = v(i,j) + dt*dv;

        end
    end

    % now that all locations have been updated, overwrite u
    u = unew;
    
    % plot every 20 time steps
    if(mod(ntime,20)==0)
        figure(1)
        pcolor(u),shading interp,caxis([-0.4 1.1]),colorbar,title(['t = ',num2str(ntime*dt)]),daspect([1 1 1]),set(gca,'xtick',[],'ytick',[]),drawnow
    end
    
end

%colormap gray
%subplot(2,1,1)
%pcolor(xx,xx,squeeze(vsave(2,:,:))),shading interp,colorbar
%subplot(2,1,2)
%pcolor(t,xx,hsave'),shading interp,xlabel('Time'),ylabel('Space'), colorbar


