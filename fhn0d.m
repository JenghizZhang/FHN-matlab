% FHN model in 0d
%     du/dt = (a-u) * (u-1) * u-v;
%     dv/dt = eps * (beta * u - gamma * v - delta);
    
% set some plotting defaults
set(0,'defaultaxesfontsize',16);
set(0,'defaultlinelinewidth',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameter values
a=0.1; beta=0.5; gamma=1; delta=0; eps=0.01; % excitable
% a=-0.1; beta=0.5; gamma=1; delta=0; eps=0.01; % oscillatory

% initial values
u = 0.7;
v = 0;

% numerical parameters
dt = 0.1;           % time step
endtime = 400;      % duration of simulation
nsteps=ceil(endtime/dt);    % calculate # of time steps

% set up arrays/vectors for saving results
usave = zeros(nsteps+1,1);  % size nsteps+1, initialized to 0
vsave = zeros(nsteps+1,1);  % size nsteps+1, initialized to 0

% set first entry for both u and v to initial values
usave(1,1) = u;
vsave(1,1) = v;

% set up a vector/array of time points to use for plotting
t = 0:dt:endtime;

% loop through all time steps (iterate!)
for ntime=1:nsteps

    % calculate derivatives
    du=(a-u)*(u-1)*u-v;
    dv=eps*(beta*u-gamma*v-delta);
    
    % update variables using forward Euler
    u = u + dt*du;
    v = v + dt*dv;
    
    % save updated variables to array for plotting
    usave(ntime+1,1) = u;
    vsave(ntime+1,1) = v;
end

%     du/dt = (a-u) * (u-1) * u-v;
%     dv/dt = eps * (beta * u - gamma * v - delta);

% plotting
% phase space on top
subplot(2,1,1)
uu=-0.5:.01:1.2;
hold off
% plot trajectory calculated
plot(usave,vsave,'r','linewidth',2)
hold on
% plot nullclines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(uu,uu.*(uu-1).*(a-uu),'k','linewidth',1)
plot(uu,beta*uu/gamma-delta/gamma,'k','linewidth',1)
hold off
xlabel('u'),ylabel('v')

% plot u and v as functions of time
subplot(2,1,2)
plot(t,usave,'k',t,vsave,'r','linewidth',2)
xlabel('Time')
legend('u','v')
legend boxoff