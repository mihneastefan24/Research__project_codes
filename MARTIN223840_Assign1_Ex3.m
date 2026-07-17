%% Spacecraft Guidance and Navigation 2023/2024
% Assignment 1 Exercise 3
% Martin Mihnea Stefan

clc; close all; clearvars;
rng default;

%% Initialization

addpath('.\kernels')
addpath('.\mice\src\mice')
addpath('.\mice\lib')

cspice_furnsh('kernels\naif0012.tls'); 
cspice_furnsh('kernels\de432s.bsp');   
cspice_furnsh('kernels\gm_de432.tpc'); 
cspice_furnsh('kernels\pck00010.tpc');

%% Part 2 
tic
% Initialization
xx0 = zeros(15,1);
t0 = cspice_str2et('2023-05-28 14:13:09.000 UTC');

r_earth = cspice_spkezr('EARTH',t0,'ECLIPJ2000','NONE','SUN');
xx0(1:6) = r_earth;

Tmax = 0.8*1e-3;
Isp = 3120;
m0  = 1000;
g0  = 9.81*1e-3;

mu = cspice_bodvrd('SUN','GM',1);
LU = cspice_convrt(1, 'AU' , 'km');
t  = (LU^3/mu) ^ 0.5;               
% Adminesionalize the variables wrt AU and the Earth's period around the
% Sun
SAVEVAL.Tmax = Tmax*t^2/LU/1000;
SAVEVAL.Isp  = Isp/t;
SAVEVAL.g0   = g0*t^2/LU;
SAVEVAL.mu   = mu*t^2/LU^3;
SAVEVAL.t0   = t0/t;
SAVEVAL.m0   = m0/1000;
SAVEVAL.t    = t;
SAVEVAL.LU   = LU;
SAVEVAL.m    = 1000;
xx0(7) = m0; 
% Create the initial Guess for the costates
lambda = zeros(8,1);
lambda(1:6) = 40*rand(6,1) - 20;
lambda(7)   = 10*abs(rand(1)) ;
lambda(8)   = rand(1)*2*pi + t0/t;

xx0(1:7) = [xx0(1:3)/LU;xx0(4:6)*t/LU;xx0(7)/1000];


%% Part 3 


iter = 0;
exitflag = 0;
   
while iter < 50  
    lambda(1:6) = 40*rand(6,1) - 20;
    lambda(7)   = 20*abs(rand(1)) ;
    lambda(8)   = rand(1)*2*pi + t0/t;
        % Solve the zero finding problem
    [SOL_FS,~,exitflag] = fsolve(@propagator,lambda,[],xx0(1:7),SAVEVAL); 
    
    iter = iter + 1;

     if exitflag == 1
         break;
     end
end
% Propagate the SC from t0 to the final time
opt = odeset('AbsTol',1e-12,'RelTol',1e-12);
[TF,XF] = ode113(@(t,y) dynamicsSC(t,y,SAVEVAL),[t0/t SOL_FS(8)],[xx0(1:7);SOL_FS(1:7)],opt);
% Check the Hamiltonian for Error
H = zeros(1,length(TF));
for j = 1:length(TF)
    x_l_f = dynamicsSC(TF(j),XF(j,:)',SAVEVAL);
    H(j) = 1 + dot(XF(j,8:end),x_l_f(1:7));
end
% Compute the states of Earth and Venus at final time
VEN = cspice_spkezr('VENUS',TF(end)*t,'ECLIPJ2000','NONE','SUN');
VEN_plot = cspice_spkpos('VENUS',TF'*t,'ECLIPJ2000','NONE','SUN')/LU;
EARTH_plot = cspice_spkpos('EARTH',TF'*t,'ECLIPJ2000','NONE','SUN')/LU;
err = [XF(end,1:3)' * LU;XF(end,4:6)' * LU/t] - VEN;
% Calculate the time of flight and use cspice to transform them into dates
tf  = cspice_et2utc(SOL_FS(8)*t,'C',3);
tof = days(seconds(SOL_FS(8)*t - t0));

% Plot the Results
figure(1)
hold on;
grid on;
plot3(VEN_plot(1,:),VEN_plot(2,:),VEN_plot(3,:),'LineWidth',2.5,'Color','r','DisplayName','Venus Trajectory')
plot3(VEN(1)/LU,VEN(2)/LU,VEN(3)/LU,'o','MarkerFaceColor',[0.44 0.55 0.66],'MarkerSize',5,'DisplayName','Venus Position at Interception')
plot3(EARTH_plot(1,:),EARTH_plot(2,:),EARTH_plot(3,:),'Color','#7E2F8E','LineWidth',2.5,'DisplayName','Earth Trajectory')
plot3(XF(:,1),XF(:,2),XF(:,3),'LineWidth',2.5,'Color',[0 0.4470 0.7410],'DisplayName','Spacecraft Trajectory')
plot3(EARTH_plot(1,1),EARTH_plot(2,1),EARTH_plot(3,1),'v','MarkerFaceColor','#7E2F8E','MarkerSize',5,'DisplayName','Earth Position at Departure')
plot3(0,0,0,'o','MarkerFaceColor','#EDB120','DisplayName','SUN')
view(30, 20)
plot3(XF(:,1),XF(:,2),XF(:,3),'LineWidth',2.5,'Color',[0 0.4470 0.7410],'DisplayName','Spacecraft Trajectory')
legend('Location','best','Orientation','vertical');
xlabel('X - Distance [AU]','FontName','Times New Roman','FontSize',12)
ylabel('Y - Distance [AU]','FontName','Times New Roman','FontSize',12)
zlabel('Z - Distance [AU]','FontName','Times New Roman','FontSize',12)
axis equal
hold off;
% Plot the Hamiltonian Error
figure(2)
plot(TF,(H - H(1))/H(1))
grid on;
xlabel('Time of Flight')
ylabel('Relative Error')
title('Hamiltonian Relative Error')
toc

%% Part 4
rng default;
% Create a vector of Thrust from the highest to smallest
Tmaxmin = linspace(800,500,10)*1e-6;

figure(3)
hold on;
%Put the before code in a for loop for each Thrust value for numerical
%continuation
% Create for loop to go through each Thrust value
for k = 1:length(Tmaxmin) 

    SAVEVAL.Tmax = Tmaxmin(k)*t^2/(m0*LU);
    
    iter = 0;
    % Use numerical continuation
    lambda = SOL_FS;
    exitflag = 0;
    while iter < 50 
        if k == 1
            lambda(1:6) = 40*rand(6,1) - 20;
            lambda(7)   = 20*abs(rand(1)) ;
            lambda(8)   = rand(1)*2*pi + t0/t;
        end
        % Solve the zero finding problem for each 
        [SOL_FS,~,exitflag] = fsolve(@(y) propagator(y,xx0(1:7),SAVEVAL),lambda); 
    
         iter = iter + 1;
         if exitflag == 1
            break;
         end
     end
    % Propagate the Values for each Thrust
    [TF,XF] = ode78(@(t,y) dynamicsSC(t,y,SAVEVAL),[t0/t SOL_FS(8)],[xx0(1:7);SOL_FS(1:7)],opt);
    % Hamiltonian Check up for each Thrust values
    H = zeros(1,length(TF));
    for j = 1:length(TF)
        OPT = dynamicsSC(TF(j),XF(j,:)',SAVEVAL);
        H(j) = 1 + dot(XF(j,8:end),OPT(1:7));
    end
    % Compute the position of Venus and Earth at final times
    VEN = cspice_spkezr('VENUS',TF(end)*t,'ECLIPJ2000','NONE','SUN');
    VEN_plot = cspice_spkpos('VENUS',TF'*t,'ECLIPJ2000','NONE','SUN')/LU;
    EARTH_plot = cspice_spkpos('EARTH',TF'*t,'ECLIPJ2000','NONE','SUN')/LU;
    
    % Compute the error between the final state and the Venus states at
    % final time
    err = [XF(end,1:3)' * LU;XF(end,4:6)' * LU/t] - VEN;
    % Compute the final time and time of flight
    tf  = cspice_et2utc(SOL_FS(8)*t,'C',3);
    tof = days(seconds(SOL_FS(8)*t - t0));
end
grid on;
% Plot the results
plot3(XF(:,1),XF(:,2),XF(:,3),'LineWidth',2.5,'Color','#77AC30','DisplayName','Spacecraft Trajectory')
hold on;
plot3(VEN_plot(1,:),VEN_plot(2,:),VEN_plot(3,:),'LineWidth',2.5,'Color','r','DisplayName','Venus Trajectory')
plot3(VEN_plot(1,1),VEN_plot(2,1),VEN_plot(3,1),'diamond','MarkerFaceColor',[0.44 0.55 0.66],'MarkerSize',5,'DisplayName','Venus Position at departure')
plot3(VEN(1)/LU,VEN(2)/LU,VEN(3)/LU,'o','MarkerFaceColor',[0.44 0.55 0.66],'MarkerSize',5,'DisplayName','Venus Position at Interception')
plot3(EARTH_plot(1,:),EARTH_plot(2,:),EARTH_plot(3,:),'Color','#7E2F8E','LineWidth',2.5,'DisplayName','Earth Trajectory')
plot3(EARTH_plot(1,1),EARTH_plot(2,1),EARTH_plot(3,1),'v','MarkerFaceColor','#7E2F8E','MarkerSize',5,'DisplayName','Earth Position at Departure')
plot3(0,0,0,'o','MarkerFaceColor','#EDB120','DisplayName','SUN')
view(31, 20)
legend('Location','best','Orientation','vertical');
xlabel('X - Distance [AU]','FontName','Times New Roman','FontSize',12)
ylabel('Y - Distance [AU]','FontName','Times New Roman','FontSize',12)
zlabel('Z - Distance [AU]','FontName','Times New Roman','FontSize',12)
axis equal
hold off;

toc

%% Functions

function [x_l_f] = dynamicsSC(~, x_l_i, data)
% DYNAMICSSC Computes the dynamics of the spacecraft in the TPBVP, in
% continous guidance (Continous Thrust) such that the final results will
% match the states of another celestial body of interest

% Inputs:
%   x_l_i : Initial state vector including the states and costates [14x1]
%   data  : Struct containing relevant data for the dynamics computation
%
% Outputs:
%   x_l_f : Final state vector including the states and the costates [14x1]

% Extract relevant data from initial state vector
xx = x_l_i(1:7);

% Extract parameters from data struct
T    = data.Tmax;
Isp  = data.Isp;
mu   = data.mu;
g0   = data.g0;
m    = x_l_i(7);
u    = 1;

% Compute velocity and mass rate
r_dot = xx(4:6);
v_dot = -mu / norm(xx(1:3))^3 * xx(1:3) - u * T / m * x_l_i(11:13) / norm(x_l_i(11:13));
m_dot = -u * T / (Isp * g0);

% Final state vector
xx_f = [r_dot; v_dot; m_dot];

% Compute adjoint variables
l_r_dot = -3 * mu * xx(1:3) * dot(xx(1:3), x_l_i(11:13)) / norm(x_l_i(1:3))^5 + mu * x_l_i(11:13) / norm(x_l_i(1:3))^3;
l_v_dot = -x_l_i(8:10);
l_m_dot = -u * norm(x_l_i(11:13)) * T / m^2;
lambda_f = [l_r_dot; l_v_dot; l_m_dot];

% Final state vector including adjoint variables
x_l_f = [xx_f; lambda_f];

end

function OPT = propagator(lambda, xx0, SAVEVAL)
% PROPAGATOR Is used in order to find the initial costate
% Inputs:
%   guess : Initial guess for the optimization problem
%   xx0   : Initial state vector
%   SAVEVAL  : Struct containing relevant data for the propagation
%
% Outputs:
%   OPT   : Vector of optimization constraints

% Extract relevant data from SAVEVAL struct
Tmax = SAVEVAL.Tmax;
Isp  = SAVEVAL.Isp;
mu   = SAVEVAL.mu;
g0   = SAVEVAL.g0;
T_conv = SAVEVAL.t;
LU   = SAVEVAL.LU;
m0   = SAVEVAL.m0;
m    = SAVEVAL.m;
t0   = SAVEVAL.t0;
tf   = lambda(end);

% Compute Venus position and velocity
Ven = cspice_spkezr('VENUS', tf * T_conv, 'ECLIPJ2000', 'NONE', 'SUN');
VEN = [Ven(1:3) / LU; Ven(4:6) * T_conv / LU];

% Append guess to initial state vector
xx0 = [xx0; lambda(1:end-1)];

% Set options for ODE solver
opt = odeset('AbsTol', 1e-11, 'RelTol', 1e-12);

% Propagate dynamics using ODE solver
[~, xf] = ode78(@(t,y) dynamicsSC(t, y, SAVEVAL), [t0, tf], xx0, opt);

% Extract Lagrange multiplier and final mass
lambda = xf(end, 8:14)';
m = xf(end, 7);

% Compute control input and mass rate
u = 1;
x_dot = xf(end, 4:6)';
v_dot = -mu * xf(end, 1:3)' / norm(xf(end, 1:3)')^3 - u * Tmax / m * xf(end, 11:13)' / norm(xf(end, 11:13)');
m_dot = -u * Tmax / Isp / g0;

% Construct state vector
x = [x_dot; v_dot; m_dot];

% Compute Hamiltonian
H = 1 + dot(lambda, x);

% Define optimization constraints
OPT = [xf(end, 1:3)' - VEN(1:3);
       xf(end, 4:6)' - VEN(4:6);
       lambda(7);
       H - dot(lambda(1:6), [VEN(4:6); -mu / norm(VEN(1:3))^3 * VEN(1:3)])];
end


% function St = STP(lambdav,lambdam,m,Isp,g0)
%  St = -norm(lambdav/m*Isp*g0 - lambdam);
% end