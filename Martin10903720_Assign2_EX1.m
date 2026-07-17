% Modeling and Simulation of Aeerospace Systems
% Assignment 2
% Martin Mihnea Stefan
% Exercise 1

clearvars;
close all;
clc;
%%
% Part 1
% Data Initialization
tic
% Electrical parameters
params.R = 0.1;        % Resistance [Ohm]
params.L = 0.001;      % Inductance [H]
params.km = 0.3;       % Motor Constant[N*m/A]
params.mr = 0.2;       % Radiator  mass [kg]
params.Lr = 0.5;       % Radiator length [m]
% Bodies number 
% 1- Main Body
% 2-3 - Solar Panels, one on each side
% 4-5 - Radiators, one on each side
% Thermal Parameters
params.Psun = 1350;                     % Sun Radiation [W/m^2]
params.C = [1.5e5 1187.5 1187.5 30 30]; % Row matrix of heat capacity for each body[J/K]
params.G = [10 10 10 10];               % Row matrix for  each thermal conductance btw the 5 bodies [J/K]
params.alpha = [0.6 0.78 0.78];         % Absorbtivity row matrix [-]
params.eps = [0.45 0.75 0.75];          % Emissivity for bodies [1-3]
params.sigma = 5.67 * 1e-8;             % Stefan Boltzamann Constant[W/m^2K^4]
params.kp = 1e-4;                       % Proportional Gain

% Areas of Each Body
params.A1e = 0.5 * 0.5 + 4 * 0.5 * 1.5; % Emmitting area of node 1
params.A2e = 0.5 * 0.95;                % Emmitting area of node 2
params.A3e = 0.5 * 0.95;                % Emmitting area of node 3
params.A4e = 0.5 * 0.5;                 % Emmitting area of node 4
params.A5e = 0.5 * 0.5;                 % Emmitting area of node 5
params.A1a = 0.5 * 0.5;                 % Absorbing area of node 1
params.A2a = 0.95 * 0.5;                % Absorbing area of node 2
params.A3a = 0.95 * 0.5;                % Absorbing area of node 3

% Limits
params.eps_max = 0.98;                  % Maximum Emissivity coefficient for radiators [-]
params.eps_min = 0.01;                  % Minimum Emissivity coefficient for radiators [-]
params.Tds = 3;                         % Deep Space Temperature [K]
params.theta_min = -0.4 * pi;           % Minimum Inclination of radiators [rad]
params.theta_max =  0;                  % Maximum Inclination of radiators [rad]
params.T1_min = 290;                    % Main Body Lower Temperature Limit [K]
params.T1_max = 300;                    % Main Body Upper Temperature Limit [K]
params.T1ref = 294.15;                  % Main Body Reference(Nominal) Temperature [K]


% Initial Values
T0 = 298.15 * ones(5,1);                % Initial Temperature for all bodies [K]

t0 = 0;                                 % Initial Time [s]
tmid = 10 * 3600;                       % Time at 10 Hours [s]
tend = 50 * 3600;                       % Final Time for simulation [s]
tspan = [t0 tend];                      % Time Interval
theta0 = -0.4 *pi;                      % Initial angle for radiators [rad]
i0 = 0;                                 % Initial current passing through circuit [A]
omega0 = 0;                             % Initial angular velocity for radiators [rad/s]
params.tend = tend;
% Column matrix for initial Values
x0 = [T0;i0;theta0;omega0];

[tt,xx] = odefunction(x0,tspan,params);


params.r1 = params.T1ref * 1.001;           % The temperature range upper limit around nominal temperature at t= 10h [K] = +0.1 % Tnominal
params.r2 = params.T1ref * 0.999;           % The temperature range lower limit around nominal temperature at t= 10h [K] = -0.1 % Tnominal
hours = 0:2.5:50;                           % Time Vector for Ploting Hours

%% Plot
% Plot the Temperatures
figure(1)
hold on;
plot(tt, xx(:,1:5),'LineWidth',2)
plot(tt, 300*ones(size(tt)),'--k','LineWidth',2)
plot(tt, 290*ones(size(tt)),'--k','LineWidth',2)
plot(tt, params.r1 *ones(size(tt)),':r','LineWidth',2)
plot(tt, params.r2 *ones(size(tt)),':r','LineWidth',2)
plot(10 * 3600 * ones(1,5),linspace(280,330,5),'--c','LineWidth',2.5)
xticks(hours * 3600);
xticklabels({hours});
ylim([284.5 325.5])
legend({'T1','','T2-T3','','T4-T5','Lower and Upper Bounds','','Lower and Upper Limits for 10H','','10H'},'Location','westoutside');
title('Temperature Evolution for each body in time','Interpreter','latex')
xlabel('Time[h]')
ylabel('Temperature [K]')
grid on;
hold off;
% zoom on the initial part of the graph where the error is larger
axes('position',[.65 .475 .25 .25])
box on
plot(tt,xx(:,1),'LineWidth',3)
hold on;
plot(tt, params.r1 *ones(size(tt)),':r','LineWidth',2)
plot(tt, params.r2 *ones(size(tt)),':r','LineWidth',2)
axis([7.5*3600 12.5*3600  params.r2-0.1 params.r1+0.1]);
xticks(hours * 3600);
xticklabels({hours});
grid on;
hold off;
% Plot the Angular velocity
figure(2)
plot(tt,xx(:,8),'LineWidth',2, 'DisplayName','$\omega$')
title('Angular Velocity','Interpreter','latex')
ylabel('$\omega$[rad/s]','Interpreter','latex')
xlabel('Time[h]')
legend('Location','best','Interpreter','latex')
grid on;
xticks(hours * 3600);
xticklabels({hours});

% Plot the radiator angles
figure(3)
plot(tt,zeros(size(tt)),'LineWidth',2,'DisplayName','Upper Limit')
hold on;
plot(tt, -0.4 *pi * ones(size(tt)),'LineWidth',2,'DisplayName','Lower Limit')
plot(tt,xx(:,7),'m','LineWidth',2,'DisplayName', '$\theta$')
title('Angular Position of Radiators','Interpreter','latex')
xticks(hours * 3600);
xticklabels({hours});
legend({'','','$\theta$'},'Interpreter','latex','Location','Best')
xlabel('Time [h]')
ylabel('$\theta$[rad]','Interpreter','latex')
grid on
hold off;

texe = toc;
disp(" Computation time is " + texe);
%% Modelica/Dymola Solution 
% Load the Data from modelica CSV to the MATLAB editor
modelica = readtable('OpenModelicaResults.csv');
time = table2array(modelica(2:end,1));   % Modelica time steps [s]
Tm = zeros(length(time),5);              
theta = zeros(length(time),1);
Tm(:,1) = table2array(modelica(2:end,3)); % Main Body Temperature [K]
Tm(:,2) = table2array(modelica(2:end,4)); % Panel 1 Temperature [K]
Tm(:,3) = table2array(modelica(2:end,5)); % Panel 2 Temperature [K]
Tm(:,4) = table2array(modelica(2:end,6)); % Radiator 1 Temperature [K]
Tm(:,5) = table2array(modelica(2:end,7)); % Radiator 2 Temperature [K]
theta(:) = table2array(modelica(2:end,2)); % Radiator angular position [rad]


% Plot the solution
% Plot the Temperature evolution

figure()
hold on;
plot(time, Tm(:,1:5),'LineWidth',2)
plot(tt, 300*ones(size(tt)),'--k','LineWidth',2)
plot(tt, 290*ones(size(tt)),'--k','LineWidth',2)
plot(tt, params.r1 *ones(size(tt)),':r','LineWidth',2)
plot(tt, params.r2 *ones(size(tt)),':r','LineWidth',2)
plot(10 * 3600 * ones(1,5),linspace(280,330,5),'--c','LineWidth',2.5)
xticks(hours * 3600);
xticklabels({hours});
ylim([284.5 325.5])
legend({'T1','','T2-T3','','T4-T5','Lower and Upper Bounds','','Lower and Upper Limits for 10H','','10H'},'Location','best');
xlabel('Time[h]')
ylabel('Temperature [K]')
grid on;
hold off;
figure()
hold on;
plot(time, Tm(:,1),'LineWidth',2)
plot(tt, xx(:,1),'LineWidth',2,'LineStyle','--')
xticks(hours * 3600);
xticklabels({hours});
ylim([290 300])
legend({'T1-acausal','T1-causal'},'Location','best');
title('Comparison Between Causal and Acausal Modeling')
xlabel('Time[h]')
ylabel('Temperature [K]')
grid on;
hold off;

% Plot the radiator angles
figure()
plot(time,zeros(size(time)),'LineWidth',2,'DisplayName','Upper Limit')
hold on;
plot(time, -0.4 *pi * ones(size(time)),'LineWidth',2,'DisplayName','Lower Limit')
plot(time,theta,'m','LineWidth',2,'DisplayName', '$\theta$')
title('Angular Position of Radiators - Acausal','Interpreter','latex')
xticks(hours * 3600);
xticklabels({hours});
legend({'','','$\theta$'},'Interpreter','latex','Location','Best')
xlabel('Time [h]')
ylabel('$\theta$[rad]','Interpreter','latex')
grid on
hold off;

% Plot the radiator angles comparison Causal vs Acausal
figure()

hold on;
plot(tt, xx(:,7),'LineWidth',2,'DisplayName','Lower Limit')
plot(time,theta,'m','LineWidth',2,'LineStyle','--','DisplayName', '$\theta$')
title('Angular Position of Radiators - Causal vs Acausal','Interpreter','latex')
xticks(hours * 3600);
xticklabels({hours});
legend({'$\theta$ Causal','$\theta$ Acausal'},'Interpreter','latex','Location','Best')
xlabel('Time [h]')
ylabel('$\theta$[rad]','Interpreter','latex')
grid on
hold off;
texe = toc;

%%
function [dxdt] = system2(~,x,params)
% Function that includes all the ODEs into the Matlab integrator,
% describing the mathematical model of the system
% Input:  x: Inital Values for each integration step
%         params:  Structure of data containing the initial fixed data
% Output: dxdt: The integrated values of the system
% Initialization
    T = x(1:5);         % Initial Temperatures [K]
    i = x(6);           % Initial Current   [A]
    theta = x(7);       % Initial Angle [rad]
    omega = x(8);       % Initial Angular Velocity [rad/s]

    % T1 just has to be checked, the limits should be imposed on the theta and see if they apply to temperature 

    % Eps(theta)
    eps = params.eps_min + (params.eps_max - params.eps_min)/(0.4*pi) * (theta + 0.4 * pi);
    
    % Mechanical System
    Jr = 1/3 * params.mr * 2 * params.Lr^2;
    domegadt = params.km * i/Jr;
    dthetadt = omega;

    % Electrical System
    Vin = params.kp * (T(1) - params.T1ref);
    didt = (Vin - params.R * i - params.km * dthetadt)/params.L;

    % Conditions
    % Imposed condition to make the temperature return to reference value
    % By fixing radiators in position 
    if  theta >= 0 && T(1) > params.T1ref       
        domegadt = 0;
        dthetadt = 0;
        didt = -i;
    end

    if  theta < -0.4*pi && T(1) < params.T1ref
        domegadt = 0;
        dthetadt = 0;
    end

    % Thermal System
    dTdt = zeros(5,1);
    dTdt(1) = 1/params.C(1) * (params.Psun * params.alpha(1) * params.A1a + params.G(1) * (T(2) - T(1)) + params.G(2) *(T(3) - T(1)) ...
             + params.G(3) * (T(4) -T(1)) + params.G(4) * (T(5) - T(1)) - params.eps(1) * params.sigma * params.A1e * (T(1)^4 - params.Tds ...
             ^4));
    dTdt(2) = 1/params.C(2) * (params.Psun * params.alpha(2) * params.A2a + params.G(1) * (T(1) - T(2)) - params.eps(2) * params.sigma ...
             * params.A2e * (T(2) ^ 4 - params.Tds ^ 4));
    dTdt(3) = 1/params.C(3) * (params.Psun * params.alpha(3) * params.A3a + params.G(2) * (T(1) - T(3)) - params.eps(3) * params.sigma ...
             * params.A3e * (T(3) ^ 4 - params.Tds ^4));
    dTdt(4) = 1/params.C(4) * (-eps * params.sigma * params.A4e * (T(4) ^ 4 - params.Tds ^ 4) + params.G(3) * (T(1) - T(4)));
    dTdt(5) = 1/params.C(5) * (-eps * params.sigma * params.A5e * (T(5) ^ 4 - params.Tds ^ 4) + params.G(4) * (T(1) - T(5)));
    
    dxdt = [dTdt; didt;dthetadt;domegadt];

end


function [tt,xx] = odefunction(x0, tspan, params)
    % Function that performs the Integration of the system
    % Input : x0: Initial Values for the Simulation explained above
    %         tspan: simulation time interval [s]  
    %         params:  Structure of data containing the initial fixed data
    % Output: tt: vector containing times for each integration
    %         xx(1:5): Temperatures of each body obtained from the integrator [K]
    %         xx(6): Current at each time step of integration [A]
    %         xx(7): Radiators'angle at each time step of integration [rad]
    %         xx(8): Radiators'angular velocity at each time step of integration [rad/s]
    
    % Set Integration Parameters
    options = odeset('AbsTol',1e-6,'RelTol',1e-6);
    [tt, xx] = ode15s(@(t,x) system2(t,x,params), tspan, x0,options);
end

