%% Modeling and Simulation of Aerospace Systems
% Assignment #2
% Martin Mihnea Stefan 10903720

clc;
close all;
clearvars;
%% Part 1
params.Msc = 300;       % Spacecraft Mass [kg]
params.ma = 0.32;       % Accelerometer Mass [kg]
params.ba = 1.5e3;      % Accelerometer damper [Ns/m][1.5e3 2e4];
params.ka = 3e-3;       % Accelerometer spring constant [N/m][5e-5 3e-3];
params.kacc = 1;        % Acc. proportional coefficient [Vs/m]
params.Rin = 0.1;       % Inverting Resistance [ohm][0.1 10];
params.Rf = 8e4;        % Feedback Resistance  [ohm][1e4 8e4];
params.mv = 0.1;        % Spool Mass [kg]
params.kv = 1e3;        % Valve Spring [N/m]
params.bv = 1e3;        % Valve Damper [N/m]
params.alpha = 2.1e-2;  % Solenoid const [1/H]
params.beta = -60;      % Solenoid Gain  [1/Hm]
params.A0 = 4.7e-12;    % Minimum Area   [m2]
params.xvmax = 1e-5;    % Maximum Extensions valve [m]
params.k = 1.66;        % Heat Ration [-]
params.pT = 2e5;        % Tank Pressure [Pa]
params.Tt = 240;        % Tank Temperature [K]
params.R = 63.32754;    % Gas Constant [J/kgK]
params.q = 1.6e-19;     % Ion Charge [C]
params.deltav = 2000;   % Voltage Drop [V]
params.mi = 2.188e-25;  % Ion mass [kg]
params.omegas = 1.658226e-6; % Secular Pulsation [rad/s]
params.omega0 = 1.160758e-3; % Orbital Pulsation [rad/s]
% Calculate the constant term in Thrust equation 
params.Tconst = sqrt(2*params.q * params.deltav/params.mi);
% Set the length for the equation of Area variation
l = params.xvmax;
% Defin the density of the gas
params.rho = params.pT/params.R/params.Tt;
% Calculate the constant term of the mass flow ewuation
params.massconst = sqrt(params.k * params.pT^2/(params.Tt * params.R) * (2/(params.k + 1))^((params.k+1)/(params.k-1)));
% Orbital Period
T0 = 2 * pi/params.omega0;

% Set initial Values for the system
x0 = 10^-8*ones(5,1);
% Define the time interval
tspan = [0 3 * T0];
% Define the settings for the integrator
options = odeset('AbsTol',1e-9,'RelTol',1e-9);
% Implement the Integrator for the ode system
[tt,xx] = ode15s(@(t,x) system(t,x,params),tspan,x0,options);
% Plot the Solutions
figure()
sgtitle('Accelerometer State Evolution','Interpreter','latex')
subplot(2,1,1)
plot(tt,xx(:,1),'LineWidth',2.5)
legend({'Position Accelerometer[m]'},'Interpreter','latex','Location','best')
xlabel('Time[s]','Interpreter','latex')
ylabel('x[m]')
grid on;
subplot(2,1,2)
plot(tt,xx(:,3),'LineWidth',2.5)
legend({'Velocity Accelerometer[m/s]'},'Interpreter','latex','Location','best')
xlabel('Time[s]','Interpreter','latex')
ylabel('v[m/s]')
ylim([-2e-10 1e-10])
grid on;
figure()
sgtitle('Solenoid State Evolution','Interpreter','latex')
subplot(2,1,1)
plot(tt,xx(:,2),'LineWidth',2.5)
legend({'Position Solenoid[m]'},'Interpreter','latex','Location','best')
xlabel('Time[s]','Interpreter','latex')
ylabel('x[m]')
grid on;
subplot(2,1,2)
plot(tt,xx(:,4),'LineWidth',2.5)
legend({'Velocity Solenoid[m/s]'},'Interpreter','latex','Location','best')
xlabel('Time[s]','Interpreter','latex')
ylabel('v[m/s]')
ylim([-1e-8 1e-8])
grid on;

figure()
plot(tt,xx(:,5),'LineWidth',2.5)
legend({'Current[A]'},'Interpreter','latex')
title('Current Evolution','Interpreter','latex')
xlabel('Time[s]')
ylabel('I[A]')
grid on;

% Calculate the Thrust and the Drag
Av = params.A0 + params.xvmax*(params.xvmax - xx(:,2));
mdot = Av * sqrt(params.k * params.pT^2/params.R/params.Tt*(2/(params.k + 1))^((params.k + 1)/(params.k - 1)));
D = (2.2 - cos(params.omegas * tt) + 1.2 * sin(params.omega0 * tt) .* cos(params.omega0 * tt)) * 10^(-3);
T = mdot * params.Tconst;
figure()
plot(tt,T,'-','LineWidth',2.5)
hold on;
plot(tt,D,'--','LineWidth',2.5)
legend({'Thrust[N]','Drag[N]'},'Interpreter','latex','Location','best')
title('Force Evolution','Interpreter','latex')
xlabel('Time[s]','Interpreter','latex')
grid on;
hold off;

% Check if the system is stiff for the choice of the integrator
J = stifftest(params);

J_func = matlabFunction(J);
n = length(tt);
eigenvals = zeros(5,n);
stiffone = zeros(1,n);
for i = 1:n
    eigenvals(:,i) = eig(J_func(xx(end,5),xx(end,3),xx(end,2)));
    % For the system to be stiff one eigenvalue must be order of magnitudes
    % higher than the others(or two)
    eigenvalsort = sort(abs(eigenvals(:,i)),'descend');
    max1 = eigenvalsort(1);
    max2 = eigenvalsort(2);
    if max1*1e-3 >= all(abs(eigenvals(3:end,i))) && max2*1e-3 >= all(abs(eigenvals(3:end,i)))
        stiffone(i) = 1;
    else
        stiffone(i) = 0;
    end
end
if all(stiffone == 1)
    disp('The system is stiff')
end
%% SimScape
% Set the inital Values for the simscape
F1 = 9000;
F2 = 12;
x1 = 0;
x2 = 0.1;
s_max = 0.1;
Vr = 0.6;
Ir = 0.1;
kc = 0.12e6;
c = 1e4;
omegas = 1.658226e-6;
omega0 = 1.160758e-3;
kacc = params.kacc;
ka = params.ka;
ba = params.ba;
ma = params.ma;
Msc = params.Msc;
%% Plot the solution of Causal (Matlab) - already ploted above together to the solution of acausal
% model, Simscape.
% Define the Simscape model name
modelName = 'Martin10903720_Assign2_EX2_Part2'; % Replace with your model name
% Load the Simulink model
load_system(modelName);

% Run the simulation
out = sim(modelName);
% Plot the Solutions
figure()
sgtitle('Accelerometer State Evolution','Interpreter','latex')
subplot(2,1,1)
plot(tt,xx(:,1),'LineWidth',2.5)
hold on;
plot(out.ScopeData9_xa{1}.Values.Time,out.ScopeData9_xa{1}.Values.Data,'LineWidth',2.5)
title('Position','Interpreter','latex')
legend({'Causal','Acausal'},'Interpreter','latex','Location','best')
xlabel('Time[s]','Interpreter','latex')
ylabel('x[m]')
grid on;
hold off;
subplot(2,1,2)
plot(tt,xx(:,3),'LineWidth',2.5)
hold on;
plot(out.ScopeData6{1}.Values.Time,out.ScopeData6{1}.Values.Data,'LineWidth',2.5)
title('Velocity','Interpreter','latex')
legend({'Causal','Acausal'},'Interpreter','latex','Location','best')
xlabel('Time[s]','Interpreter','latex')
ylabel('v[m/s]')
ylim([-2e-10 2e-10])
grid on;
hold off;
figure()
sgtitle('Solenoid State Evolution','Interpreter','latex')
subplot(2,1,1)
plot(tt,xx(:,2),'LineWidth',2.5)
hold on;
plot(out.ScopeData11_xv{1}.Values.Time,out.ScopeData11_xv{1}.Values.Data,'LineWidth',2.5,'LineStyle','--')
title('Position','Interpreter','latex')
legend({'Causal','Acausal'},'Interpreter','latex','Location','best')
xlabel('Time[s]','Interpreter','latex')
ylabel('x[m]')
grid on;
hold off;
subplot(2,1,2)
plot(tt,xx(:,4),'LineWidth',2.5)
hold on;
plot(out.ScopeData12_vv{1}.Values.Time,out.ScopeData12_vv{1}.Values.Data,'LineWidth',2.5,'LineStyle','--')
title('Velocity','Interpreter','latex')
legend({'Causal','Acausal'},'Interpreter','latex','Location','best')
xlabel('Time[s]','Interpreter','latex')
ylabel('v[m/s]')
ylim([-4e-9 4e-9])
grid on;
hold off;

figure()
plot(tt,xx(:,5),'LineWidth',2.5)
hold on;
plot(out.ScopeData7_i{1}.Values.Time,out.ScopeData7_i{1}.Values.Data,'LineWidth',2.5)
legend({'Causal','Acausal'},'Interpreter','latex')
title('Current Evolution','Interpreter','latex')
xlabel('Time[s]')
ylabel('I[A]')
grid on;
hold off;

figure()
plot(tt,T-D,'LineWidth',2.5)
hold on;
plot(out.ScopeData3_TD{1}.Values.Time,out.ScopeData3_TD{1}.Values.Data,'LineWidth',2.5)
legend({'Causal','Acausal'},'Interpreter','latex')
title('T-D Evolution','Interpreter','latex')
xlabel('Time[s]')
ylabel('I[A]')
grid on;
hold off;
%%

function [dxdt] = system(t,x,params)
% Function that sets the mathemoatical model of the system
% Input: t: time
%        x: initial data od the system
%        params: Structure containing the fixed data od the system
% Output: dxdt: ODE for the ode15s
    xa = x(1);      % Position of Accelerometer [m]
    xv = x(2);      % Position of Solenoid      [m]
    va = x(3);      % Velocity of Accelerometer [m/s]
    vv = x(4);      % Velocity of Solenoid      [m/s]
    I = x(5);       % Current [A]
% Reinitialize data - This step could be ommited, it is done just for
% clarity in the equations
Msc = params.Msc;
ma = params.ma;
ba = params.ba;
ka = params.ka;
kacc = params.kacc;
Rin = params.Rin;
Rf = params.Rf;
mv = params.mv;
kv = params.kv;
bv = params.bv;
alpha = params.alpha;
beta = params.beta;
A0 = params.A0 ;
xvmax = params.xvmax;
l = xvmax;
k = params.k;
pt = params.pT;
Tt = params.Tt;
R = params.R;
q = params.q;
deltav = params.deltav;
mi = params.mi;
omegas = params.omegas;
omega0 = params.omega0;
    % Define the Dynamical Systems
    % Velocities
     dxadt = va;
     dxvdt = vv;

   % Intermediate Equations
   Av = A0 + l*(l - xv);
   D = (2.2 - cos(omegas * t) + 1.2 * sin(omega0 * t) * cos(omega0 * t)) * 10^(-3);
   T = sqrt(2 * q * deltav/ mi) * sqrt(k * pt/(Tt * R) * pt * (2/(k + 1))^((k+1)/(k-1))) * Av;
   dLdxv = -beta/(alpha + beta * xv)^2;
    
   % Acclerations
    dvvdt = (- kv * xv - bv * vv + 1/2 * I^2 * dLdxv)/mv;
    dvadt =  (T - D)/Msc - (ba * va + ka * xa)/ma;

   % Current
    dIdt = (alpha + beta * xv) *(- Rf/Rin) * kacc * va;

    dxdt = [dxadt; dxvdt; dvadt; dvvdt; dIdt];
end

function [J] = stifftest(params)
% Function To check the stiffness of the system
% Input: params: Structure containing the fixed data od the system
% Output: J: Jacobian Matrix of the system
% Initialize Data
Msc = params.Msc;
ma = params.ma;
ba = params.ba;
ka = params.ka;
kacc = params.kacc;
Rin = params.Rin;
Rf = params.Rf;
mv = params.mv;
kv = params.kv;
bv = params.bv;
alpha = params.alpha;
beta = params.beta;
A0 = params.A0 ;
xvmax = params.xvmax;
l = xvmax;
k = params.k;
pt = params.pT;
Tt = params.Tt;
R = params.R;
q = params.q;
deltav = params.deltav;
mi = params.mi;
omegas = params.omegas;
omega0 = params.omega0;

syms va vv xa xv I t % Symbolic Representation
% Create the system of ODEs
f = [va;
    vv;
    ((sqrt(2 * q * deltav/ mi) * sqrt(k * pt/(Tt * R) * pt * (2/(k + 1))^((k+1)/(k-1))) * (A0 + l*(l - xv))) - ((2.2 - cos(omegas * t) + 1.2 * sin(omega0 * t) * cos(omega0 * t)) * 10^(-3)))/Msc - (ba * va + ka * xa)/ma;
    (- kv * xv - bv * vv + 1/2 * I^2 * (-beta/(alpha + beta * xv)^2))/mv;
    (alpha + beta * xv) *(- Rf/Rin) * kacc * va]; % System of equations
% Based on the Control Matrix calculate the Jacobian
    J = jacobian(f, [xa xv va vv I]);  % Jacobian matrix
end