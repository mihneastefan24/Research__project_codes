% Spacecraft Guidance and Navigation
% Assignment 1 
% Exercise 1 
% Martin Mihnea Stefan 223840
% 2023/2024

clc; close all; clearvars;


%% Exercise 1. Find the x coordinates of L1
% Initialize data and function for fsolve
mu = 0.012150;
dUdx = @(x) x- (x+mu)*(1-mu)./abs(x+mu).^3 + (1-mu-x).*mu./abs(x+mu-1).^3;
options = optimoptions("fsolve",'FunctionTolerance',1e-10,'OptimalityTolerance',1e-10);
L1 = fsolve(dUdx,-2,options);
L2 = fsolve(dUdx,0,options);
L3 = fsolve(dUdx,1,options);
% Initialize the data for plots
x = [-3.5:0.01:3.5];
x1 = [-3.5:0.01:-0.028];
x2 = [0:0.01:0.98];
x3 = [1:0.01:3.5];
p = zeros(size(x));
figure(1)
for i = 1:length(x)
     p(i) = dUdx(x(i)); 
end
% Plot the Lagrange points
plot(x1,dUdx(x1),x2,dUdx(x2),x3,dUdx(x3), "Color","#0072BD","LineWidth",2);
axis([-3.5 3.5 -30 30]);
xlabel("x",'Interpreter','latex','FontWeight','bold','FontSize',14)
ylabel("$\frac{dU}{dx}$",'Interpreter','latex','FontWeight','bold','FontSize',20)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
hold on;
% Plot of planets 
plot(-mu,0,'o', 'MarkerFaceColor', 'c',"MarkerSize",10)
axis([-2 2 -30 30])
title('Lagrange Points')
hold on;
plot(1-mu,0,'o', 'MarkerFaceColor', "#D95319","MarkerSize",10)
text([1-mu-0.1 mu-0.1], [1 1], {'P1','P2'},'VerticalAlignment','bottom' )
hold on;
grid on;
plot([L1 L2 L3],[dUdx(L1) dUdx(L2) dUdx(L3)],'o',"MarkerFaceColor","#7E2F8E","MarkerSize",10)
text([L1-0.1 L2-0.1 L3-0.1],[dUdx(L1)-1.5, dUdx(L2)-1.5, dUdx(L3)-1.5],{' L3 ', ' L1 ' ' L2 '},'VerticalAlignment','top')
legend({'$\frac{dU}{dx}$','','', 'First Celestial Body', 'Second Celestial Body', 'Lagrange Points'}, 'Interpreter', 'latex');
hold on;
yL3 = dUdx(L3);
%% Part 2
syms x y z vx vy vz mu            
% Compute distances from bodies 1 and 2
r1 = sqrt((x + mu)^2 + y^2 + z^2);
r2 = sqrt((x + mu - 1)^2 + y^2 + z^2);
    
% Compute derivative of the potential
dUdx_sym = x - (1-mu)/(r1^3)*(mu+x) + mu/(r2^3)*(1-mu-x);
dUdy_sym = y - (1-mu)/(r1^3)*y - mu/(r2^3)*y;
dUdz_sym = -(1-mu)/(r1^3)*z - mu/(r2^3)*z;

% Assemble the matrix A(t)=dfdx 6x6 matrix
f = [vx; vy; vz; dUdx_sym + 2*vy; dUdy_sym - 2*vx; dUdz_sym];
dfdx_sym = jacobian(f,[x,y,z,vx,vy,vz]);
simplify(dfdx_sym);
dfdx = matlabFunction(dfdx_sym); 
dUdx = matlabFunction(dUdx_sym);
dUdy = matlabFunction(dUdy_sym);
dUdz = matlabFunction(dUdz_sym);
functions_cell = {dfdx, dUdx, dUdy, dUdz};

tic
cspice_furnsh('kernels\naif0012.tls'); % (LSK)
cspice_furnsh('kernels\de432s.bsp');   % (SPK)
cspice_furnsh('kernels\pck00010.tpc'); % (PCK)
cspice_furnsh('kernels\gm_de432.tpc'); % (PCK)

%Initalize values
mu = 0.012150;

x0 = 1.08892819445324;
y0 = 0;
z0 = 0.0591799623455459;
vx0 = 0;
vy0 = 0.257888699435051;
vz0 = 0;
xx = [x0; y0; z0;vx0;vy0;vz0];

tf = 2.0;
%Propagate the initial values
[~,~,~,xx_P]  = propagate(0,xx,-tf,mu,functions_cell,false);
[~,~,~,xx_F]  = propagate(0,xx, tf,mu,functions_cell,false);

figure(2)
hold on
grid on
plot3(xx_F(:,1),xx_F(:,2),xx_F(:,3),'LineWidth',2)
plot3(xx_P(:,1),xx_P(:,2),xx_P(:,3),'LineWidth',2)
xlabel('x [-]','Interpreter','latex','FontWeight','bold','FontSize',14)
ylabel('y [-]','Interpreter','latex','FontWeight','bold','FontSize',14)
legend(["Forward","Backward"],'Location','best',FontSize=14)

% Initial guess for the corrected velocity
vy0_new = vy0;
x0_new  = x0;

% Set parameters for the correction
err_vx0 = 1;  % Set to high value initially
err_vz0 = 1; 
Nmax = 50;    % Maximum number of iterations
tol = 1e-10;  % Tolerance for convergence
iter = 0;     % Iteration counter

% Perform pseudo-Newton differential correction
while (abs(err_vx0) > tol && abs(err_vz0)) > tol && iter < Nmax 
    % Perform propagation up to event time
    [xf, PHI, ~,xx] = propagate(0,[x0_new; y0; z0; vx0; vy0_new; vz0], tf, mu,functions_cell);
    
    % Compute the deviation in the final state (x velocity)
    err_vx0 = xf(4);
    err_vz0 = xf(6);
    delta =  [PHI(4,1),PHI(4,5);PHI(6,1),PHI(6,5)] \ [xf(4);xf(6)];
    % Compute the corrections using pseudo-Newton method
    vy0_new = vy0_new - delta(2);
    x0_new  = x0_new - delta(1);
    % Update iteration counter
    iter = iter + 1;
end
%Propagate the corrected values into an orbit
xx0_new = [x0_new; y0; z0; vx0;vy0_new;vz0];
fprintf('\nv_x(t_e)= %+.3e \n v_z(t_e) = %+.3e \n Iterations: %f \n',xf(4),xf(6),iter);
[~,~,te, xx_F, tt_F] = propagate(0,xx0_new, tf,mu,functions_cell);
[~,~,~, xx_B, tt_B] = propagate(0,xx0_new, -tf,mu,functions_cell);

% Plot the Results
figure(3)
box on
hold on
grid on;
plot3(xx_F(:,1),xx_F(:,2),xx_F(:,3),'LineWidth',3,'Color','#0072BD')
plot3(xx_B(:,1),xx_B(:,2),xx_B(:,3),'LineWidth',3,'Color','#0072BD')
xlabel('x [-]','Interpreter','latex','FontWeight','bold','FontSize',14)
ylabel('y [-]','Interpreter','latex','FontWeight','bold','FontSize',14)
zlabel('z [-]','Interpreter','latex','FontWeight','bold','FontSize',14)
hZLabel = get(gca,'ZLabel');
set(hZLabel,'rotation',0,'VerticalAlignment','middle')
view([45 45 45])
hold on;
plot3(L3,yL3,0,'o','MarkerSize',10,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
plot3(1-mu,0,0,'o', 'MarkerFaceColor', "#D95319","MarkerSize",10)
title('Halo Orbit')
legend(["Halo Orbit","","L2","Moon"], 'Interpreter', 'latex','Location','best',FontSize=14)

toc
%% Calculate the families of Halo Orbits
tic

%% Initialization
mu = 0.012150;
x0 = 1.08892819445324;
y0 = 0;
z0 = 0.0591799623455459;
vx0 = 0;
vy0 = 0.257888699435051;
vz0 = 0;
tf = 2.0;

% Define initial conditions
z0_start = z0; % Initial value of z0
z0_end = 0.34;    % Target value of z0
z0_step = (0.34-z0)/10;   % Step size for incrementing z0
z0_new  = linspace(z0,0.034,10);
% Initialize arrays to store halo orbits
halo_orbits = cell(1, 1);
xx0_new = [x0,y0,z0,vx0,vy0,vz0];

% Perform pseudo-Newton differential correction for increasing z0
for i = 1:length(z0_new)
    x0 = xx0_new(1);
    y0 = xx0_new(2);
    xx0_new(3) = z0_new(i);
    vx0 = xx0_new(4);
    vy0 = xx0_new(5);
    vz0 = xx0_new(6);

    % Reset initial guess for velocity correction
    vy0_new = vy0;
    x0_new = x0;
    
    % Set parameters for the correction
    err_vx0 = 1; % Set to high value initially
    err_vz0 = 1; 
    Nmax = 50;   % Maximum number of iterations
    tol = 1e-8;  % Tolerance for convergence
    iter = 0;    % Iteration counter
    
    % Perform pseudo-Newton differential correction
    while (abs(err_vx0) > tol && abs(err_vz0)) > tol && iter < Nmax 
        % Perform propagation up to event time
        [xf, PHI, ~, ~] = propagate(0, [x0_new; y0; z0_new(i); vx0; vy0_new; vz0], tf, mu, functions_cell);
        
        % Compute the deviation in the final state (x and z velocities)
        err_vx0 = xf(4);
        err_vz0 = xf(6);
        
        % Compute the corrections using pseudo-Newton method
        delta = [PHI(4,1), PHI(4,5); PHI(6,1), PHI(6,5)] \ [xf(4); xf(6)];
        vy0_new = vy0_new - delta(2);
        x0_new = x0_new - delta(1);
        
        % Update iteration counter
        iter = iter + 1;
    end
    
    % Store halo orbit trajectories for the current value of z0
    xx0_new = [x0_new; y0; z0_new(i); vx0; vy0_new; vz0];
    [~, ~, ~, xx_F, ~] = propagate(0, xx0_new, tf, mu, functions_cell);
    [~, ~, ~, xx_B, ~] = propagate(0, xx0_new, -tf, mu, functions_cell);
    halo_orbits{end+1} = {xx_F,xx_B};
end

% Visualize the resulting halo orbits

figure;
hold on;
grid on;
color = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F','r','b','m','k'};
plotleg = cell(2*length(z0_new)+2,1);
plotleg(1,1) = {'Halo Family'}; 
plotleg(2:2*length(z0_new),1) = {''};
plotleg(2*length(z0_new)+1,1) = {'L2'};
plotleg(2*length(z0_new)+2,1) = {'Moon'};
for i = 2:length(halo_orbits)
    plot3(halo_orbits{i}{1}(:,1), halo_orbits{i}{1}(:,2), halo_orbits{i}{1}(:,3),'LineWidth',1.75,'Color',color{i});
    plot3(halo_orbits{i}{2}(:,1), halo_orbits{i}{2}(:,2), halo_orbits{i}{2}(:,3),'--','LineWidth',1.75,'Color',color{i});
end

plot3(L3,yL3,0,'o','MarkerSize',10,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
plot3(1-mu,0,0,'o', 'MarkerFaceColor', "#D95319","MarkerSize",10)
xlabel('x [-]','Interpreter','latex','FontWeight','bold','FontSize',14);
ylabel('y [-]','Interpreter','latex','FontWeight','bold','FontSize',14);
zlabel('z [-]','Interpreter','latex','FontWeight','bold','FontSize',14);
hZLabel = get(gca,'ZLabel');
set(hZLabel,'rotation',0,'VerticalAlignment','middle')
title('Family of Halo Orbits');
legend(string(plotleg), 'Interpreter', 'latex','Location','best',FontSize=14)
hold off;
box on;
view([45 45 45])
toc

%% Auxiliary Functions
% Circular Restricted 3 Body Problem 
function [dxdt] = xyCR3BP_STM(~, xx, mu, functions_cell)
% XYCR3BP_STM  Computes the state derivative for the CR3BP with State Transition Matrix (STM)
% Inputs:
%   ~              : Unused time variable (can be ignored)
%   xx             : [42x1 vector] State vector, including position, velocity, and flattened STM
%   mu             : CR3BP parameter (mass ratio)
%   functions_cell : {cell array} Contains function handles for partial derivatives of potential function
%
% Outputs:
%   dxdt : [42x1 vector] State derivative

% Extract CR3BP parameter
mu1 = mu;

% Extract flattened State Transition Matrix (STM)
Phi = reshape(xx(7:end), 6, 6);

% Extract function handles for partial derivatives of potential function
dfdxv = functions_cell{1}(mu1, xx(1), xx(2), xx(3));
dUdxv = functions_cell{2}(mu1, xx(1), xx(2), xx(3));
dUdyv = functions_cell{3}(mu1, xx(1), xx(2), xx(3));
dUdzv = functions_cell{4}(mu1, xx(1), xx(2), xx(3));

% Compute the derivative of the STM
Phidot = dfdxv * Phi;

% Assemble right-hand side
dxdt = zeros(42, 1);
dxdt(1:3) = xx(4:6);
dxdt(4) = dUdxv + 2 * xx(5);
dxdt(5) = dUdyv - 2 * xx(4);
dxdt(6) = dUdzv;
dxdt(7:end) = Phidot(:);

end

% Propagate Function
function [xf, PHIf, tf, xx, tt] = propagate(t0, x0, tf, mu, dfdx, varargin)
% PROPAGATE  Propagates the state vector and State Transition Matrix using ODE integration
% Inputs:
%   t0        : [scalar] Initial time
%   x0        : [6x1 vector] Initial state vector
%   tf        : [scalar] Final time
%   mu        : [scalar] Gravitational parameter
%   dfdx      : [function handle] Function handle for the derivative of the state vector
%   varargin  : [optional] Additional arguments for event handling
%
% Outputs:
%   xf    : [6x1 vector] Final state vector
%   PHIf  : [6x6 matrix] Final State Transition Matrix
%   tf    : [scalar] Final time
%   xx    : [Nx6 matrix] State vectors at each integration step
%   tt    : [Nx1 vector] Time at each integration step

% Check for optional argument for event handling
if nargin > 5
    evtFlag = varargin{1};
else
    evtFlag = true;
end

% Compute time of flight
tof = tf - t0;

% Initialize State Transition Matrix at t0
Phi0 = eye(6);

% Append to initial conditions the conditions for the STM
x0Phi0 = [x0; Phi0(:)];

% Perform integration
options_STM = odeset('reltol', 1e-9, 'abstol', 1e-9, 'Events', @(t, x) y_axis_crossing(t, x, evtFlag));
[tt, xx] = ode78(@(t, x) xyCR3BP_STM(t, x, mu, dfdx), [0 tof], x0Phi0, options_STM);

% Extract state vector and State Transition Matrix
xf = xx(end, 1:6)';
PHIf = reshape(xx(end, 7:end), 6, 6);
tf = tt(end);

end

% Event function - x-axis crossing
function [value, isterminal, direction] = y_axis_crossing(~, xx, isTerminal)
% Y_AXIS_CROSSING  Event function to detect crossing of the y-axis
%
% Inputs:
%   ~           : Unused time variable (can be ignored)
%   xx          : [6x1 vector] State vector
%   isTerminal  : [logical] Flag indicating if the event is terminal
%
% Outputs:
%   value       : Value to be monitored for zero crossing (here, y-coordinate)
%   isterminal  : Flag indicating if integration should terminate when the event is triggered
%   direction   : Direction of zero crossing to detect (here, crossing from positive to negative)

% Extract y-coordinate from state vector
value = xx(2);

% Set isterminal and direction based on input flag
isterminal = isTerminal;
direction = 0; % Detect zero crossings in either direction

end

