% Spaceraft Guidance and Navigation 2023/2024
% Assignement 2 Exercise 1
% Martin Mihnea Stefan

clearvars; close all; clc;
%de asta sunt sigur acum ca totul este in regula scrie comentariile si scrie modifica catefa functii si linii si numele constantelor
addpath('.\kernels\')
addpath('.\sgp4\')
addpath('.\tle\')
addpath('.\mice\lib')
addpath('.\mice\src\mice')

cspice_furnsh('assignment02.tm');
%% Initialization

% Convert launch and separation time strings to ephemeris time (ET)
et_launch = cspice_str2et('2010-07-15');
et_sep = cspice_str2et('2010-08-12 05:27:39.114 UTC');

% Define constants
alpha = 0.1;
beta = 2;
N = 10;

% Initial positions and velocities for the two spacecraft
r_01 = [4622.232026629, 5399.3369588058, -0.0212138165769957];
r_02 = [4621.69343340281, 5399.26386352847, -3.09039248714313];
v_01 = [0.812221125483763, -0.721512914578826, 7.42665302729053]; 
v_02 = [0.813960847513811, -0.719449862738607, 7.42706066911294];

% Combine positions and velocities into state vectors
x_01 = [r_01, v_01];
x_02 = [r_02, v_02];

% Initial covariance matrix
P0 = [5.6e-7, 3.5e-7, -7.1e-8, 0, 0, 0;
      3.5e-7, 9.7e-7, 7.6e-8, 0, 0, 0;
      -7.1e-8, 7.6e-8, 8.1e-8, 0, 0, 0;
      0, 0, 0, 2.8e-11, 0, 0 ;
      0, 0, 0, 0, 2.7e-11, 0; 
      0, 0, 0, 0, 0, 9.6e-12];

% Convert separation time string to ephemeris time (ET)
t_sep = cspice_str2et('2010-08-12T05:27:39.114');

%% Part 1

% Get the gravitational constant of Earth
mu_E = cspice_bodvrd('EARTH','GM',1);

% Calculate magnitudes of position and velocity vectors
rn_1 = norm(r_01);
rn_2 = norm(r_02);
vn_1 = norm(v_01);
vn_2 = norm(v_02);

% Calculate semi-major axis for spacecraft MANGO
an_1 = mu_E/(-vn_1^2 + 2*mu_E/rn_1);
a_sat1 = 1/(2/rn_1 - vn_1^2/mu_E);
% Calculate orbital period for spacecraft MANGO
T1 = 2*pi*sqrt(an_1^3/mu_E);

% Create time grid for propagation
t_grid = t_sep:T1:t_sep + N*T1;

% Initialize arrays to store propagated states and covariance matrices for
% MANGO and TANGO
P1_prop = zeros(6,6,length(t_grid));
P2_prop = zeros(6,6,length(t_grid));
x1_prop = zeros(6,length(t_grid));
x2_prop = zeros(6,length(t_grid));

% Set initial values for states and covariance matrices
P1_prop(:,:,1) = P0;
P2_prop(:,:,1) = P0;
x1_prop(:,1) = x_01';
x2_prop(:,1) = x_02';

% Propagate states and update covariance matrices
for i = 2:length(t_grid)
    [x1_prop(:,i), tt_1, xx_1, PHIf1] = keplerian_propagator_STM(t_grid(i-1), x1_prop(:,i-1), t_grid(i), mu_E);

    P1_prop(:,:,i) = PHIf1 * P1_prop(:,:,i-1) * PHIf1';

    [x2_prop(:,i), tt_2, xx_2, PHIf2] = keplerian_propagator_STM(t_grid(i-1), x2_prop(:,i-1), t_grid(i), mu_E);

    P2_prop(:,:,i) = PHIf2 * P2_prop(:,:,i-1) * PHIf2';
end

% Unscented Transformation (UT)

n = 6;          
k = 0;

lambda = alpha^2 * n - n;

% Define the Weights for (mean and covariance weights)
W0_m   = lambda/(n+lambda);
W0_c   = lambda/(n + lambda) + (1 - alpha^2 + beta);
Wi   = 1/(2 * (n + lambda));
Wm = [W0_m;Wi*ones(2*n,1)]';
Wc = [W0_c;Wi*ones(2*n,1)]';

% Define the Sigma Matrix
sigma_matrix = sqrtm((n + lambda)*P0) ;

% Define the Sigma Points
sigma_points_01 = [x1_prop(:,1), x1_prop(:,1) + sigma_matrix(:,1:6),x1_prop(:,1) - sigma_matrix(:,1:6)];
sigma_points_02 = [x2_prop(:,1), x2_prop(:,1) + sigma_matrix(:,1:6),x2_prop(:,1) - sigma_matrix(:,1:6)];

% Initialize the Unscented Transform covariance matrices for TANGO and
% MANGO
P1_prop_UT = zeros(6,6,length(t_grid));
P2_prop_UT = zeros(6,6,length(t_grid));
P1_prop_UT(:,:,1) = P0;
P2_prop_UT(:,:,1) = P0;

% Define the propagated sigma points and weighted sample mean and covariance
sigma_points_prop1 = sigma_points_01;
sigma_points_prop2 = sigma_points_02;
weight_mean1 = zeros(6,11);
weight_mean2 = zeros(6,11);
weight_mean1(:,1) = x_01';
weight_mean2(:,1) = x_02';
weight_cov1 = P0;
weight_cov2 = P0;


for i = 2:length(t_grid) 
    weight_cov1 = zeros(6,6);
    weight_cov2 = zeros(6,6);
    % Propagate the sigma points and find the weighted sample mean  for MANGO
    for k = 1: 2*n+1
        sigma_points_prop1(:,k) =  keplerian_propagator(t_grid(i-1), sigma_points_prop1(:,k), t_grid(i), mu_E);
        weight_mean1(:,i) = weight_mean1(:,i) + Wm(k)*sigma_points_prop1(:,k);
    end
    % Calculate the weighted sample covariance for MANGO
    for k = 1:2*n+1
        diff1 = sigma_points_prop1(:,k) - weight_mean1(:,i);
        weight_cov1 = weight_cov1 + Wc(k)*(diff1*diff1');     
    end
    P1_prop_UT(:,:,i) = weight_cov1;
    % Propagate the sigma points and find the weighted sample mean  for TANGO
    for k = 1: 2*n+1                     
        sigma_points_prop2(:,k) = keplerian_propagator(t_grid(i-1), sigma_points_prop2(:,k), t_grid(i), mu_E);
        weight_mean2(:,i) = weight_mean2(:,i) + Wm(k)*sigma_points_prop2(:,k);
    end
    % Calculate the weighted sample covariance for MANGO
    for k = 1:2*n+1
        diff2 = sigma_points_prop2(:,k) - weight_mean2(:,i);
        weight_cov2 = weight_cov2 + Wc(k)*(diff2*diff2');
    end
    P2_prop_UT(:,:,i) = weight_cov2;
end

x1_prop_UT = weight_mean1;
x2_prop_UT = weight_mean2;

%% Part 2
% Initialize the the delta, P_sum and the vector saving the otbits index
% for both methods
delta_r = zeros(1,length(t_grid))';
P_sum = zeros(3,3,length(t_grid));
delta_r_UT = zeros(1,length(t_grid))';
P_sum_UT =  zeros(3,3,length(t_grid));
N = zeros(1,length(t_grid));
N_UT = zeros(1,length(t_grid));

% Initialize a vector through it is saved the evolution of the
% triggering equation
 delta_r_lim = zeros(length(t_grid), 1);
 delta_r_lim_UT = zeros(length(t_grid), 1);
% Loop for each iteration in the numeber of orbits
for k = 1 : length(t_grid)
    % Calculate delta r for both methods
    delta_r(k) = norm(x1_prop(1:3,k) - x2_prop(1:3,k));
    delta_r_UT(k) = norm(x1_prop_UT(1:3,k) - x2_prop_UT(1:3,k));
    % Calculate P_sum for both methods
    P_sum(:,:,k) = P1_prop(1:3,1:3,k) + P2_prop(1:3,1:3,k);
    P_sum_UT(:,:,k) = P1_prop_UT(1:3,1:3,k) + P2_prop_UT(1:3,1:3,k);    
    %Calculate delta r for both methods
    delta_r_lim(k) = 3 * sqrt(max(eig(P_sum(:, :, k))));
    delta_r_lim_UT(k) = 3 * sqrt(max(eig(P_sum_UT(:, :, k))));
end
% Check at which orbits is the condition triggered fot LinCov
for k = 1 : length(t_grid)
    if delta_r(k) < 3*sqrt(eigs(P_sum(:,:,k),1)) 
        N(k) = 1;
    else
        N(k) = 0;
    end
end
% Check at which orbits is the condition triggered for UT
for k = 1 : length(t_grid)
    if delta_r_UT(k) < 3*sqrt(eigs(P_sum_UT(:,:,k),1)) 
        N_UT(k) = 1;
    else
        N_UT(k) = 0;
    end
end
% Save the index of the first orbit at which the condition occurs for
% LinCov
if nnz(N ~= zeros(size(N))) ~= 0             % Calculating for the 10 orbits
    firstN = find(N,1,'first') - 1;
else
    firstN = 0;
end
% Save the index of the first orbit at which the condition occurs for UT
if nnz(N_UT ~= zeros(size(N))) ~= 0 
    firstN_UT = find(N_UT,1,'first') - 1;     % Calculating for the 10 orbits
else
    firstN_UT = 0;
end
% Display the results
if firstN > 0
    disp(['The condition through LinCov  occurs at NC = ', num2str(firstN)]);
else
    disp('The condition does not occur for LinCov');
end
if firstN_UT > 0
    disp(['The condition through UT occurs at NC = ', num2str(firstN_UT)]);
else
    disp('The condition does not occur for UT.');
end
if firstN == firstN_UT && firstN > 0
    disp('The condition occurs in the same revolution for both methods');
elseif firstN > 0 
    disp('The condition occurs in different revolutions.');
else 
    disp('The condition did not occur for one or both methods.');
end
% Plot the differences
figure()
hold on
grid on
plot(0:10,delta_r,'DisplayName','$\Delta$r LinCov','LineWidth',2.5)
plot(0:10,delta_r_UT,'--','DisplayName','$\Delta$r UT','LineWidth',2.5)
plot(0:10,delta_r_lim,'g ','DisplayName','$\Delta$r limit LinCov','LineWidth',2.5)
plot(0:10,delta_r_lim_UT,'--','Color','r','DisplayName','$\Delta$r limit UT','LineWidth',2.5)
legend('Location','best','Interpreter','latex')
xlabel('Number of orbits [-]')
ylabel('$\Delta$r [km]','Interpreter','latex')

%% Part III - Monte Carlo 
% Set the number simulations for Monte Carlo
samples = 750;

% Generate random initial data for MANGO
R_1 = mvnrnd(x_01,P0,samples);
% Generate random initial data for TANGO
R_2 = mvnrnd(x_02,P0,samples);
% Initialize the Storing value of Monte Carlo simulation propagation for
% update inside the loop
MC_1 = zeros(6,samples);
MC_2 = zeros(6,samples);

% Initialize mean and covariance for Monte Carlo
P1_prop_MC = zeros(6,6,11);
P2_prop_MC = zeros(6,6,11);
x1_prop_MC = zeros(6,11);
x2_prop_MC = zeros(6,11);

% Define the first mean and covariance instance
P1_prop_MC(:,:,1) = P0;
P2_prop_MC(:,:,1) = P0;

x1_prop_MC(:,1) = x_01';
x2_prop_MC(:,1) = x_02';

% Loop for Monte Carlo propagation
for i = 2:length(t_grid) 
    for k = 1:samples
        [MC_1(:,k),~,~] = keplerian_propagator(t_grid(i-1), R_1(k,:)', t_grid(i), mu_E);
    end
    % Define the mean and covariance in Monte Carlo for MANGO at each time
    % step
    x1_prop_MC(:,i) = mean(MC_1,2);
    P1_prop_MC(:,:,i)  = cov(MC_1');

    % Loop for propagating for each sample for TANGO
    for k = 1:samples
        [MC_2(:,k),~,~] = keplerian_propagator(t_grid(i-1), R_2(k,:)', t_grid(i), mu_E);
    end
    % Define the mean and Covariance in Monte Carlo for TANGO at each time
    % step
    x2_prop_MC(:,i) = mean(MC_2,2);
    P2_prop_MC(:,:,i)  = cov(MC_2');
    % Save the value of propagation of Monte Carlo for update
    R_1 = MC_1';
    R_2 = MC_2';

end


% Display results for Satelllite 1
disp('Monte Carlo simulation Results for Satellite 1: ');
disp('Sample Mean: ');
disp(x1_prop_MC);
disp('Sample Covariance: ');
disp(P1_prop_MC);

% Display results for Satellite 2
disp('Monte Carlo simulation Results for Satellite 2: ')
disp('Sample Mean: ');
disp(x2_prop_MC);
disp('Sample Covariance: ');
disp(P2_prop_MC);


delta_r_MC = zeros(1,length(t_grid))';
P_sum_MC =  zeros(3,3,length(t_grid));
N_MC = zeros(1,length(t_grid));

% Initialize a vector through it is saved the evolution of the
% triggering equation
 delta_r_lim_MC = zeros(length(t_grid), 1);
% Loop for each iteration in the numeber of orbits
for k = 1 : length(t_grid)
    % Calculate delta r for MC
    delta_r_MC(k) = norm(x1_prop_MC(1:3,k) - x2_prop_MC(1:3,k));
    % Calculate P_sum for MC
    P_sum_MC(:,:,k) = P1_prop_MC(1:3,1:3,k) + P2_prop_MC(1:3,1:3,k);    
    %Calculate delta r for MC
    delta_r_lim_MC(k) = 3 * sqrt(max(eig(P_sum_MC(:, :, k))));
end
% Check at which orbits is the condition triggered for MC
for k = 1 : length(t_grid)
    if delta_r_MC(k) < 3*sqrt(max(eig(P_sum_MC(:,:,k)))) 
        N_MC(k) = 1;
    else
        N_MC(k) = 0;
    end
end
figure()
hold on
grid on;
plot(0:10,delta_r_MC,'-','Color','b','DisplayName','$\Delta$r MC','LineWidth',2.5)
plot(0:10,delta_r_lim_MC,'-','Color','c','DisplayName','$\Delta$r limit MC','LineWidth',2.5)
plot(0:10,delta_r_UT,'--','Color','m','DisplayName','$\Delta$r UT','LineWidth',2.5)
plot(0:10,delta_r_lim_UT,'--','Color','r','DisplayName','$\Delta$r limit UT','LineWidth',2.5)
legend('Location','best','Interpreter','latex')
xlabel('Number of orbits [-]')
ylabel('$\Delta$r [km]','Interpreter','latex')
hold off

%% Plot The Results
% Define 
% Initialize arrays to store critical condition
CC_r_1 = zeros(1,length(t_grid))';     % Critical condition in position for spacecraft MANGO
CC_r_UT_1 = zeros(1,length(t_grid))';  % Critical condition in position for spacecraft MANGO (UT)
CC_r_MC_1 = zeros(1,length(t_grid))';  % Critical condition in position for spacecraft MANGO (MC)

CC_r_2 = zeros(1,length(t_grid))';     % Critical condition in position for spacecraft TANGO
CC_r_UT_2 = zeros(1,length(t_grid))';  % Critical condition in position for spacecraft TANGO (UT)
CC_r_MC_2 = zeros(1,length(t_grid))';  % Critical condition in position for spacecraft TANGO (MC)

CC_v_1 = zeros(1,length(t_grid))';     % Critical condition in velocity for spacecraft MANGO
CC_v_UT_1 = zeros(1,length(t_grid))';  % Critical condition in velocity for spacecraft MANGO (UT)
CC_v_MC_1 = zeros(1,length(t_grid))';  % Critical condition in velocity for spacecraft MANGO (MC)

CC_v_2 = zeros(1,length(t_grid))';     % Critical condition in velocity for spacecraft TANGO
CC_v_UT_2 = zeros(1,length(t_grid))';  % Critical condition in velocity for spacecraft TANGO (UT)
CC_v_MC_2 = zeros(1,length(t_grid))';  % Critical condition in velocity for spacecraft TANGO (MC)

% Compute critical condition for each time step
for k = 1 : length(t_grid)
    % Critical condition in position for spacecraft MANGO
    CC_r_1(k)    = 3*sqrt(max(eig(P1_prop(1:3, 1:3,k))));
    CC_r_UT_1(k) = 3*sqrt(max(eig(P1_prop_UT(1:3, 1:3,k))));
    CC_r_MC_1(k) = 3*sqrt(max(eig(P1_prop_MC(1:3, 1:3,k))));

    % Critical condition in velocity for spacecraft MANGO
    CC_v_1(k)    = 3*sqrt(max(eig(P1_prop(4:6, 4:6,k))));
    CC_v_UT_1(k) = 3*sqrt(max(eig(P1_prop_UT(4:6, 4:6,k))));
    CC_v_MC_1(k) = 3*sqrt(max(eig(P1_prop_MC(4:6, 4:6,k))));

    % Critical condition in position for spacecraft TANGO
    CC_r_2(k)    = 3*sqrt(max(eig(P2_prop(1:3, 1:3,k))));
    CC_r_UT_2(k) = 3*sqrt(max(eig(P2_prop_UT(1:3, 1:3,k))));
    CC_r_MC_2(k) = 3*sqrt(max(eig(P2_prop_MC(1:3, 1:3,k))));

    % Critical condition in velocity for spacecraft TANGO
    CC_v_2(k)    = 3*sqrt(max(eig(P2_prop(4:6, 4:6,k))));
    CC_v_UT_2(k) = 3*sqrt(max(eig(P2_prop_UT(4:6, 4:6,k))));
    CC_v_MC_2(k) = 3*sqrt(max(eig(P2_prop_MC(4:6, 4:6,k))));
end

% Plot the results
figure()
hold on
grid on
subplot(2,2,1)
grid on
hold on
plot(0:10, CC_r_1, '-','LineWidth',2.5)
plot(0:10,  CC_r_UT_1,'--','LineWidth',2.5)
plot(0:10,  CC_r_MC_1,'.','LineWidth',2.5,'Color','#00FF00')
title('Mango Position Covariance')
legend('LinCov','UT','MC')
xlabel('Number of orbits')
ylabel('$\sigma_r$','Interpreter','latex')
subplot(2,2,2)
hold on
grid on
plot(0:10, CC_v_1,'-','LineWidth',2.5)
plot(0:10, CC_v_UT_1,'--','LineWidth',2.5)
plot(0:10, CC_v_MC_1,'.','LineWidth',2.5,'Color','#00FF00')
title('Mango Velocity Covariance')
legend('LinCov','UT','MC')
xlabel('Number of orbits [-]')
ylabel('$\sigma_v$','Interpreter','latex')
subplot(2,2,3)
hold on
grid on
plot(0:10, CC_r_2,'-','LineWidth',2.5)
plot(0:10, CC_r_UT_2,'--','LineWidth',2.5)
plot(0:10, CC_r_MC_2,'.','LineWidth',2.5,'Color','#00FF00')
title('Tango Position Covariance')
legend('LinCov','UT','MC')
xlabel('Number of orbits [-]')
ylabel('$\sigma_r$','Interpreter','latex')
subplot(2,2,4)
hold on
grid on
plot(0:10, CC_v_2, '-','LineWidth',2.5)
plot(0:10, CC_v_UT_2,'--','LineWidth',2.5)
plot(0:10, CC_v_MC_2,'.','LineWidth',2.5,'Color','#00FF00')
xlabel('Number of orbits [-]')
ylabel('$\sigma_v$','Interpreter','latex')
title('Tango Velocity Covariance')
legend('LinCov','UT','MC')

orbitplot(r_01,v_01,samples,R_1,x1_prop,x1_prop_MC,x1_prop_UT,P1_prop,P1_prop_MC,P1_prop_UT)
orbitplot(r_02,v_02,samples,R_2,x2_prop,x2_prop_MC,x2_prop_UT,P2_prop,P2_prop_MC,P2_prop_UT)

%% Functions
function [dxdt] = keplerian_rhs(~, x, GM)
%KEPLERIAN_RHS  Evaluates the right-hand-side of a 2-body (keplerian) propagator
%   Evaluates the right-hand-side of a newtonian 2-body propagator.
%
%
% Author
%   Name: ALESSANDRO 
%   Surname: MORSELLI
%   Research group: DART
%   Department: DAER
%   University: Politecnico di Milano 
%   Creation: 24/10/2021
%   Contact: alessandro.morselli@polimi.it
%   Copyright: (c) 2021 A. Morselli, Politecnico di Milano. 
%                  All rights reserved.
%
%
% Notes:
%   This material was prepared to support the course 'Satellite Guidance
%   and Navigation', AY 2021/2022.
%
%
% Inputs:
%   t   : [ 1, 1] epoch (unused)
%   x   : [6, 1] cartesian state vector wrt Solar-System-Barycentre and
%                 State Transition Matrix elements
%   GM  : [ 1, 1] gravitational constant of the body
%
% Outputs:
%   dxdt   : [6,1] RHS, newtonian gravitational acceleration only
%

% Initialize right-hand-side
dxdt = zeros(6,1);

% Extract positions
rr = x(1:3);

% Compute square distance and distance
dist2 = dot(rr, rr);
dist = sqrt(dist2);

% Position detivative is object's velocity
dxdt(1:3) = x(4:6);   
% Compute the gravitational acceleration using Newton's law
dxdt(4:6) = - GM * rr /(dist*dist2);

end

function [xf, tt, xx] = keplerian_propagator(et0, x0, etf, attractor)
%KEPLERIAN_PROPAGATOR Propagate a Keplerian Orbit and the STM
 

% Initialize propagation data
if isfloat(attractor)
    GM = attractor;
else
    GM = cspice_bodvrd(attractor, 'GM', 1);
end

tof = etf-et0;

options = odeset('reltol', 1e-12, 'abstol', [ones(3,1).*1e-8; ones(3,1).*1e-11]);

% Perform integration
[tt, xx] = ode78(@(t,x) keplerian_rhs(t,x,GM), [0 tof], x0, options);

% Extract state vector 
xf = xx(end,1:6)';

end

function [xf, tt, xx,PHIf] = keplerian_propagator_STM(et0, x0, etf, attractor)
%KEPLERIAN_PROPAGATOR Propagate a Keplerian Orbit and the STM
%   Detailed explanation goes here

% Initialize propagation data
if isfloat(attractor)
    GM = attractor;
else
    GM = cspice_bodvrd(attractor, 'GM', 1);
end

Phi = eye(6);
x0Phi = [x0;Phi(:)];
options = odeset('reltol', 1e-12, 'abstol', 1e-13);

% Perform integration
[tt, xx] = ode78(@(t,x) keplerian_STM_rhs(t,x,GM), [et0 etf], x0Phi, options);

% Extract state vector 
xf = xx(end,1:6)';
PHIf = reshape(xx(end,7:end),6,6);

end
 function [dxdt] = keplerian_STM_rhs(t, x, GM)
%KEPLERIAN_RHS  Evaluates the right-hand-side of a 2-body (keplerian)
%               propagator with STM
%   Evaluates the right-hand-side of a newtonian 2-body propagator with STM.
%
%
% Author
%   Name: ALESSANDRO 
%   Surname: MORSELLI
%   Research group: DART
%   Department: DAER
%   University: Politecnico di Milano 
%   Creation: 11/10/2023
%   Contact: alessandro.morselli@polimi.it
%   Copyright: (c) 2023 A. Morselli, Politecnico di Milano. 
%                  All rights reserved.
%
%
% Notes:
%   This material was prepared to support the course 'Satellite Guidance
%   and Navigation', AY 2023/2024.
%
%
% Inputs:
%   t   : [ 1, 1] epoch (unused)
%   x   : [42, 1] cartesian state vector wrt Solar-System-Barycentre and
%                 State Transition Matrix elements
%   GM  : [ 1, 1] gravitational constant of the body
%
% Outputs:
%   dxdt   : [42,1] RHS, newtonian gravitational acceleration only
%

% Initialize right-hand-side
dxdt = zeros(42,1);

% Extract positions
rr = x(1:3);
Phi = reshape(x(7:end),6,6);

% Compute square distance and distance
r_norm = norm(rr);
dist2 = dot(rr, rr);
dist = sqrt(dist2);

% Compute the gravitational acceleration using Newton's law
aa_grav =  - GM * rr /r_norm^3;

% Compute the derivative of the flow
dfdv = 3*GM/dist^5*(rr*rr')-GM/dist^3*eye(3);

% Assemble the matrix A(t)=dfdx
dfdx = [zeros(3), eye(3); dfdv, zeros(3)];
% Compute the derivative of the state transition matrix
Phidot = dfdx*Phi;

dxdt(1:3) = x(4:6);   % Position detivative is object's velocity
dxdt(4:6) = aa_grav;  % Sum up acceleration to right-hand-side
dxdt(7:end) = Phidot(:);

 end

function orbitplot(r, v, samples, R, x, x_MC, x_UT, P, P_MC, P_UT)
% ORBITPLOT Plot the orbit and covariance ellipses
% Inputs:
%   r       : Initial position vector
%   v       : Initial velocity vector
%   samples : Number of samples
%   R       : Rotation matrices
%   x       : State vector
%   x_MC    : Monte Carlo state vector
%   x_UT    : Unscented transform state vector
%   P       : Linear covariance matrix
%   P_MC    : Monte Carlo covariance matrix
%   P_UT    : Unscented transform covariance matrix


n = 500;
angle = 2 * pi / n * (0:n);
circle = [cos(angle); sin(angle)];

% Calculate local frame
i = r / norm(r);
k = cross(r, v) / norm(cross(r, v));
j = cross(k, i);
R_frame = [i; j; k];

% Transform positions to local frame
Pos_new = zeros(3, samples);
for i = 1:samples
    Pos_new(:, i) = R_frame * R(i, 1:3)';
end

% Rotate mean and covariance estimates
x_rotated_MC = R_frame * x_MC(1:3, end);
P_rotated_MC = R_frame * P_MC(1:3, 1:3, end) * R_frame';
x_rotated    = R_frame * x(1:3, end);
P_rotated    = R_frame * P(1:3, 1:3, end) * R_frame';
x_rotated_UT = R_frame * x_UT(1:3, end);
P_rotated_UT = R_frame * P_UT(1:3, 1:3, end) * R_frame';
P_rot(:, :, 1) = P_rotated(1:2, 1:2);
P_rot(:, :, 2) = P_rotated_UT(1:2, 1:2);
P_rot(:, :, 3) = P_rotated_MC(1:2, 1:2);
x_rot = [x_rotated(1:2), x_rotated_UT(1:2), x_rotated_MC(1:2)];

% Plot orbit and covariance ellipses
figure()
hold on
grid on
plot(Pos_new(1, :), Pos_new(2, :), 'k.');
for i = 1:3
    [R, D] = svd(P_rot(:, :, i));
    ellipse = 3 * R * sqrt(D) * circle;
    x = x_rot(1, i) + ellipse(1, :);
    y = x_rot(2, i) + ellipse(2, :);
    plot(x, y,'LineWidth',1.2)
end
% Plot points
plot(x_rotated(1), x_rotated(2), 'diamond', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')
plot(x_rotated_UT(1), x_rotated_UT(2), '*', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g')
plot(x_rotated_MC(1), x_rotated_MC(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'none','MarkerEdgeColor','c')
xlabel('x [km]')
ylabel('y [km]')
legend('Propagated points', 'LinCov ellipse', 'UT ellipse', 'MC ellipse', ...
    'LinCov mean', 'UT mean', 'MC mean')
end

