% Spacecraft Guidance and Navigation (2023/2024)
% Assignment 1 Exercise 2
% Martin Mihnea Stefan 10903720
clc; clearvars; close all;cspice_kclear;
rng default;

 %% Initialize
% Kernels
addpath('.\kernels')
addpath('.\mice\src\mice')
addpath('.\mice\lib')

cspice_furnsh('kernels\de432s.bsp');
cspice_furnsh('kernels\gm_de432.tpc');
cspice_furnsh('kernels\naif0012.tls');
cspice_furnsh('kernels\pck00010.tpc');
cspice_furnsh('kernels\20099942_Apophis.bsp');
cspice_furnsh('ex02.tm');

% Part one: initialize data for the propagator
% Define list of celestial bodies:
labels = {'Sun';
          'Mercury';
          'Venus';
          'Earth';
          'Moon';
          'Mars Barycenter';
          'Jupiter Barycenter';
          'Saturn Barycenter';
          'Uranus Barycenter';
          'Neptune Barycenter';
          'Pluto Barycenter'};
components = {'Sun','Earth','Moon'};
%Initialize propagation data (same as regular n-body)
bodies = nbody_init(labels);
mu = bodies{1}.GM;

%Time Epoch
et0 = cspice_str2et('2029-Jan-01 00:00:00.0000 TDB');
etf = cspice_str2et('2029-Jul-31 23:59:59.9999 TDB');
time = linspace(et0,etf,1e+6);
time_plot =  cspice_et2utc( time, 'C',0);
date_time = datetime(time_plot, 'InputFormat', 'yyyy MMM dd HH:mm:ss');
%Take the data from kernels
EARTH = cspice_spkpos('20099942', time,'IAU_EARTH','NONE','EARTH');
MOON = cspice_spkpos('20099942', time,'IAU_EARTH','NONE','MOON');
SUN = cspice_spkpos('20099942', time,'IAU_EARTH','NONE','SUN');
EARTH_SUN = cspice_spkpos('SUN',time,'IAU_EARTH','NONE','EARTH');
r_earth = vecnorm(EARTH,2,1);
r_moon = vecnorm(MOON,2,1);
r_sun = vecnorm(SUN,2,1);
r_sun_earth = vecnorm(EARTH_SUN,2,1);

% Plot the Distances
axis tight;
figure(1)
ax = NaN(1,3);
dist = [r_earth;r_moon;r_sun];
for i=1:3
    ax(i) = subplot(3,1,i);
    plot(time/cspice_spd ,dist(i,:));
    xlabel('Time since Epoch [MJD2000]','FontSize',11,'FontName','Arial','FontWeight','bold')
    ylabel('Distance','Interpreter','tex','FontSize',11,'FontName','Arial','FontWeight','bold')
    title(components{i},'Interpreter','tex')
    axis tight;

end
%linkaxes(ax,'x');
%Calculate the angle Earth-Apophis-Sun
distSE=ones(1,length(time))*149597870700/1000;
alpha=acos((r_earth.^2+r_sun.^2-r_sun_earth.^2)./(2*r_earth.*r_sun))*180/pi;
% Plot the Angle
figure(2)
hold on;
plot(date_time,alpha,'LineWidth',3)
hold off;
grid on;
xlabel('Time since Epoch [MJD2000]', 'FontSize',14,'FontName','Arial','FontWeight','bold')
ylabel('Angle [deg]','FontSize',14,'FontName','Arial','FontWeight','bold')       
title('Earth-Apopsis-Sun Angle','FontName','Arial','FontSize',14)
components = {'Earth' 'Moon' 'Sun'};
set(gca, 'FontSize', 8);
axis tight;
%Obtain the minimum positionn between Earth and Apophis and its index
[dmin,minpos]=min(r_earth);
%Calulate the time when the position occurs
tmin = time(minpos);

%Define the time window
time_window = tmin - 3600*6:tmin + 3600*6;
%Extract position of Apophis with Earth in the time window from kernels
rr_Apophis_spkpos = cspice_spkpos('20099942',time_window,'IAU_EARTH','NONE','EARTH');
%Extract the longitude and the latitude from kernels at the obtained
%position
[radapo, lonapo, latapo] = cspice_reclat(rr_Apophis_spkpos);
lonapo = lonapo*cspice_dpr;
latapo = latapo*cspice_dpr;
[~,minlonapo] = min(lonapo);
b = time_window - tmin;
[~,index] = min(abs(time_window - tmin));
%Plot the ground Track
figure(3)
hold on;
c = imread('2k_earth_daymap.jpg');%importing the image into the plot
image(c,'XData',[-180 180],'YData',[90 -90]);
plot(lonapo(1:minlonapo),latapo(1:minlonapo),'-','Color','#D95319',"LineWidth",4)
plot(lonapo(minlonapo+1:end),latapo(minlonapo+1:end),'-','Color','#D95319',"LineWidth",4)
plot(lonapo(1),latapo(1),'o','MarkerSize',10,'MarkerFaceColor',[0.6113 0.115, 0.765]);
plot(lonapo(end),latapo(end),'s','MarkerSize',10,'MarkerFaceColor',[0.0113 0.415, 0.765]);
plot(lonapo(21602),latapo(21602),'o','MarkerSize',10,'MarkerFaceColor',[0.111 0.744 0.348]);
axis([-180 180 -90 90]);
xlabel('Longitude');
ylabel('Latitude');
legend('','GroundTrack','Starting Point','Ending Point','Minimum Distance','Location','best')
%% Part 2

clearvars -except tmin; cspice_kclear;
cspice_furnsh('kernels\de432s.bsp');
cspice_furnsh('kernels\20099942_Apophis.bsp');
cspice_furnsh('kernels\gm_de432.tpc');
cspice_furnsh('kernels\pck00010.tpc');
rng default;

labels = {'Sun';
          'Mercury';
          'Venus';
          'Earth';
          'Moon';
          'Mars Barycenter';
          'Jupiter Barycenter';
          'Saturn Barycenter';
          'Uranus Barycenter';
          'Neptune Barycenter';
          'Pluto Barycenter'};

%Initialize propagation data (same as regular n-body)
bodies = nbody_init(labels);
mu = bodies{1}.GM;


% Load the kernel
cspice_furnsh('ex02.tm');
conv = 86400;
%Initial dates for the Time windows
ref_epoch_str_imp2 = '2028-Aug-01 00:00:00.0000 TDB'; % Impact Window Open
et_impact_WO = cspice_str2et(ref_epoch_str_imp2);
ref_epoch_fin_imp2 = '2029-Feb-28 23:59:59.0000 TDB'; % Impact Window Close
et_impact_WC = cspice_str2et(ref_epoch_fin_imp2);
ref_epoch_LWO = '2024-Oct-01 00:00:00.0000 TDB';      % Launch Window Open
et_LWO = cspice_str2et(ref_epoch_LWO);
ref_epoch_LWC = '2025-Feb-01 23:59:59.0000 TDB';      % Launch Window Close
et_LWC = cspice_str2et(ref_epoch_LWC);
et_DSM_WO = cspice_str2et('2025-Apr-01 00:00:00.0000 TDB');
et_DSM_WC = cspice_str2et('2026-Aug-01 23:59:59.0000 TDB');
%Initial guesses for time; Different trials have been used
guess_departure = et_LWO  + 1/3 * (et_LWC - et_LWO); 
guess_DSM = et_DSM_WO  + 1/3 * (et_DSM_WC - et_DSM_WO);    
guess_impact = et_impact_WO   + 1/3 * (et_impact_WC - et_impact_WO);

%Define the upper and lower boundary for kernels
e_min = -1e+21*ones(6,1);           %values low enough that the Apophis would not reach
e_plus = 1e+21*ones(6,1);           %values high enough that the Apophis would not reach                
lb = [e_min; e_min; e_min; et_LWO;et_DSM_WO;et_impact_WO];
ub = [e_plus;e_plus;e_plus;et_LWC;et_DSM_WC;et_impact_WC];


% Retrieve Earth state vector at the guessed departure time
x_EARTH = cspice_spkezr('EARTH', guess_departure, 'ECLIPJ2000', 'NONE', 'SSB');
%Define the initial deltav implied in each manouver
 rand_dv = ceil(abs(randn()) * 2.5);
 rand1 = abs(randn(3,1));
 DVL = rand_dv * rand1/norm(rand1);
 rand2 = abs(randn(3,1));
 DSM = (5 - rand_dv) * rand2/norm(rand2);


% Propagate states using keplerian propagator and add the deltavs
x1 = x_EARTH + [zeros(3,1); DVL];
[x12, ~,~] = keplerian_propagator(guess_departure, x1, guess_DSM, mu);
x2 = x12   + [zeros(3,1); DSM];
[x3, ~, ~] = keplerian_propagator(guess_DSM, x2, guess_impact, mu);
% Initial guess for the optimization
xx0 = [x1; x2; x3; guess_departure ; guess_DSM ; guess_impact ];


% Optimization options
optfmincon = optimoptions('fmincon','Display','iter');
optfmincon.MaxIterations = 150;
optfmincon.StepTolerance = 1e-12;
optfmincon.ConstraintTolerance = 1e-12;
optfmincon.FunctionTolerance = 1e-6;

%Define the objective function for the fmincon
objective = @(y) SCobjective(y,bodies(1),mu);
% Perform optimization using fmincon
[SOL_FMINCON,FVAL] = fmincon(objective,xx0,[],[],[],[],lb,ub,@(y) constraints(y,bodies(1),mu),optfmincon);
%Retrieve the data obtained from the fmincon
[DISTANCE_E_A,TIME_E_A,DISTANCE_E_A_xx,tt_sol] = SCobjective(SOL_FMINCON(1:21),bodies(1));

% Calculate deltav launch
x_earth_sol = cspice_spkezr('EARTH', SOL_FMINCON(19) ,'ECLIPJ2000','NONE','SSB');
deltav_launch_sol = SOL_FMINCON(4:6) - x_earth_sol(4:6); 
dv_launch = norm(deltav_launch_sol);

% Calculation for deltav dsm
[xf12, ~, ~] = keplerian_propagator(SOL_FMINCON(19),SOL_FMINCON(1:6),SOL_FMINCON(20), mu);
deltav_dsm_sol = SOL_FMINCON(10:12) - xf12(4:6);
dv_dsm = norm(deltav_dsm_sol);
dv = dv_dsm + dv_launch;

fprintf('Launch [UTC]: %s\n', cspice_et2utc(SOL_FMINCON(19) ,'C',3));
fprintf('DSM  [UTC]:%s\n',cspice_et2utc(SOL_FMINCON(20) ,'C',3));
fprintf('Impact [UTC]: %s\n', cspice_et2utc(SOL_FMINCON(21) ,'C',3));
fprintf('TCA [UTC]: %s\n', cspice_et2utc(TIME_E_A,'C',3));
fprintf('Delta launch [km/s] %.2f\n', deltav_launch_sol);
fprintf('Delta DSM [km/s] %.2f\n', deltav_dsm_sol);
fprintf('DCA[Re] %.2f\n', DISTANCE_E_A);

%% Ploting

Re = cspice_bodvrd('399','RADII',3);
time_vec = linspace(SOL_FMINCON(19) , tmin, 500);
dreal = vecnorm(cspice_spkpos('20099942',time_vec,'IAU_EARTH','NONE','EARTH'),2,1)./Re(1);

time_vec_plot = linspace(time_vec(1),time_vec(end),6);
time_plot =  cspice_et2utc(time_vec_plot, 'C',0);
position = linspace(1,length(time_vec),10);

% Plot the data; distance of apophis wrt to earth in Earth RADII
figure()
axis tight;
hold on
grid on
plot(time_vec,dreal,'-','LineWidth',2,'Color','#0072BD')
hold on
plot(tt_sol,-DISTANCE_E_A_xx,'r--','LineWidth',2)
xticks(time_vec_plot);
xticklabels(time_plot(:,1:11));
ylabel('DistanceEarth Apophis')
hold off;

firstxlim = tt_sol(-DISTANCE_E_A_xx < 40);
time_vec_plot = linspace(firstxlim(1),tt_sol(end),6);
time_plot =  cspice_et2utc(time_vec_plot, 'C',0);

figure()
hold on;
grid on;
plot(time_vec,dreal,'-','LineWidth',2,'Color','#0072BD')
hold on
plot(tt_sol,-DISTANCE_E_A_xx,'r--','LineWidth',2)
xticks(time_vec_plot);
xticklabels(time_plot(:,6:17));
ylabel('DistanceEarth Apophis')

xlim([firstxlim(1) tt_sol(end)]);
ylim([5 40]);
hold off;
%% Functions
function [DISTANCE_E_A,TIME_E_A,DISTANCE_E_A_xx, tt_apo] = SCobjective(y,bodies,~)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function SCobjective is the objective function of the spacecraft,
% that is needed in fmincon for trajectory optimization 
% Inputs:
%        - x[21x1] is the vector of states and times: it includes the
% values that need to be optimized
%        - bodies [cell array] is the cell array containing informationn on the celestial bodies
%        - ~ - is the GM,graviational parameter, can be introduced if needed depending on the type
%        of propagation
% Outputs: 
%           - DISTANCE_E_A[1X1] is the optimal value obtained by the fmincon
%           - TIME_E_A [1x1] is the time values respective for the DISTANCE_E_A
%           - DISTANCE_E_A_xx[nx1] is the vector of points on the Apophis trajectory,
%           resulted from the initial value to the optimal value obtained.
%           - tt_apo[nx1] the vector containing the times respective to each value in
%           vector DISTANCE_E_A_xx

    Re = cspice_bodvrd('Earth','RADII',3);
    xapo = cspice_spkezr('20099942', y(21), 'ECLIPJ2000','NONE','SSB') + [zeros(3,1);y(16:18)*5e-5];   
    if nargout < 4
        [x_apo,~,~,tt_apo]= propagate(y(21) ,xapo,y(21)*10,bodies,'ECLIPJ2000','True');        
        x_earth = cspice_spkpos('EARTH', tt_apo(end), 'ECLIPJ2000','NONE','SSB');
        % Calculate the distance between earth and apophis for fmincon
        DISTANCE_E_A = -norm(x_earth(1:3) - x_apo(1:3))/Re(1);
    else
        %Propagate the data of apophis after the impact
        [x_apo,~,xx_apo,tt_apo] = propagate_tf(y(21) ,xapo,y(21)*10,bodies,'ECLIPJ2000','True');
        x_earth = cspice_spkpos('EARTH', tt_apo(end), 'ECLIPJ2000','NONE','SSB');   
        %Calculate the final distance between Eath and Apophis
        DISTANCE_E_A = -norm(x_earth(1:3) - x_apo(1:3))/Re(1);
        %Take the vector of distance evolution with time for plot 
        DISTANCE_E_A_xx =  -vecnorm(cspice_spkpos('EARTH', tt_apo', 'ECLIPJ2000','NONE','SSB') - xx_apo(:,1:3)')./Re(1);
        %Take the time vector for plot
       TIME_E_A = tt_apo(end);

    end

end


 function [CC, CC_eq] = constraints(xx,~,mu)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that inputs the nonlinear constraintd for the fmincon, to obtain
% the trajectory optimization
% Inputs: - xx[21,1]: Vector containg the optimization parameters that are
%         placed under constraints
%         - bodies: bodies [cell array] is the cell array containing information
%                 on the celestial bodies that can be needed based on the
%                 type of propagator
%           - mu[1,1]: GM,graviational parameter, can be introduced if needed 
%             depending on the type of propagation
% Outputs: - CC[1]: Inequality constraint value
%          - CC_eq[15,1]: Equality constraint values.
%Initialize data
t1 = xx(19);
t2 = xx(20);
t3 = xx(21);
x1 = xx(1:6);
x2 = xx(7:12);
x3 = xx(13:18);


%Propagate the initial data
[x12, ~, ~] = keplerian_propagator(t1,x1,t2, mu);
[x23, ~, ~] = keplerian_propagator(t2,x2,t3, mu);
%Retrieve data from kernels
x_earth = cspice_spkezr('EARTH',    t1 ,'ECLIPJ2000','NONE','SSB');
x_apo   = cspice_spkezr('20099942', t3 ,'ECLIPJ2000','NONE','SSB');
%Calculate the deltavs that are placed under non linea constraints
deltaV1 = x1(4:6) - x_earth(4:6);
deltaV2 = x12(4:6) - x2(4:6);

% Define the equality constraint values vector
CC_eq = [x1(1:3) - x_earth(1:3);
         x12(1:3) - x2(1:3);
         x23(1:6) - x3(1:6);
         -x_apo(1:3) + x3(1:3); 
         ];
%Define the inequality constraint vectors
CC = norm(deltaV1) + norm(deltaV2) - 5;
end

function [bodies] = nbody_init(labels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NBODY_INIT Initialize planetary data for n-body propagation
%   Given a set of labels of planets and/or barycentres, returns a
%   cell array populated with structures containing the body label and the
%   associated gravitational constant.
% Inputs:
%   labels : [1,n] cell-array with object labels
%
% Outputs:
%   bodies : [1,n] cell-array with struct elements containing the following
%                  fields
%                  |
%                  |--bodies{i}.name -> body label
%                  |--bodies{i}.GM   -> gravitational constant [km**3/s**2]


%Initialize output 
bodies= cell(size(labels));

%Loop over labels
for i=1:length(labels)
    %Store body labels
    bodies{i}.name=labels{i};
    %Store body gravitational constant
    bodies{i}.GM = cspice_bodvrd(labels{i}, 'GM' , 1);
end

end
% Function nbody RHS
function [dx] = nbody_rhs(t, x, bodies, frame)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the right-hand side of the ordinary differential 
% equation (ODE) for the n-body problem.
%
% INPUTS:
% - t[1]: Time variable.
% - x[6x1]: State vector containing position and velocity components.
% - bodies[cell array]: Cell array containing information about celestial 
%                      bodies (name, gravitational parameter GM).
% - frame[char]: Reference frame for computing positions (e.g., 'ECLIPJ2000').
%
% OUTPUT:
% - dx[6x1]: Derivative of the state vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize right-hand-side
dx = zeros(6,1);
% Extract the object position from state x
rr_ssb_obj = x(1:3);
%r = norm(rr);
% Loop over all bodies
for i = 1:length(bodies)

    % Retrieve position and velocity of i-th celestial body wrt Solar
    % System Barycentre in inertial frame
    %bodies{i}.name
    rv_ssb_body = cspice_spkezr(bodies{i}.name,t,frame,'NONE','SSB');
    % Extract object position wrt. i-th celestial body
    rr_body_obj = rr_ssb_obj - rv_ssb_body(1:3);
    % Compute square distance and distance
    %r_i = norm(rr_i);
    dist2 = dot(rr_body_obj, rr_body_obj);
    dist = sqrt(dist2);
    % Compute the gravitational acceleration using Newton's law
    % GM = bodies{i}.GM;
    aa_grav = - bodies{i}.GM*rr_body_obj / (dist*dist2);

    % Position derivative is object's velocity
    dx(1:3) = x(4:6);
    % Sum up acceleration to right-hand-side
    dx(4:6) = dx(4:6) + aa_grav;

end
end

function [xf, tt, xx] = keplerian_propagator(et0, x0, etf, attractor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%KEPLERIAN_PROPAGATOR Propagate a Keplerian Orbit and the STM
% Inputs:  et0: the initial epoch of the propagation [1x1]
%          etf: the final epoch of the propagation [1x1]
%          x0: the initial states to be propagates [6x1]
%          attractor: contains the data for the celestial bodies that
%          influences the trajectory of the spacecraft
% Output: xf: the propagated states at final time step [6x1]
%         xx: the matrix containing the states at each time step
%         tt: the vector containing all time steps 

% Initialize propagation data
if isfloat(attractor)
    GM = attractor;
else
    GM = cspice_bodvrd(attractor, 'GM', 1);
end

tof = etf-et0;

options = odeset('reltol', 1e-11, 'abstol', 1e-12);

% Perform integration
[tt, xx] = ode78(@(t,x) keplerian_rhs(t,x,GM), [0 tof], x0, options);

% Extract state vector 
xf = xx(end,1:6)';

end

% Propagete function
function [xf,tf,xx,tt] = propagate(et0,x0,etf,bodies,frame,varargin) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function propagates the states 
% Inputs: et0: the initial epoch of the propagation [1x1]
%          etf: the final epoch of the propagation [1x1]
%          x0: the initial states to be propagates [6x1]
%          bodies: the data of the celestial bodies which may affect the propagation [cell]
%          frame: the inertial frame [char]
%          varargin: it contains the choice to use or not the Event Function in the propagator [1x1]
% Output: xf: the propagated states at final time step [6x1]
%         tf: the final time step [1x1]
%         xx: the matrix containing the states at each time step
%         tt: the vector containing all time steps 

    evtFlag = varargin{1};
    % Perform integration
    options = odeset('RelTol', 1e-8,'AbsTol',1e-9,'Events',@(t,s) EventFct(t,s,evtFlag));
    [tt,xx] = ode78(@(t,x) nbody_rhs(t,x,bodies,frame),[et0 etf], x0, options);

    % extract state and State Transition Matrix
    xf = xx(end,1:6)';

    tf = tt(end);
end

function [xf,tf,xx,tt] = propagate_tf(et0,x0,etf,bodies,frame,varargin) 
    
    evtFlag = varargin{1};
    time_vec = linspace(et0,etf,10^7);
    
    % Perform integration
   options = odeset('RelTol', 1e-8,'AbsTol',1e-9,'Events',@(t,s) EventFct(t,s,evtFlag));
%     options = odeset('RelTol', 1e-8,'AbsTol',1e-9);
   [tt,xx] = ode78(@(t,x) nbody_rhs(t,x,bodies,frame),time_vec, x0, options);
   
    % extract state and State Transition Matrix
    xf = xx(end,1:6)';

    tf = tt(end);
end

function [value,isterminal,direction] = EventFct(t,xx,isTerminal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The event function needed to stop the integration at minimum deltar
% Inputs:    t: integration time 
%           xx: the state vector
%           isTerminal: value ordering the termination of the integration(1 or 0)
% Outputs: value: the variable that starts the propagation
%          isTerminal: 1
%          direction:always 1
     x_EARTH = cspice_spkezr('EARTH',t,'ECLIPJ2000','NONE','SSB');
     value = dot(x_EARTH(1:3)-xx(1:3), x_EARTH(4:6) - xx(4:6));
     isterminal = isTerminal;
     direction  = 1;
end


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
