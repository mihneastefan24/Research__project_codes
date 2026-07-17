%Spacecraft Guidance and Navigation 2023/2024
% Assignment 2 Exercie 3
% Martin Mihnea Stefan

clearvars; close all; clc;
rng default
%% Initialiaze

addpath('.\kernels\')
addpath('.\sgp4\')
addpath('.\tle\')
addpath('.\mice\lib')
addpath('.\mice\src\mice')

cspice_furnsh('assignment02.tm');
%% Exercise 1

% Satellite 1's mean state
r_MAN = [4622.232026629, 5399.3369588058, -0.0212138165769957];
v_MAN = [0.812221125483763, -0.721512914578826, 7.42665302729053];

% Satellite 2's mean state
r_TAN = [4621.69343340281, 5399.26386352847, -3.09039248714313];
v_TAN = [0.813960847513811, -0.719449862738607, 7.42706066911294];

x_TAN = [r_TAN,v_TAN]';
x_MAN = [r_MAN,v_MAN]';

% Covariance P0

P0 = [5.6e-7 3.5e-7 -7.1e-8 0 0 0;
      3.5e-7 9.7e-7 7.6e-8  0 0 0; 
      -7.1e-8 7.6e-8 8.1e-8  0 0 0;
      0 0 0 2.8e-11 0 0;
      0 0 0 0 2.7e-11 0;
      0 0 0 0 0 9.6e-12];

tsep = cspice_str2et('August 12 05:27:39.114 UTC 2010');
t0   = cspice_str2et('2010-08-12T05:30:00.000');
tf   = cspice_str2et('2010-08-12T06:30:00.000');

t_grid = [t0:5:tf];

mu_E = cspice_bodvrd('EARTH','GM',1);
J2 = 0.0010826269;
Re = cspice_bodvrd('Earth','RADII',3);
SVAL = struct('name','SVALBARD','LAT',78.229772,'LONG',15.407786,'ALT',  458 , 'minel',5, 'TYPE','RADAR', 'R', diag([0.01^2 0.125^2 0.125^2]), 'TOPO','SVALBARD_TOPO', 'W', diag([1/0.01 1/0.125 1/0.125]));
KOUROU = struct('name','KOUROU','LAT',5.25144,  'LONG',-52.80466,'ALT' ,-14.67,'minel',10,'TYPE','RADAR', 'R', diag([0.01^2 0.1^2 0.1^2]),      'TOPO','KOUROU_TOPO',   'W', diag([1/0.01 1/0.1 1/0.1]));

MANGO = struct('name', 'MANGO','ID', 36599, 'Mango',x_MAN, 'COV',P0,'TLE1','1 36599U 10028B   10224.22752732 -.00000576  00000-0 -16475-3 0  9998', ...
                                                                    'TLE2','2 36599 098.2803 049.5758 0043871 021.7908 338.5082 14.40871350  8293');
TANGO = struct('name', 'TANGO','ID', 36827, 'Tango',x_TAN, 'COV',P0,'TLE1','1 36827U 10028F   10224.22753605  .00278492  00000-0  82287-1 0  9996', ...
                                                                    'TLE2','2 36827 098.2797 049.5751 0044602 022.4408 337.8871 14.40890217    55');

[xf,~,~,~] = keplerian_propagator(tsep, x_MAN, t0, mu_E);
x_MAN_span(:,1) = xf; 
for k = 2:length(t_grid)
    [xf,~,~,~] = keplerian_propagator(t_grid(k-1), x_MAN_span(:,k-1),t_grid(k),mu_E);
    x_MAN_span(:,k) = xf;
end
 

[SVAL.ECI,SVAL.ECI2TOPO] = stationcoordinates(SVAL,t_grid);
SVAL.Mango_coord = SCGScoord(SVAL, x_MAN_span);
SVAL.vis_time = t_grid(SVAL.Mango_coord(3,:) >= SVAL.minel );
SVAL.vis_position = SVAL.Mango_coord(:,SVAL.Mango_coord(3,:) >= SVAL.minel);

% b) Simultate Measurements Mango Only

[MANGO.coord,r_ECI,v_ECI,et0,SVAL, vis_time,measures] = measurements(MANGO, SVAL,SVAL.vis_time);

x_MAN_ECI = [r_ECI;v_ECI];
MANGO.sgp4 = x_MAN_ECI(:,logical(MANGO.coord (end,:)));

% c)From 1st exercise
MANGO.measure = [[t0,0,0,0];[vis_time,measures]];

%Initialize variables
P_UKF = zeros(6,6,length(MANGO.measure(:,1)));
x_UKF = zeros(6,length(MANGO.measure(:,1)));
sigma_points = zeros(6,13,length(MANGO.measure(:,1)));
sigma_r_UKF = zeros(length(MANGO.measure(:,1))-1,1);
sigma_v_UKF = zeros(length(MANGO.measure(:,1))-1,1);
%Set the initial values of the variable at t0
P0 = 1e4*P0;            
SVAL.ECI_UKF = SVAL.ECI(:,1);
SVAL.ECI2TOPO_UKF = SVAL.ECI2TOPO(:,:,1);
valueoftruth = 0;
% Propagate the initial value with J2
[sigma_points(:,:,1),x_UKF(:,1),P_UKF(:,:,1)] = UKF(x_MAN,P0,SVAL,MANGO.measure,J2,[tsep t0],valueoftruth);
                                                
% initialize performance parameters
sigma_r = zeros(length(MANGO.measure(:,1)) - 1,1);
sigma_v = zeros(length(MANGO.measure(:,1)) - 1,1);
% Sequential filter 
valueoftruth = 1;
for k = 2:length(MANGO.measure(:,1))
    SVAL.ECI_UKF = SVAL.ECI(:,k-1);
    SVAL.ECI2TOPO_UKF = SVAL.ECI2TOPO(:,:,k-1); 
    [sigma_points(:,:,k),~,~,x_UKF(:,k),P_UKF(:,:,k)] =  UKF(x_UKF(:,k-1),P_UKF(:,:,k-1),SVAL,MANGO.measure(k,2:4) ,J2,[MANGO.measure(k-1,1) MANGO.measure(k,1)],valueoftruth); 
end
for i = 2 : length(MANGO.measure(:,1))
     sigma_r(i - 1) = 3*sqrt(trace(P_UKF(1:3,1:3,i)));
     sigma_v(i - 1) = 3*sqrt(trace(P_UKF(4:6,4:6,i)));
end
% Plot
time_plot = datetime(cspice_timout(MANGO.measure(2:end,1)',...
    'YYYY-MM-DD HR:MN:SC.###'),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');

error_r = vecnorm(MANGO.sgp4(1:3,:) - x_UKF(1:3,2:end),2,1);
error_v = vecnorm(MANGO.sgp4(4:6,:) - x_UKF(4:6,2:end),2,1);
figure()
subplot(1,2,1)
semilogy(time_plot,error_r,'LineWidth',1.5)
hold on
semilogy(time_plot,sigma_r,'LineWidth',1.5)
grid on
xlabel('Time')
ylabel('Range [km]')
legend('error','3$\sigma$','Interpreter','latex')

subplot(1,2,2)
semilogy(time_plot,error_v,'LineWidth',1.5)
hold on
semilogy(time_plot,sigma_v,'LineWidth',1.5)
grid on
ylabel('Velocity [km/s]')
xlabel('Time')
legend('error','3$\sigma$','Interpreter','latex')


%% Exercise 2

% a)
FFRF =  diag([(1e-5)^2 1^2 1^2]);
Mango.TLE{1} = MANGO.TLE1;
Mango.TLE{2} = MANGO.TLE2;  
%Calculate the position of Mango in ECI
[~,r_ECI_M,v_ECI_M,~,~] = measurements(MANGO, SVAL,t0);     

n_CW = sqrt(mu_E/norm(r_ECI_M)^3);
%Calculate the transformation matrix from ECI to LVLH
R = LVLHTrans(r_ECI_M,v_ECI_M,n_CW);                       
% Transform the Mango states in LVLH from ECI
x_LVLH_M = R * [r_ECI_M;v_ECI_M];
[~,r_ECI_T,v_ECI_T,~,~] = measurements(TANGO, SVAL,t0);
% Transform the relative data between the two satellites in the new reference frame
MAN_TAN_state = -R*[ (r_ECI_M - r_ECI_T); (v_ECI_M - v_ECI_T)];

% b) 

x_CW = zeros(6,length(t_grid));
x_CW(:,1) = MAN_TAN_state;
% Use the Clohessy-Wiltshire (CW) equations to propagate the states 
for i = 2:length(t_grid)
      [~,x_CW(:,i)] = CW(x_CW(:,i-1),t_grid(i-1),t_grid(i) ,n_CW);
end 

sat_coord = cspice_xfmsta(x_CW, 'RECTANGULAR', 'LATITUDINAL', 'Earth');
Az = wrapTo360(sat_coord(2, :) * cspice_dpr);
El = sat_coord(3, :) * cspice_dpr;
range = sat_coord(1, :);
Coord_CW = [range; Az; El];

%Implement FFRF noise

measures_CW = mvnrnd(Coord_CW',FFRF);
measures_2b = [t_grid', measures_CW]';  

% c)
t0_LVLH = vis_time(end)+5; 
tf_LVLH = t0_LVLH + 20*60;
t_grid1 = t0_LVLH:5:tf_LVLH;
% Save positions of the new time grid
ind_1 = 1:length(measures_2b);
ind_2 = 1:length(measures_2b);
ind_1 = ind_1(t0_LVLH ==  t_grid);
ind_1 = ind_1(1);
ind_2 = ind_2(tf_LVLH == t_grid);
ind_2 = ind_2(1);
% Save the updated states and measures for the new time grid wrt LVLH
TANGO.states = x_CW(:,ind_1:ind_2);
TANGO.sat_measurements = (measures_2b(:,ind_1:ind_2));

TANGO.measure = [[t_grid1(1),0,0,0]',TANGO.sat_measurements]';
% Initialize the data
P_UKF_CW = zeros(6,6,length(TANGO.measure(:,1)));
x_UKF_CW = zeros(6,length(TANGO.measure(:,1)));
sigma_points_UKF_CW = zeros(6,13,length(TANGO.measure(:,1)));
sigma_r_UKF_CW = zeros(length(TANGO.measure(:,1))-1,1);
sigma_v_UKF_CW = zeros(length(TANGO.measure(:,1))-1,1);

x_UKF_CW(:,1) = TANGO.states(:,1);
P_UKF_CW(:,:,1) = diag([0.01, 0.01, 0.1, 0.0001, 0.0001, 0.001]); 

% Propagate the data through Unsceented Kalman Filter
valueoftruth = 1;
for k = 2:length(TANGO.measure(:,1))
    [sigma_points_UKF_CW(:,:,k),~,~,x_UKF_CW(:,k),P_UKF_CW(:,:,k)] =  UKF(x_UKF_CW(:,k-1),P_UKF_CW(:,:,k-1),SVAL,TANGO.measure(k,2:4) ,0,[TANGO.measure(k-1,1) TANGO.measure(k,1)],valueoftruth,n_CW); 
end
% Calculate the errors
for i = 2 : length(TANGO.measure(:,1))
     sigma_r_UKF_CW(i-1) = 3*sqrt(trace(P_UKF_CW(1:3,1:3,i)));
     sigma_v_UKF_CW(i-1) = 3*sqrt(trace(P_UKF_CW(4:6,4:6,i)));
end
% Calculate the error
error_r = vecnorm(TANGO.states(1:3,:) - x_UKF_CW(1:3,2:end),2,1);
error_v = vecnorm(TANGO.states(4:6,:) - x_UKF_CW(4:6,2:end),2,1);

% Plot
time_plot = datetime(cspice_timout(TANGO.measure(2:end,1)',...
    'YYYY-MM-DD HR:MN:SC.###'),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');

figure()
subplot(1,2,1)
semilogy(time_plot,error_r,'LineWidth',2)
hold on
semilogy(time_plot,sigma_r_UKF_CW,'LineWidth',2)
grid on
xlabel('Time')
ylabel('Range [km]')
legend('error','3$\sigma$','Interpreter','latex')

subplot(1,2,2)
semilogy(time_plot,error_v,'LineWidth',2)
hold on
semilogy(time_plot,sigma_v_UKF_CW,'LineWidth',2)
grid on
ylabel('Velocity [km/s]')
xlabel('Time')
legend('error','3$\sigma$','Interpreter','latex')


%% Exercise 3
% a) 
% Initialize data
t_grid2 = [t_grid1(1) - 5, t_grid1];

P_TAN_rec = zeros(6,6,length(t_grid2));
x_TAN_rec = zeros(6,length(t_grid2));

x_TAN_rec(:,1) = x_UKF(:,end);
P_TAN_rec(:,:,1) = P_UKF(:,:,end);
% Propagate the states
valueoftruth = 1;
for k=2:length(t_grid2)
    [~,x_TAN_rec(:,k),P_TAN_rec(:,:,k),~,~] =  UKF(x_TAN_rec(:,k-1),P_TAN_rec(:,:,k-1),SVAL,TANGO.measure(k,2:4) ,J2,[t_grid2(k-1) t_grid2(k)],valueoftruth); 
end
% b) Provide the covariance in ECI ref frame
P_rot_ECI = zeros(6,6,length(t_grid2));
for k = 2 : length(t_grid2) 
    R = LVLHTrans(x_TAN_rec(1:3,k),x_TAN_rec(4:6,k),n_CW);
    P_rot_ECI(:,:,k) = R' * P_UKF_CW(:,:,k) * R;
end

% c)
P_TAN_ABS = zeros(6,6,length(t_grid2)-1);
sigma_r_TAN_ABS = zeros(1,length(t_grid2)-1);
sigma_v_TAN_ABS = zeros(1,length(t_grid2)-1);
% Calculate the sigmas
for k = 2:length(t_grid2)
    P_TAN_ABS(:,:,k-1) = P_rot_ECI(:,:,k) + P_TAN_rec(:,:,k-1);
    sigma_r_TAN_ABS(k-1) = 3*sqrt(trace(P_TAN_ABS(1:3,1:3,k-1)));
    sigma_v_TAN_ABS(k-1) = 3*sqrt(trace(P_TAN_ABS(4:6,4:6,k-1)));
end


% Plot
time_plot = datetime(cspice_timout(t_grid2(2:end),...
    'YYYY-MM-DD HR:MN:SC.###'),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');

figure()
subplot(1,2,1)
semilogy(time_plot,sigma_r_TAN_ABS,'LineWidth',2)
grid on
xlabel('Time')
ylabel('Range [km]')
legend('3$\sigma$','Interpreter','latex')

subplot(1,2,2)
semilogy(time_plot,sigma_v_TAN_ABS,'LineWidth',2)
grid on
ylabel('Velocity [km/s]')
xlabel('Time')
legend('3$\sigma$','Interpreter','latex')

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
function [xf, tt, xx,tf, PHIf] = keplerian_propagator(et0, x0, etf, attractor)
%KEPLERIAN_PROPAGATOR Propagate a Keplerian Orbit and the STM
%   Detailed explanation goes here

% Initialize propagation data
if isfloat(attractor)
    GM = attractor;
else
    GM = cspice_bodvrd(attractor, 'GM', 1);
end

Phi = eye(6);

tof = etf-et0;

options = odeset('reltol', 1e-13, 'abstol',1e-13); 

% Perform integration
[tt, xx] = ode78(@(t,x) keplerian_rhs(t,x,GM), [0 tof], x0, options);

% Extract state vector 
xf = xx(end,1:6)';
PHIf = Phi;
tf = tt(end);

end


function [xf, tt, xx,tf] = keplerian_propagator_J2(et0, x0, etf, attractor,J2,Re)
%KEPLERIAN_PROPAGATOR Propagate a Keplerian Orbit and the STM
%   Detailed explanation goes here

% Initialize propagation data
if isfloat(attractor)
    GM = attractor;
else
    GM = cspice_bodvrd(attractor, 'GM', 1);
end

options = odeset('reltol', 1e-13, 'abstol',1e-13);
x0Phi0 = x0;
% Perform integration
[tt, xx] = ode113(@(t,x) keplerian_rhs_J2(t,x,GM,J2,Re), [et0 etf], x0Phi0, options);

% Extract state vector 
xf = xx(end,1:6)';
tf = tt(end);

end

function [dxdt] = keplerian_rhs_J2(et, x, GM,J2,Re)
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
rotm = cspice_pxform('J2000','ITRF93',et);
r_ECEF = rotm * rr;
a_J2_ecef = 1.5*GM*J2*r_ECEF/norm(rr)^3*(Re(1)/norm(rr))^2 .* (5*(r_ECEF(3)/norm(rr))^2 - [1;1;3]);

a_J2 = rotm' * a_J2_ecef; 
% Position detivative is object's velocity
dxdt(1:3) = x(4:6);   
% Compute the gravitational acceleration using Newton's law
dxdt(4:6) = - GM * rr /(dist*dist2) + a_J2;

end

function [coord,r_ECI,v_ECI,et0,Station, varargout] = measurements(SC, Station, varargin)
% Function that simulates the measurements
% Input: SC[struct]: The structure containing all the data regarding the Satellite
%        Station[struct]:  The structure containing all the data regarding the Station
%        varargin: vector containing the reference times for the measurements
% Output: coord: The matrix containing the Angular measurements. It can contain also
%         vector indexes at which the Minimum Elevation conditiion is respected
%        r_ECI: Satellite position vector in ECI ref frame
%        v_ECI: Satellite velocity vector in ECI ref frame
%        et0[1]: Reference time obtained from the TLE
%        Station[struct]: The updated information structure of the Station 
%        varargout{1}: Vector containing the visibility Time
%        varargout{2}: Matrix containing the measurements that respect the minimum 
%        elevation constraint
% Simulate measurements
typerun = 'u';           % From the slides
typeinput = 'e';
whichconst = 72;
opsmode = 'a';

[Station.ECI, Station.ECI2TOPO] = stationcoordinates(Station, Station.vis_time);
%Extract the data from TLE
[satrec, ~, ~,~] = twoline2rv(SC.TLE1, SC.TLE2, typerun, typeinput, opsmode, whichconst);
%Extract the exact date from the data retrieved from TLE
[year1, mon1, day1, hr1, min1, sec1] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
if hr1 < 10 a = '0';else a = ''; end 
if min1 < 10 b = '0';else b = '';end
if sec1 < 10 c = '0';else c = '';end
if mon1 < 10 d = '0';else d = '';end
if day1 < 10 e = '0';else e = '';end
% Set the reference time value
time_str = strcat(string(year1),'-',d,string(mon1),'-',e,string(day1),"T",a,string(hr1),":",b,string(min1),":",c,string(sec1));
et0 = cspice_str2et(cellstr(time_str));

% Compute s/c states during visibility windows
if nargin == 2
    t0 = et0;
    t_span = 0;
else
    if length(varargin{1}) == 1
        t0 = varargin{1};
    else
        t0 = Station.vis_time;
    end
    t_span = (t0 - et0) / 60;  % Convert time span to minutes
end

n = length(t_span);

% For the2eci you need velocity distance and acceleration in teme, but we
% don't care abot velocity thus I take acceleration as zero
r_TEME = zeros(3,n);
v_TEME = zeros(3,n);
a_TEME = zeros(3,n);

% Nutation correction
dpsi = -0.073296 * pi/(180 * 3600); % arcseconds into [rad]
deps = -0.009373 * pi/(180 * 3600);

r_ECI = zeros(3, n);
v_ECI = zeros(3, n);

for i = 1:n
    a_TEME(:,i) = zeros(3,1);
    % Retrieve the data from sgp4 into TEME ref frame
    [~, r_TEME(:, i), v_TEME(:, i)] = sgp4(satrec, t_span(i));
    
    t = cspice_unitim(t0(i), 'ET', 'TDT') / cspice_jyear() / 100;
    %Convert the data from TEME into ECI, with the respective nutation
    [r_ECI(:, i), v_ECI(:, i), ~] = teme2eci(r_TEME(:, i), v_TEME(:, i), a_TEME(:,i), t, dpsi, deps);
end

[coord] = SCGScoord(Station, [r_ECI; v_ECI]);

if nargout >= 5
    rng('default')
    % Create the initial measurements including the noise erros
    measures = mvnrnd(coord(1:3, :)', Station.R);
    %Save the positions where error is larger than the limit elevation
    saveposition = measures(:, 3) >= Station.minel;
    %Update the measures
    measures = measures(measures(:, 3) >= Station.minel, :);
    %Calculate the visibility time wrt minimum Elevation
    vis_time = Station.vis_time(saveposition);
    %Output the data
    varargout{1} = vis_time';
    varargout{2} = measures;
    coord = [coord; saveposition'];
end

end

function [Coord] = SCGScoord(Station, x_MAN_span)
% Function Transforming the measurements in angular measurements 
% Input: x_MAN_span: The initial states that must be converted in angluar coordinates
%        Station[struct]: The structure containing all the information
%        regarding the Station
% Output: Coord: The initial coordinated in angular coordinates
    SC_ECI = x_MAN_span - Station.ECI;
    SC_TOPO = zeros(6,length(SC_ECI(1,:)));
    % Obtain the station's position in Topocentric reference frame
    for i = 1:length(SC_ECI(1,:))
        SC_TOPO(:,i) = Station.ECI2TOPO(:,:,i) * SC_ECI(:,i);
    end
    % Retrive the angular measurements data from the cspice kernels    
    SC_ang_coord = cspice_xfmsta(SC_TOPO,'RECTANGULAR','LATITUDINAL','EARTH');
    SC_ang_coord(2,:) = wrapTo360(SC_ang_coord(2,:)*180/pi);
    Coord = [SC_ang_coord(1,:);SC_ang_coord(2,:);SC_ang_coord(3,:)*180/pi];
end

function [Coord] = SCGScoord1(Station, x_MAN_span)
% Function Transforming the measurements in angular measurements 
% Input: x_MAN_span: The initial states that must be converted in angluar coordinates
%        Station[struct]: The structure containing all the information
%        regarding the Station
% Output: Coord: The initial coordinated in angular coordinates
    SC_ECI = x_MAN_span - Station.ECI_UKF;
    SC_TOPO = zeros(6,length(SC_ECI(1,:)));
    % Obtain the station's position in Topocentric reference frame
    for i = 1:length(SC_ECI(1,:))
        SC_TOPO(:,i) = Station.ECI2TOPO_UKF(:,:,i) * SC_ECI(:,i);
    end
    % Retrive the angular measurements data from the cspice kernels
    SC_ang_coord = cspice_xfmsta(SC_TOPO,'RECTANGULAR','LATITUDINAL','EARTH');
    SC_ang_coord(2,:) = wrapTo360(SC_ang_coord(2,:)*180/pi);
    Coord = [SC_ang_coord(1,:);SC_ang_coord(2,:);SC_ang_coord(3,:)*180/pi];
end

function [Coord,varargout] = SCGScoord2(states,varargin)
% Function Transforming the measurements in angular measurements used for
% FFRF case
% Input: states: The initial states that must be converted in angluar coordinates
%        varargin: The noise matrix required into simulating the measurements
% Output: Coord: The initial coordinated in angular coordinates
%         varargout: The measurements obtained with the noise matrix

        sat_coord = cspice_xfmsta(states, 'RECTANGULAR', 'LATITUDINAL', 'Earth');
        Az = wrapTo360(sat_coord(2, :) * cspice_dpr);
        El = sat_coord(3, :) * cspice_dpr;
        range = sat_coord(1, :);
        Coord = [range; Az; El];

        if nargout == 2
            %Implement FFRF noise
            R = varargin{1};
            % Simulate the Gaussian noise
            rng('default')
            varargout{1} = mvnrnd(Coord',R);
        end
end

function [ECI, ECI2TOPO] = stationcoordinates(Station,tspan)
% Input: Station[Struct]: the structure containing all the Station's data
%        tspan: the time vector for which the ECI should be calculated and ECI2TOPO
% Output: ECI: Position matrix for the desired time grid at of the Station
%         ECI2TOPO: Transformation matrix from ECI into Topocentric reference frame
%         for the entire period of the time grid with respect to J2000
    % Extract the ECI positionn of the station and the Transformation
    % matrix from ECI to TOPO reference frame
    ECI2TOPO = cspice_sxform('J2000',Station.TOPO,tspan);
    ECI  = cspice_spkezr(Station.name,tspan,'J2000','NONE','EARTH');

end

function [sigma_points,x_k_min,P_k_min,varargout]  = UKF(x,Px,station,measure,J2,t,varargin)
% The UKF function has the purpose of propagating the data throught the
% Unscented Kalman Filter
% Input:  x[6x1]: The initial mean state of the satellite that wants to be
%         propagated
%         Px[6x6]: The initial covariance of the satellite that wants to be
%         propagated
%         station[struct]: The structure containing all the information
%         regarding the Station
%         measure []:
%         J2[1]: The J2 perturbation value
%         t: The propagation time vector, containing the initial epoch of
%         the propagation(et0) and the final propagation epoch(etf)
%         varargin{1}: Contains the value of the stage of the propagation:
%                         - 0 it stops the propagation after computing the apriori covariance and mean state
%                         - 1 it continues the propagation to compute the aposteriori covariance and mean state
%           varargin{2}: The value of the Satellite's mean motion
% Output:   sigma_points[6, 2n]: The matric containing the sigma points
%             x_k_min[6xm]: The apriori mean states, saved in a matrix
%             P_k_min[6,6,m]: The apriori covariance, saved in set of m 6x6 matrices
%             varargout{1}[6xm]: The aposterori mean states, saved in a matrix
%             varargout{2}[6,6,m]: The aposteriori covariance, saved in a matrix
% Initialize the data
t0 = t(1);
tf = t(2);
n = 6;
Re = cspice_bodvrd('EARTH','RADII',3);
Re = Re(1);
mu = cspice_bodvrd('EARTH','GM',1);
alpha = 0.1;
beta  = 2;
lambda = alpha^2 * n - n;

% Initialize the weight matrices
W0_m   = lambda/(n+lambda);
W0_c   = lambda/(n + lambda) + (1 - alpha^2 + beta);
Wi   = 1/(2 * (n + lambda));
Wm = [W0_m;Wi*ones(2*n,1)]';
Wc = [W0_c;Wi*ones(2*n,1)]';

% Initialzie the sigma points
sigma_matrix = sqrtm((n+lambda)*Px);
sigma_points = [x, x + sigma_matrix(:,1:6), x - sigma_matrix(:,1:6)];
if t0 < tf
    for i = 1:(2 * n + 1)
         if J2 ~=0   
             % Compute the propagated sigma points through Keplerian J2
             % propagato
                [x_ode, ~, ~,~] = keplerian_propagator_J2(t0, sigma_points(:,i), tf, mu,J2,Re);
            
         else
             % Compute propagates sigma points using the Clohessy-Wiltshire equations
                n_CW = varargin{2};
                [~,x_ode] = CW(sigma_points(:,i),t0,tf,n_CW);
        end
        sigma_points(:, i) = x_ode;
    end
    % Compute apriori mean state
    x_k_min = zeros(6, 1);
    for i = 1:size(sigma_points, 2)
        x_k_min = x_k_min + Wm(i) * sigma_points(:, i);
    end
    % Compute a priori covariance
    P_k_min = zeros(6, 6);
    for i = 1:size(sigma_points, 2)
        diff = sigma_points(:,i) - x_k_min;
        P_k_min = P_k_min + Wc(i) * diff * diff';
    end
else
    x_k_min = x;
    P_k_min = Px;
end

if ~varargin{1}
    return;
end

% Initialize the initial value of the meaurements covariance
Pyy = station.R;

if J2 == 0
 Pyy = diag([0.00001^2 1 1]);
end
n = 6;

% Simulate measurements using propagated sigma points
why_k = zeros(3, 2 * n + 1);
for i = 1:size(why_k, 2)

    if J2 ~= 0
            why_k(:, i) = SCGScoord1(station, sigma_points(:, i));
    end
    if J2 == 0
            why_k(:, i) = SCGScoord2(sigma_points(:, i));
    end
end

% Compute the predicted measurement
y = zeros(3,1);
for k = 1: 2 * n + 1
    y = y + Wm(k)*why_k(:,k);
end
% Compute the measurements covariance 
diff_y = zeros(3,size(why_k,2));
for k = 1:size(why_k,2)
     diff_y(:,k) =  [why_k(1, k)' - y(1); ...
         180/pi * angdiff(y(2:3) * pi/180, why_k(2:3,k) * pi/180)];

    Pyy = Pyy + Wc(k) * diff_y(:,k) * diff_y(:,k)';
end

% Compute the measurements cross_covaariance with the states
Pxy    = zeros(6,3);
diff_x = zeros(6,13); 
for k = 1:size(why_k,2)
    diff_x(:,k) = sigma_points(:,k) - x_k_min;
    Pxy = Pxy + Wc(k) * diff_x(:,k) * diff_y(:,k)';
end

% Perform correction if Pyy is invertible
if cond(Pyy) < 1e10

    %Kalman gain
    K = Pxy / Pyy;

    %Correction of the state
     diff = [measure(1)' - y(1); ...
         180/pi * angdiff(y(2:3) * pi/180, measure(2:3)' * pi/180)];

    x_k_plus = x_k_min + K * diff;

    %Propagation of covariance matrix
    P_k_plus = P_k_min - K * Pyy * K';
end

varargout{1} = x_k_plus;
varargout{2} = P_k_plus;
end

function R = LVLHTrans(r_eci,v_eci,n)
% This function provides the transformationn matrix from ECI frame to LVLH
% Input:  r_eci[3x1]: the position vector of the satellite
%         v_eci[3x1]: the velocity vector of the satellite
%         n[1]:       the satellite's mean motion
h_vec = cross(r_eci,v_eci);
h = h_vec/norm(h_vec);
i = r_eci / norm(r_eci);
k = h;
j = cross(k,i);

% Assemble the transformation matrix from ECI to LVLH
R = [ [i'; j'; k']             zeros(3);
      [0 n 0; -n 0 0; 0 0 0]*[i';j';k']  [i'; j'; k'] ];

end

function [tf,xf] = CW(x0,et0,etf,n)
% Propagator that is based on the Clohessy-Wiltshire (CW)
% Input: x0 [6x1] the initial state of the Satellite that is propagated
%        et0[1] the initial epoch of the propagation
%        etf[1] the final epoch of the propagation
%        n[1] Satellite's mean motion
options = odeset('reltol', 1e-13, 'abstol',1e-13);
% Perform integration
[tt, xx] = ode113(@(t,x) CW_rhs(t,x,n), [et0 etf], x0, options);

% Extract state vector 
xf = xx(end,1:6)';
tf = tt(end);
end

function dxdt = CW_rhs(t,x,n)
%This function is used in Ode to propagate the satellite under the Clohessy-
%Wiltshire (CW) equations
% Input: t[1]: not used
%        x[6x1]: the states of the satellite
%        n[1]:mean motion of satellite
% Output: dxdt[6x1]: the data input for the ODE, composed of the velocities
% and the positions

dxdt = zeros(6,1);
dxdt(1:3) = x(4:6);
dxdt(4)   = 3 * n^2 *x(1) + 2*n*x(5);
dxdt(5)   = -2 * n * x(4);
dxdt(6)   = - n^2 * x(3);

end
