% Spacecraft Guidance and Navigation 2023/2024
% Assignment 2 Exercise 2
% Martin Mihnea Stefan

clc;clearvars;close all;

addpath('.\kernels\')
addpath('.\sgp4\')
addpath('.\tle\')
addpath('.\mice\lib')
addpath('.\mice\src\mice')

cspice_furnsh('assignment02.tm');

% Initial states of Satellite 1 (MAN)
r_MAN = [4622.232026629, 5399.3369588058, -0.0212138165769957];
v_MAN = [0.812221125483763, -0.721512914578826, 7.42665302729053];
x_MAN = [r_MAN,v_MAN]';

% Initial states of Satellite 2 (TAN)
r_TAN = [4621.69343340281, 5399.26386352847, -3.09039248714313];
v_TAN = [0.813960847513811, -0.719449862738607, 7.42706066911294];
x_TAN = [r_TAN,v_TAN]';


% Initial covariance matrix
P0 = [5.6e-7 3.5e-7 -7.1e-8 0 0 0;
      3.5e-7 9.7e-7 7.6e-8 0 0 0; 
      -7.1e-8 7.6e-8 8.1e-8 0 0 0;
      0 0 0 2.8e-11 0 0;
      0 0 0 0 2.7e-11 0;
      0 0 0 0 0 9.6e-12];

% Convert time strings to ephemeris time (ET)
tsep = cspice_str2et('August 12 05:27:39.114 UTC 2010');
t0   = cspice_str2et('2010-08-12T05:30:00.000');
tf   = cspice_str2et('2010-08-12T11:00:00.000');

% Define Time Span
t_grid = t0:60:tf;

% Ground station definitions
SVAL = struct('name','SVALBARD','LAT',78.229772,'LONG',15.407786,'ALT',  458 , 'minel',5, 'TYPE','RADAR', 'R', diag([0.01^2 0.125^2 0.125^2]), 'TOPO','SVALBARD_TOPO', 'W', diag([1/0.01^2 1/0.125^2 1/0.125^2]));
KOUROU = struct('name','KOUROU','LAT',5.25144,  'LONG',-52.80466,'ALT' ,-14.67,'minel',10,'TYPE','RADAR', 'R', diag([0.01^2 0.1^2 0.1^2]),     'TOPO','KOUROU_TOPO',   'W', diag([1/0.01^2 1/0.1^2 1/0.1^2]));

% Satellite information
MANGO = struct('name', 'MANGO','ID', 36599, 'Mango',x_MAN, 'COV',P0,'TLE1','1 36599U 10028B   10224.22752732 -.00000576  00000-0 -16475-3 0  9998', ...
                                                                    'TLE2','2 36599 098.2803 049.5758 0043871 021.7908 338.5082 14.40871350  8293');
TANGO = struct('name', 'TANGO','ID', 36827, 'Tango',x_TAN, 'COV',P0,'TLE1','1 36827U 10028F   10224.22753605  .00278492  00000-0  82287-1 0  9996', ...
                                                                    'TLE2','2 36827 098.2797 049.5751 0044602 022.4408 337.8871 14.40890217    55');

% Earth's gravitational parameter
mu_E = cspice_bodvrd('EARTH','GM',1);

% Propagate mean states of satellites from Separation to t0
[xf,~,~,~] = keplerian_propagator(tsep, x_MAN, t0, mu_E);
x_MAN_span(:,1) = xf; 
[xf,~,~,~] = keplerian_propagator(tsep, x_TAN, t0, mu_E);
x_TAN_span(:,1) = xf; 

% Propgate mean states of satellites for entire period
for i = 1:length(t_grid)-1
    [x_MAN_span(:,i+1),~,~,~] = keplerian_propagator(t_grid(i),x_MAN_span(:,i),t_grid(i+1),mu_E);
    [x_TAN_span(:,i+1),~,~,~] = keplerian_propagator(t_grid(i),x_TAN_span(:,i),t_grid(i+1),mu_E);
end

% Earth's J2 coefficient
J2 = 0.0010826269;

% Earth's mean radius
Re = cspice_bodvrd('EARTH','RADII',3);
Re = Re(1);


%% Part 1 

% Calculate ECI and ECI2TOPO coordinates for KOUROU and SVALBARD stations
% at each time step
[KOUROU.ECI,KOUROU.ECI2TOPO] =  COORDStation(KOUROU,t_grid);
[SVAL.ECI,SVAL.ECI2TOPO] =  COORDStation(SVAL,t_grid);

% Compute satellite coordinates with respect to KOUROU and SVALBARD stations
KOUROU.Mango_coord = SCGScoord(KOUROU, x_MAN_span);
SVAL.Mango_coord = SCGScoord(SVAL, x_MAN_span);

% Determine visibility time and position for KOUROU station
KOUROU.vis_time = t_grid(KOUROU.Mango_coord(3,:) >= KOUROU.minel );
KOUROU.vis_position = KOUROU.Mango_coord(:,KOUROU.Mango_coord(3,:) >= KOUROU.minel);

% Determine visibility time and position for SVALBARD station
SVAL.vis_time = t_grid(SVAL.Mango_coord(3,:) >= SVAL.minel );
SVAL.vis_position = SVAL.Mango_coord(:,SVAL.Mango_coord(3,:) >= SVAL.minel);

% Generate a vector of evenly spaced times for plotting
time_VEC_NR_SVAL = linspace(SVAL.vis_time(1),SVAL.vis_time(end),10);
time_VEC_NR_KOUROU = linspace(KOUROU.vis_time(1),KOUROU.vis_time(end),4);
time_vec_SVAL = cspice_et2utc(linspace(SVAL.vis_time(1),SVAL.vis_time(end),10),'C',3);
time_vec_KOUROU = cspice_et2utc(linspace(KOUROU.vis_time(1),KOUROU.vis_time(end),4),'C',3);
time_plot_tot_SVAL = string(time_vec_SVAL(:,6:20));
time_plot_tot_KOUROU = string(time_vec_KOUROU(:,6:20));

% Plot SVALBARD

figure(1)
hold on;
plot(SVAL.vis_time, SVAL.vis_position(2,:),'o','MarkerFaceColor','#EDB120','DisplayName','Azimuth[deg]')
plot(SVAL.vis_time, SVAL.vis_position(3,:),'s','MarkerFaceColor','#D95319','DisplayName','Elevation[deg]')
xlabel('Time','FontSize',14,'FontName','Arial')
ylabel('Svalbard','FontSize',14,'FontName','Arial')
xticks(time_VEC_NR_SVAL)
xticklabels(time_plot_tot_SVAL)
xtickangle(0);
ax = gca;
ax.XAxis.FontSize = 7.5;
ax.XLabel.FontSize = 14;
legend
grid on;
hold off;
% Plot KOUROU
figure(2)
hold on;
grid on;
plot(KOUROU.vis_time, KOUROU.vis_position(2,:),'o','MarkerFaceColor','b','DisplayName','Azimuth[deg]')
plot(KOUROU.vis_time, KOUROU.vis_position(3,:),'s','MarkerFaceColor','#A2142F','DisplayName','Elevation[deg]')
xlabel('Time','FontSize',14,'FontName','Arial')
ylabel('Kourou','FontSize',14,'FontName','Arial')
xticks(time_VEC_NR_KOUROU)
xticklabels(time_plot_tot_KOUROU)
ax = gca;
ax.XAxis.FontSize = 7.5;
ax.XLabel.FontSize = 14;
xtickangle(0);
legend
hold off;

%% Part 2 
% Obtain the Update data for SVALBARD, as the measurements(Angular
% measurements), and the States measurements obtained from SVALBARD
[coord_SVAL,~,~,~,SVAL,SVAL.MANGOmeasure] = measurements(MANGO, SVAL,SVAL.vis_time);
% Obtain the Update data for KOUROU, as the measurements(Angular
% measurements), and the States measurements obtained from KOUROU
[coord_KOUROU,~,~,~,KOUROU,KOUROU.MANGOmeasure] = measurements(MANGO,KOUROU,KOUROU.vis_time);

%% Part 3
% a) Define the lsqnonlin to obtain the measurements
[~,r_ECI,v_ECI,~,~] = measurements(MANGO, KOUROU,t0);
opt = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','iter');%,'MaxIterations',75);

[x_K,resnorm_K,residual_K,exitflag_K,~,~,jac_K] = lsqnonlin(@(x) costfunction(x, {KOUROU},t0,mu_E),x_MAN_span(:,1),[],[],opt);
Jac = full(jac_K);                                               
P_ls = resnorm_K/(length(residual_K) - length(x_MAN_span)) .*inv(Jac'*Jac);
%%
% b) 
% Repeat the process for both stations
BOTH_STA = {KOUROU;SVAL};
[x_SK,resnorm_SK,residual_SK,exitflag_SK,~,~,jac_SK] = lsqnonlin(@(x) costfunction(x, BOTH_STA,t0,mu_E),x_MAN_span(:,1),[],[],opt);
Jac_SK = full(jac_SK);
P_ls_SK = resnorm_SK/(length(residual_SK) - length(x_MAN_span)) .*inv(Jac_SK'*Jac_SK);

%% 
% c) 
% Repeat the process for both stations, including the J2 perturbation
[xf_MAN_span_J2(:,1),~,~,~,~] = keplerian_propagator_J2(tsep, x_MAN, t0, mu_E,J2,Re);

BOTH_STA = {KOUROU;SVAL};
[x_SK_J2,resnorm_SK_J2,residual_SK_J2,exitflag_SK_J2,~,~,jac_SK_J2] = lsqnonlin(@(x) costfunction(x, BOTH_STA,t0,mu_E,'True'),xf_MAN_span_J2(:,1),[],[],opt);
Jac_SK_J2 = full(jac_SK_J2);
P_ls_SK_J2 = resnorm_SK_J2/(length(residual_SK_J2) - length(xf_MAN_span_J2(:,1))) .*inv(Jac_SK_J2'*Jac_SK_J2);

%% Part 5
% Solution for Kourou with J2 perturbation
[xf_MAN_span_K_J2(:,1),~,~,~,~] = keplerian_propagator_J2(tsep, x_MAN, t0, mu_E,J2,Re);
[x_K_J2,resnorm_K_J2,residual_K_J2,exitflag_K_J2,~,~,jac_K_J2] = lsqnonlin(@(x) costfunction(x, {KOUROU},t0,mu_E,'True'),xf_MAN_span_K_J2(:,1),[],[],opt);
Jac_K_J2 = full(jac_K_J2);
P_ls_K_J2 = resnorm_K_J2/(length(residual_K_J2) - length(xf_MAN_span_K_J2)) .*inv(Jac_K_J2'*Jac_K_J2);

% Solution for Svalbard
[~,r_ECI_S,v_ECI_S,~,~] = measurements(MANGO, SVAL,t0);
% Calculate the solution for Svalbard without J2
[x_S,resnorm_S,residual_S,exitflag_S,~,~,jac_S] = lsqnonlin(@(x) costfunction(x, {SVAL},t0,mu_E),x_MAN_span(:,1),[],[],opt);
Jac_S = full(jac_S);                                               
P_ls_S = resnorm_K/(length(residual_S) - length(x_MAN_span)) .*inv(Jac_S'*Jac_S);

% SOLUTION FOR SVALBARD IN J2
[KOUROU.ECI,KOUROU.ECI2TOPO] =  COORDStation(KOUROU,t_grid);
[SVAL.ECI,SVAL.ECI2TOPO] =  COORDStation(SVAL,t_grid);
% Use the J2 propagator
[xf_MAN_span_S_J2(:,1),~,~,~,~] = keplerian_propagator_J2(tsep, x_MAN, t0, mu_E,J2,Re);
% Lsqnonlin for J2 for Svalbard station
[x_S_J2,resnorm_S_J2,residual_S_J2,exitflag_S_J2,~,~,jac_S_J2] = lsqnonlin(@(x) costfunction(x, {SVAL},t0,mu_E,'True'),xf_MAN_span_S_J2(:,1),[],[],opt);
Jac_S_J2 = full(jac_S_J2);
P_ls_S_J2 = resnorm_S_J2/(length(residual_S_J2) - length(xf_MAN_span_S_J2(:,1))) .*inv(Jac_S_J2'*Jac_S_J2);

%Calculate the errors for each method
error_K = [norm(r_ECI - x_K(1:3)),norm(v_ECI-x_K(4:6))];
error_SK = [norm(r_ECI - x_SK(1:3)),norm(v_ECI - x_SK(4:6))];
error_SK_J2 = [norm(r_ECI - x_SK_J2(1:3)),norm(v_ECI - x_SK_J2(4:6))];
error_K_J2 = [norm(r_ECI - x_K_J2(1:3)),norm(v_ECI - x_K_J2(4:6))];
error_S = [norm(r_ECI - x_S(1:3)),norm(v_ECI-x_S(4:6))];
error_S_J2 = [norm(r_ECI - x_S_J2(1:3)),norm(v_ECI - x_S_J2(4:6))];

% The best solution is considering both stations and J2 model
%%
[KOUROU.Tango_coord,~,~,et0P,~,KOUROU.MANGOmeasure ] = measurements(TANGO, KOUROU, t_grid);
[SVAL.Tango_coord,~,~,~,~,SVAL.MANGOmeasure ] = measurements(TANGO, SVAL, t_grid);

[~,r_ECI_SK,v_ECI_SK,~,~] = measurements(TANGO,SVAL,t0);              

[xf_TAN_span_J2(:,1),~,~,~,~] = keplerian_propagator_J2(tsep, x_TAN, t0, mu_E,J2,Re);% ok

BOTH_STA = {KOUROU;SVAL};
[x_SK_J2_TAN,resnorm_SK_J2_TAN,residual_SK_J2_TAN,exitflag_SK_J2_TAN,~,~,jac_SK_J2_TAN] = lsqnonlin(@(x) costfunction(x, BOTH_STA,t0,mu_E,'True'),xf_TAN_span_J2(:,1),[],[],opt);
Jac_SK_J2_TAN = full(jac_SK_J2_TAN);
P_ls_SK_J2_TAN = resnorm_SK_J2_TAN/(length(residual_SK_J2_TAN) - length(xf_TAN_span_J2(:,1))) .*inv(Jac_SK_J2_TAN'*Jac_SK_J2_TAN);

error_TAN = [norm(r_ECI_SK - x_SK_J2_TAN(1:3)),norm(v_ECI_SK - x_SK_J2_TAN(4:6))];
sigma_r_TAN = 3*sqrt(trace(P_ls_SK_J2_TAN(1:3,1:3)));
sigma_v_TAN = 3*sqrt(trace(P_ls_SK_J2_TAN(4:6,4:6))); 
%to here
%% Functions

function [xf, tt, xx,tf, PHIf] = keplerian_propagator(et0, x0, etf, attractor)
%KEPLERIAN_PROPAGATOR Propagate a Keplerian Orbit and the STM

% Initialize propagation data
if isfloat(attractor)
    GM = attractor;
else
    GM = cspice_bodvrd(attractor, 'GM', 1);
end

Phi = eye(6);


options = odeset('reltol', 1e-12, 'abstol',1e-13);
x0Phi0 = x0;
% Perform integration
[tt, xx] = ode78(@(t,x) keplerian_rhs(t,x,GM), [et0 etf], x0Phi0, options);

% Extract state vector 
xf = xx(end,1:6)';
PHIf = Phi;
tf = tt(end);

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

function [ECI, ECI2TOPO] =  COORDStation(Station,tspan)
% This function takes the data from the Statiion and provides its position
% in eci and the matrix to Transfrom from ECI reference frame into
% Topocentric reference frame
% Input: Station [struct]: contains all the data of the Station
%        tspan: [1xn]: the time vector for which the transformation must be
%        calculated
% Output: ECI [6X1]: the ECI states of the Station
%         ECI2TOPO [6x6]: the tranformation matrix of the Stationn from ECI
%        to Topocentric
    ECI2TOPO = cspice_sxform('J2000',Station.TOPO,tspan);
    ECI  = cspice_spkezr(Station.name,tspan,'J2000','NONE','EARTH');

end

function [Coord] = SCGScoord(Station, x_MAN_span)
% Function that converts the states into Angular coordinates
% INPUT: Station: The data that associates the Station
%        x_MAN_span: The state of the satellite that must be converted into
%        angular
    % Compute satellite coordinates in ECI with respect to the station
    SC_ECI = x_MAN_span - Station.ECI;
    % Initialize a Topocentric vector
    SC_TOPO = zeros(6, length(SC_ECI(1, :)));

    % Convert ECI coordinates to TOPO coordinates for each time step
    for i = 1:length(SC_ECI(1, :))
        SC_TOPO(:, i) = Station.ECI2TOPO(:, :, i) * SC_ECI(:, i);
    end
    
    % Convert TOPO coordinates to angular coordinates (azimuth, range and elevation)
    SC_ang_coord = cspice_xfmsta(SC_TOPO, 'RECTANGULAR', 'LATITUDINAL', 'EARTH');
    
    % Adjust Azimuth values to be within [0, 360] degrees
    for i = 1:length(SC_ang_coord(3, :))
        if SC_ang_coord(2, i) < 0
            SC_ang_coord(2, i) = 2*pi + SC_ang_coord(2, i);
        end
    end
    
    % Output the satellite coordinates 
    Coord = [SC_ang_coord(1, :); SC_ang_coord(2, :) * 180 / pi; SC_ang_coord(3, :) * 180 / pi];
end


function [coord,r_ECI,v_ECI,et0,Station, varargout] = measurements(SC, Station, varargin)
%Provides the coordinates of the satellite wrt the station, in Angular
%measurements, the states of the Satellite in ECI frame and position for
%time window
% Input: SC[struct]: the structure containing the data for the Satellite
%        Station[struct]: the structure containing the data for the Station
%        varargin{1} : the time vector for which the measurement it is
%        wanted
% Output: Coord: The coordintes and the saved position for the visibility
% time window
%        r_ECI [3,length(varargin{1})]: vector containing the position of
%        the Satellite in ECI referance frame
%        v_ECI [3,length(varargin{1})]: vector containing the veloocity of
%        the Satellite in ECI referance frame
%        et0: The time contained in TLE
%        Station[struct]: Updated Station structure
%        varargout{1}: The matrix containing the measurements including the
%        noise and the visibility time window values
% Simulate measurements
typerun = 'u';           %PDF
typeinput = 'e';
whichconst = 72;
opsmode = 'a';

[Station.ECI, Station.ECI2TOPO] =  COORDStation(Station, Station.vis_time);

[satrec, ~, ~,~] = twoline2rv(SC.TLE1, SC.TLE2, typerun, typeinput, opsmode, whichconst);
[year1, mon1, day1, hr1, min1, sec1] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
if hr1 < 10 a = '0';else a = ''; end 
if min1 < 10 b = '0';else b = '';end
if sec1 < 10 c = '0';else c = '';end
if mon1 < 10 d = '0';else d = '';end
if day1 < 10 e = '0';else e = '';end

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
    t_span = (t0 - et0) / 60;  
end

n = length(t_span);

% For the2eci you need velocity distance and acceleration in teme, but we
% don't care about velocity thus I take acceleration as zero
r_TEME = zeros(3,n);
v_TEME = zeros(3,n);
a_TEME = zeros(3,n);

% Nutation correction
dpsi = -0.0733 * pi/(180 * 3600); % arcseconds into [rad]
deps = -0.00937 * pi/(180*3600);


r_ECI = zeros(3, n);
v_ECI = zeros(3, n);

for i = 1:n
    a_TEME(:,i) = zeros(3,1);
    % Extract Teme position of the Satellite from sgp4 
    [~, r_TEME(:, i), v_TEME(:, i)] = sgp4(satrec, t_span(i));
    % Convert TEME into ECI
    [r_ECI(:, i), v_ECI(:, i), ~] = teme2eci(r_TEME(:, i), v_TEME(:, i), a_TEME(:,i), t0(i)/3155760000, dpsi, deps);
end
% Extract the Angular coordianted from ECI
[coord] = SCGScoord(Station, [r_ECI; v_ECI]);

if nargout == 6
    % Simulate the measurements including the noise
    measures = mvnrnd(coord(1:3, :)', Station.R);
    % Save the positions for the time window 
    saveposition = measures(:, 3) >= Station.minel;
    % Ectract the measurements that are inside the time window
    measures = measures(measures(:, 3) >= Station.minel, :);
    vis_time = Station.vis_time(saveposition);
    
    varargout{1} = [vis_time', measures];
    
    coord = [coord; saveposition'];
end
end

function residual = costfunction(x, Station, t0, mu_E, varargin)
% COSTFUNCTION Compute the residuals for the least squares optimization (lsqnonlin)
%   This function calculates the residuals between the predicted and actual
%   measurements of a satellite's position, used for optimization with lsqnonlin.
%
%   Inputs:
%       x       - Initial state vector of the satellite
%       Station - Cell array containing measurement data from ground stations
%       t0      - Initial epoch time
%       mu_E    - Gravitational parameter of Earth
%       varargin - Optional parameter to include J2 perturbation in the propagation
%
%   Outputs:
%       residual - Residuals between the predicted and actual measurements

    % Constants for J2 perturbation and Earth radius
    J2 = 0.0010826269;  % J2 perturbation constant
    Re = cspice_bodvrd('EARTH', 'RADII', 3);  % Earth radius from SPICE
    Re = Re(1);  % Extract equatorial radius

    % Initialize variables
    n = 1e+4;  % Maximum number of measurements
    measure = zeros(3, n)';  % Matrix to store measurements
    t_span = zeros(n, 1);  % Vector to store measurement times
    W = zeros(3, 3, n);  % Weight matrix
    i = 0;

    % Compile measurements and weight matrices from all ground stations
    for k = 1:size(Station, 1)
        j = i + size(Station{k, 1}.MANGOmeasure, 1);
        t_span(i+1:j) = Station{k, 1}.MANGOmeasure(:, 1);  % Measurement times
        measure(i+1:j, :) = Station{k, 1}.MANGOmeasure(:, 2:4);  % Measurement values
        W(:, :, i+1:j) = Station{k, 1}.W .* ones(3, 3, j-i);  % Weight matrix
        [St_ECI(:, i+1:j), St_ECI2TOPO(:, :, i+1:j)] = COORDStation(Station{k, 1}, Station{k, 1}.MANGOmeasure(:, 1)');
        i = j;
    end
    
    % Sort measurements by time
    t_span = t_span(1:i);
    [t_span, ind] = sort(t_span);
    measure = measure(ind, :);
    residual = zeros(size(measure));
    W = W(:, :, ind);
    StationFUNC.ECI = St_ECI(:, ind);
    StationFUNC.ECI2TOPO = St_ECI2TOPO(:, :, ind);

    % Propagate the initial state to the measurement times using Keplerian motion
    t = [t0; t_span];
    xf = zeros(6, length(t));
    xf(:, 1) = x;

    for k = 2:length(t_span) + 1
        [xf(:, k), ~, ~, ~, ~] = keplerian_propagator(t(k-1), xf(:, k-1), t(k), mu_E);
    end

    % Compute the angular measurements from the propagated states
    ang_final = SCGScoord(StationFUNC, xf(:, 2:end));
    for k = 1:length(t_span)
        residual(k, :) = W(:, :, k) * [ang_final(1, k) - measure(k, 1), 180/pi * angdiff(measure(k, 2:3) * pi/180, ang_final(2:3, k)' * pi/180)]';
    end

    % If additional arguments are provided, use J2-perturbed propagation
    if nargin >= 5
        xf = zeros(6, length(t));
        xf(:, 1) = x;
        t = [t0; t_span];
        for k = 2:length(t_span) + 1
            [xf(:, k), ~, ~, ~, ~] = keplerian_propagator_J2(t(k-1), xf(:, k-1), t(k), mu_E, J2, Re);
        end
        ang_final = SCGScoord(StationFUNC, xf(:, 2:end));
        for k = 1:length(t_span)
            residual(k, :) = W(:, :, k) * [ang_final(1, k) - measure(k, 1), 180/pi * angdiff(measure(k, 2:3) * pi/180, ang_final(2:3, k)' * pi/180)]';
        end
    end
end


function [xf, tt, xx,tf, PHIf] = keplerian_propagator_J2(et0, x0, etf, attractor,J2,Re)
%KEPLERIAN_PROPAGATOR Propagate a Keplerian Orbit and the STM
%   Detailed explanation goes here

% Initialize propagation data
if isfloat(attractor)
    GM = attractor;
else
    GM = cspice_bodvrd(attractor, 'GM', 1);
end

Phi = eye(6);

options = odeset('reltol', 1e-12, 'abstol',1e-13); %[ones(3,1).*1e-8; ones(3,1).*1e-11; 1e-12*ones(36,1)]);
x0Phi0 = x0;%Phi(:)];
% Perform integration
[tt, xx] = ode78(@(t,x) keplerian_rhs_J2(t,x,GM,J2,Re), [et0 etf], x0Phi0, options);

% Extract state vector 
xf = xx(end,1:6)';
PHIf = Phi;%PHIf = reshape(xx(end,7:end),6,6); 
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
a_J2_ecef = 1.5*GM*J2*r_ECEF/norm(rr)^3*(Re(1)/norm(rr))^2 .* (5*(r_ECEF(3)/norm(rr)).^2 - [1;1;3]);

a_J2 = rotm' * a_J2_ecef; 
% Position detivative is object's velocity
dxdt(1:3) = x(4:6);   
% Compute the gravitational acceleration using Newton's law
dxdt(4:6) = - GM * rr /(dist*dist2) + a_J2;

end
