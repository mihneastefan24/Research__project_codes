% Modelling and Simulation of Aerospace Systems(2024/2025)
% Assignment #1
% Author: Martin Mihnea Stefan - 10903720

%% Exercise 1

clc;
clearvars;
close all;

% Exercise 1.1.
% Initialize constants
a1 = 10;
a2 = 13;
a3 = 8;
a4 = 10;
ph = 1;
% Set tolerance for Newton Solver(NS)
tol = 1e-5;
% Set the tolerance for fzero to be the same asfor NS
options = optimset('TolFun', 10^-5, 'TolX',10^-5);
% Initialize the initial value for alpha 
alpha0 = [-0.1 2/3*pi];
% Create Function for alpha - Calculate the derivative and the function
% symbolically
syms alpha beta
f = a1/a2*cos(beta) - a1/a4*cos(alpha)-cos(beta-alpha) + (a1^2 + a2^2 - a3^2 + a4^2)/(2*a2*a4);
f_prime = diff(f, alpha);
% Create the range for the beta angle
beta0 = 0:0.1:2/3*pi;    
%Initialize the vectors alpha_NS and a
a = zeros(length(alpha0),length(beta0));
alpha_NS = zeros(length(alpha0),length(beta0));
% Solve the problem for each value of beta
for j = 1:length(beta0)
    % Define the function and its derivative that are going to be used in the Newton Solver
    f_sub = @(alpha) a1/a2*cos(beta0(j)) - a1/a4*cos(alpha)-cos(beta0(j)-alpha) + (a1^2 + a2^2 - a3^2 + a4^2)/(2*a2*a4);
    f_prime_sub = @(alpha) sin(alpha) + sin(alpha - beta0(j));
    % For each value of alpha0 apply Newton Solver
    for i = 1:length(alpha0)

        alpha = alpha0(i);
        %Use the function Newton Solver created at the end of the script
        [alpha, iter] = NS(alpha,f_sub,f_prime_sub,tol);
        % Save the values obtained from NS into a matrix of [2x21]
        alpha_NS(i,j) = alpha;
        % Use the fzero to check the values obtained from NS 
        a(i,j) = fzero(f_sub,alpha0(i),options);

    end
end
% Calculate the difference between the NS and fzerp
errNSf0 = abs(a - alpha_NS);
%Plot the results of Newton Solver
figure()
set(gcf,'Position',[100 100 800 500])
sgtitle('\textbf{$\alpha$ evolution with respect to $\beta$}','Interpreter','latex')
subplot(1,2,1)
plot(beta0,alpha_NS,'LineWidth',2)
xticks([0 pi/6 pi/3 pi/2 2*pi/3 150*pi/180 pi])
xticklabels([0 30 60 90 120 150 180])
yticks([0 pi/6 pi/3 pi/2 2*pi/3])
yticklabels([0 30 60 90 120])
title('Newton Solver solution')
xlabel('$\beta$ [deg]','Interpreter','latex');
ylabel('$\alpha$ [deg]','Interpreter','latex');
legend({'Initial guess of $\alpha = -0.1$','Initial guess of $\alpha = \frac{2\pi}{3}$'},'Location','northwest','Interpreter','latex','FontSize',10)
grid on;
subplot(1,2,2)
plot(beta0,a,'LineWidth',2)
xticks([0 pi/6 pi/3 pi/2 2*pi/3 150*pi/180 pi])
xticklabels([0 30 60 90 120 150 180])
yticks([0 pi/6 pi/3 pi/2 2*pi/3])
yticklabels([0 30 60 90 120])
title('Fzero solution')
xlabel('$\beta$ [deg]','Interpreter','latex');
ylabel('$\alpha$ [deg]','Interpreter','latex');
legend({'Initial guess of $\alpha = -0.1$','Initial guess of $\alpha = \frac{2\pi}{3}$'},'Interpreter','latex','Location','northwest','FontSize',10)
grid on;
% Plot the error
figure()
set(gcf,'Position',[100 100 1000 500])
subplot(1,2,1)
semilogy(beta0,errNSf0(1,:),'LineWidth',1.5);
title('Difference between NS and Fzero')
xlabel('$\beta$ [rad]','Interpreter','latex');
ylabel('$\alpha$ [rad]','Interpreter','latex');
legend({'Initial guess of $\alpha = -0.1$'},'Location','northwest','Interpreter','latex','FontSize',10)
grid on;
subplot(1,2,2)
plot(beta0,errNSf0(2,:),'LineWidth',1.5);
title('Difference between NS and Fzero')
xlabel('$\beta$ [rad]','Interpreter','latex');
ylabel('$\alpha$ [rad]','Interpreter','latex');
legend({'Initial guess of $\alpha = \frac{2\pi}{3}$'},'Location','northwest','Interpreter','latex','FontSize',10)
grid on;


%% PART III
% Initialize the matrices for Finite Difference
alphaf = zeros(length(alpha0), length(beta0));
alphac = zeros(length(alpha0), length(beta0)); 
alphab = zeros(length(alpha0), length(beta0));

for j = 1:length(beta0)
    % Define function for current beta0(j)
    f_sub = @(alpha) a1/a2 * cos(beta0(j)) - a1/a4 * cos(alpha) - cos(beta0(j) - alpha) + (a1^2 + a2^2 - a3^2 + a4^2) / (2 * a2 * a4);
    
    % Step size (h) for finite difference
    h = 0.1; % Or set h = abs(alpha0(i) * eps);

    for i = 1:length(alpha0)
        % Initialize alpha for each method
        alpha = alpha0(i);

        % Backward difference approximation
        backward_diff = @(alpha_un) (f_sub(alpha_un) - f_sub(alpha_un - h)) / h;

        % Central difference approximation
        central_diff = @(alpha_un)(f_sub(alpha_un + h) - f_sub(alpha_un - h)) / (2 * h);

        % Forward difference approximation
        forward_diff = @(alpha_un) (f_sub(alpha_un + h) - f_sub(alpha_un)) / h;

        % Initialize iteration parameters
        iter = 1;
        err = 1;
        max_iter = 50;
        tol = 1e-5;

        % Newton-Raphson using Backward Difference
        alpha_b = alpha;
        % Calculate until the desired tolerance is reached or the max nr of
        % iteration is reached
        while iter <= max_iter && err > tol
            alpha_new_b = alpha_b - f_sub(alpha_b) / backward_diff(alpha_b);
            % Calculate the error between the last two instances
            err = abs(alpha_new_b - alpha_b);
            % Save the new alpha value obtained
            alpha_b = alpha_new_b;
            iter = iter + 1;
        end
        alphab(i, j) = alpha_b;

        % Reset for central difference iteration
        iter = 1;
        err = 1;

        % Newton-Raphson using Central Difference
        alpha_c = alpha;
        % Calculate until the desired tolerance is reached or the max nr of
        % iteration is reached
        while iter <= max_iter && err > tol
            alpha_new_c = alpha_c - f_sub(alpha_c) / central_diff(alpha_c);
             % Calculate the error between the last two instances
            err = abs(alpha_new_c - alpha_c);
             % Save the new alpha value obtained
            alpha_c = alpha_new_c;
            iter = iter + 1;
        end
        alphac(i, j) = alpha_c;

        % Reset for forward difference iteration
        iter = 1;
        err = 1;

        % Newton-Raphson using Forward Difference
        alpha_f = alpha;
        while iter <= max_iter && err > tol
            alpha_new_f = alpha_f - f_sub(alpha_f) / forward_diff(alpha_f);
            % Calculate the error between the last two instances
            err = abs(alpha_new_f - alpha_f);
             % Save the new alpha value obtained
            alpha_f = alpha_new_f;
            iter = iter + 1;
        end
        alphaf(i, j) = alpha_f;

    end
end

% Plotting the results
figure()
set(gcf,'Position',[100 100 1000 500])
subplot(1,3,1)
plot(beta0, alphab,'LineWidth',2);
xticks([0 pi/6 pi/3 pi/2 2*pi/3 150*pi/180 pi])
xticklabels([0 30 60 90 120 150 180])
yticks([0 pi/6 pi/3 pi/2 2*pi/3])
yticklabels([0 30 60 90 120])
title('\textbf{Backward Difference}','Interpreter','latex')
ylabel('$\alpha$ [deg]','Interpreter','latex')
xlabel('$\beta [deg]$','Interpreter','latex')
legend({'Initial guess of $\alpha = -0.1$','Initial guess of $\alpha = \frac{2\pi}{3}$'},'Location','northwest','Interpreter','latex','FontSize',10)
grid on;
subplot(1,3,2)
plot(beta0, alphac,'LineWidth',2);
xticks([0 pi/6 pi/3 pi/2 2*pi/3 150*pi/180 pi])
xticklabels([0 30 60 90 120 150 180])
yticks([0 pi/6 pi/3 pi/2 2*pi/3])
yticklabels([0 30 60 90 120])
title('\textbf{Central Difference}','Interpreter','latex')
ylabel('$\alpha$ [deg]','Interpreter','latex')
xlabel('$\beta [deg]$','Interpreter','latex')
legend({'Initial guess of $\alpha = -0.1$','Initial guess of $\alpha = \frac{2\pi}{3}$'},'Location','northwest','Interpreter','latex','FontSize',10)
grid on;
subplot(1,3,3)
plot(beta0, alphaf,'LineWidth',2);
title('\textbf{Forward Difference}','Interpreter','latex')
xticks([0 pi/6 pi/3 pi/2 2*pi/3 150*pi/180 pi])
xticklabels([0 30 60 90 120 150 180])
yticks([0 pi/6 pi/3 pi/2 2*pi/3])
yticklabels([0 30 60 90 120])
ylabel('$\alpha$ [deg]','Interpreter','latex')
xlabel('$\beta [deg]$','Interpreter','latex')
legend({'Initial guess of $\alpha = -0.1$','Initial guess of $\alpha = \frac{2\pi}{3}$'},'Location','northwest','Interpreter','latex','FontSize',10)
grid on;

% Plotting the results - error wrt NS
figure()
set(gcf,'Position',[50 100 1350 500])
subplot(1,3,1)
plot(beta0, abs(alphab - alpha_NS),'LineWidth',2);
title('\textbf{Error-Backward Difference}','Interpreter','latex')
ylabel('$\alpha$ [rad]','Interpreter','latex')
xlabel('$\beta [rad]$','Interpreter','latex')
legend({'Initial guess of $\alpha = -0.1$','Initial guess of $\alpha = \frac{2\pi}{3}$'},'Location','northwest','Interpreter','latex','FontSize',10)
grid on;
subplot(1,3,2)
plot(beta0, abs(alphac - alpha_NS),'LineWidth',2);
title('\textbf{Error-Central Difference}','Interpreter','latex')
ylabel('$\alpha$ [rad]','Interpreter','latex')
xlabel('$\beta [rad]$','Interpreter','latex')
legend({'Initial guess of $\alpha = -0.1$','Initial guess of $\alpha = \frac{2\pi}{3}$'},'Location','northwest','Interpreter','latex','FontSize',10)
grid on;
subplot(1,3,3)
plot(beta0, abs(alphaf - alpha_NS),'LineWidth',2);
title('\textbf{Error-Forward Difference}','Interpreter','latex')
ylabel('$\alpha$ [rad]','Interpreter','latex')
xlabel('$\beta [rad]$','Interpreter','latex')
legend({'Initial guess of $\alpha = -0.1$','Initial guess of $\alpha = \frac{2\pi}{3}$'},'Location','northwest','Interpreter','latex','FontSize',10)
grid on;

%% PART V.
% Initialization
clearvars -except f_prime ph alpha0 options beta a1 a2 a3 a4
tol = 1e-5;
beta0 = linspace(0,pi,500);     


for j = 1:length(beta0)    
    % Define the functions
    f_beta_sub = @(alpha) a1/a2*cos(beta0(j)) - a1/a4*cos(alpha)-cos(beta0(j)-alpha) + (a1^2 + a2^2 - a3^2 + a4^2)/(2*a2*a4);
    f_prime_beta_sub = @(alpha) sin(alpha) + sin(alpha - beta0(j));
    % Iterate for each initial value of alpha
        for i = 1:length(alpha0)
            alphav = alpha0(i);
            % Check is the function is constant for the values taken
            if subs(f_prime,beta,beta0(j)) ~= 0 
                % If it is not constant go on and use NS and fzero
                [alphav, iter] = NS(alphav,f_beta_sub,f_prime_beta_sub, tol);
                a(i,j) = fzero(f_beta_sub,alpha0(i),options);
            else
                % Because if it is considered the analytical expression,
                % the derivative is 0 because the f(alpha, pi) = const
                alphah_new = Inf;
                iter = iter + 1;
                % For infinity values of alpha ignor the error(which would be infinity otherwise) 
                %err = NaN; 
                % Display message that the fraction of f/f_prime is infinity
                disp(fprintf('The value of alpha at Beta = %.2f is Infinity and alpha_0 = %.2f\n', beta0(j),alpha0(i)));
            end
            alpha_NS(i,j) = alphav;
        end
end

errNSf0 = abs(a - alpha_NS(:,1:end-1));
% Plot the results: error, fzero and NS
figure()
plot(beta0(1:end-1),errNSf0(1,:),'LineWidth',3.5);
title('\textbf{Error for $\beta$ $\in$ [0 pi]}','Interpreter','latex','FontWeight','bold')
ylim([ -0.3491    2.0944]);
xticks([0 pi/6 pi/3 pi/2 2*pi/3 150*pi/180 pi])
xticklabels([0 30 60 90 120 150 180])
yticks([0 pi/6 pi/3 pi/2 2*pi/3])
yticklabels([0 30 60 90 120])
ylabel('$\alpha [deg]$','Interpreter','latex')
xlabel('$\beta [deg]$','Interpreter','latex')
hold on;
plot(beta0(1:end-1),errNSf0(2,:),'LineStyle','--','LineWidth',3.5);
grid on;
legend({'Initial guess of $\alpha = -0.1$','Initial guess of $\alpha = \frac{2\pi}{3}$'},'Location','northwest','Interpreter','latex','FontSize',10)
hold off;

figure()
plot(beta0,alpha_NS,'LineWidth',1)
ylim([ -0.3491    2.0944]);
xlim([0 pi]);
set(gca,'Clipping','on')
xticks([0 pi/6 pi/3 pi/2 2*pi/3 150*pi/180 pi])
xticklabels([0 30 60 90 120 150 180])
yticks([0 pi/6 pi/3 pi/2 2*pi/3])
yticklabels([0 30 60 90 120])
ylabel('$\alpha [deg]$','Interpreter','latex')
xlabel('$\beta [deg]$','Interpreter','latex')
title('\textbf{Numerical Solution for $\beta$ $\in$ [0 pi]}','Interpreter','latex','FontWeight','bold')
grid on;
legend({'Initial guess of $\alpha = -0.1$','Initial guess of $\alpha = \frac{2\pi}{3}$'},'Location','northwest','Interpreter','latex','FontSize',10)

figure()
ph = ph + 1;
plot(beta0(1:end-1), a,'LineWidth',2)
ylim([ -0.3491    2.0944]);
xlim([0 pi]);
set(gca,'ClippingStyle','rectangle')
xticks([0 pi/6 pi/3 pi/2 2*pi/3 150*pi/180 pi])
xticklabels([0 30 60 90 120 150 180])
yticks([0 pi/6 pi/3 pi/2 2*pi/3])
yticklabels([0 30 60 90 120])
title('\textbf{Analytical Solution for $\beta$ $\in$ [0 pi]}','Interpreter','latex','FontWeight','bold')
grid on;
legend({'Initial guess of $\alpha = -0.1$','Initial guess of $\alpha = \frac{2\pi}{3}$'},'Location','northwest','Interpreter','latex','FontSize',10)
ylabel('$\alpha [deg]$','Interpreter','latex')
xlabel('$\beta [deg]$','Interpreter','latex')

%% Exercise 2 

clearvars;
close all;
clc;
%% 2.1 Heun's method (RK2)

% Initialize Values
data.m0 = 20;                % mass[kg]
data.cm = 0.1;               % mass coefficient [kg/s]
data.f_t = 1;                % f(t) - is a function [N]
data.alpha = 0.01;           % velocity coefficient [Ns/m]
data.F = 0;                  % Drag Force [N]
% Heun's method 
% Define the function  
f = @(v,t) (data.F + data.f_t - data.alpha*v)/(data.m0 - data.cm*t);  % acceleration function [m/s^2]

% Initialize the time steps
h = [50,20,10,1];       % given step sizes
v_RK2 = zeros(200,length(h));     % initialize velocity vector for RK2
v_an = zeros(200,length(h));      % initialize velocity vector foa analytical
tic                               % begin time measurement of RK2
for j = 1:length(h)
    % Define the time interval for each step size
    t_HEUN = 0:h(j):160;
    v = zeros(ceil(160/h(j)),1);
    % Initial condition of velocity
    v(1) = 0;
    % Heun's method (equation)
    for i = 1:(length(t_HEUN)-1)
        v_pred = v(i) + 1 * h(j) * f(v(i),t_HEUN(i)); % Predictor
        v(i+1) = v(i) + h(j) * (1/2 * f(v(i),t_HEUN(i)) + 1/2 * f(v_pred,t_HEUN(i + 1)));       % Corrector
    end
    % Plot the results obtained for each step size
    figure(1)
    hold on;
    subplot(2,2,j)
    plot(t_HEUN,v,'b','LineWidth',2,'DisplayName','RK2 Numerical')
    title(sprintf('\\textbf{RK2 method for h =} %.0f', h(j)),'Interpreter','latex');  
    xlabel('Time [s]')
    ylabel('Velocity v [m/s]')
    grid on;
    % Save the velocity obtained through RK2 for future purposes
    v_RK2(1:length(v),j) = v; 
    
end
sgtitle('\textbf{Analytical vs Numerical for RK2}','Interpreter','latex');
hold off;
t_RK2 = toc;                    % save the GPU time needed by RK2 in variable t_RK2

% Compare numerical vs analytical solution
% Begin analytical
clear v dvdt
% Set performance characteristic for ode
options = odeset('RelTol',10^-9,'AbsTol',10^-9);
% Use integrator as analytical method for the overall interval [0 160], ode
% sets its one stepsize
[t1,v] = ode78( @(t,v) (data.F + data.f_t - data.alpha*v)/(data.m0 - data.cm*t),[0 160],0,options);
% Calculate the analytical solution using ode for the RK2 step sizes 
for j = 1:length(h)
    % Define the time interval wrt stepsize for ode
    t_ode = 0:h(j):160;
    % Ignore exit time it is equal to input time interval and calculate for
    % each time interval
    [~,v_eval] = ode78(@(t,v) (data.F + data.f_t - data.alpha*v)/(data.m0 - data.cm*t),t_ode,0,options);
    % Plot the analytical results
    figure(1)
    set(gcf,'Position',[100 100 800 500])
    hold on;
    subplot(2,2,j)
    plot(t1,v,'Color',"#D95319",'LineStyle',':','LineWidth',2,'DisplayName','Analytical ')
    legend('Location', 'best','FontSize',8);
    % Calculate the analytical solution using the given formula for
    % velocity
    v_sol = 1/data.alpha - (1/data.alpha)*(ones(1,length(t_ode)) - data.cm .* t_ode ./ data.m0) .^ (data.alpha / data.cm); 
    % Plot the solution obtained from the equation
    plot(t_ode,v_sol,'Color',"g",'LineStyle','--','LineWidth',2,'DisplayName','Eq. Analytical')
    % Calculate the difference between the RK2 and analytical solution
    err = abs(v_RK2(1:length(v_eval),j) - v_sol');
    % Plot the error
    figure(4)
    set(gcf,'Position',[100 100 800 500])
    sgtitle('\textbf{Error between numerical and analytical}','Interpreter','latex')
    hold on;
    subplot(2,2,j)
    semilogy(t_ode,err,'LineWidth',2,'DisplayName','Error RK2')
    title(sprintf('\\textbf{Error for h =} %.0f', h(j)),'Interpreter','latex');
    grid on;
    xlabel('Time[s]')
    ylabel('Velocity [m/s]')
end
hold off;
% Plot the difference betweent the ode and given formula for analytial
% method
v_sol1 = 1/data.alpha - (1/data.alpha)*(ones(1,length(t1)) - data.cm .* t1 ./ data.m0) .^ (data.alpha / data.cm); 
err2 = abs(v_sol1 - v);
figure(12)
plot(t1,err2,'LineWidth',2)
title('\textbf{Error between Exact solution and ode78}','Interpreter','latex')
xlabel('Time[s]')
ylabel('$\Delta$ v [m/s]','Interpreter','latex')
grid on;


%% Point 3 RK4
clearvars -except v_RK2 t_RK2 options data f

% RK4 methods
% Define the function
% f =@(v,t) (F + f_t - alpha*v)/(m0 - cm*t);
% Define the step size
h = [50,20,10,1];
% Initialize the RK4 velocity vector
v_RK4 = zeros(165,length(h));      
% Begin timming the RK4 method
tic
for j = 1:length(h)
    % Define the time vector for each step size
    t = 0:h(j):160;
    % Initialize the velocity vector
    v = zeros(length(t),1);
    % Set the initial condition
    v(1) = 0;
    % Begin the RK4 method
    for i = 1:length(t) - 1
        vp1 = v(i) + 1/2*h(j)*f(v(i),t(i));             % First Predictor
        vp2 = v(i) + 1/2*h(j)*f(vp1,t(i) + 1/2*h(j));   % Second Predictor
        vp3 = v(i) + 1*h(j)*f(vp2,t(i) + 1/2*h(j));     % Third Predictor
        % Corrector
        v(i+1) = v(i) + 1*h(j)*(1/6 * f(v(i),t(i)) + 1/3*f(vp1, t(i) + 1/2*h(j)) + 1/3*f(vp2,t(i) + 1/2*h(j)) + 1/6*f(vp3,t(i+1)));
    end
    % Plot the results obtained from the RK4
    figure(3)
    set(gcf,'Position',[100 100 800 500])
    hold on;
    subplot(2,2,j)
    plot(t,v,'b','LineWidth',2,'DisplayName','RK4 Numerical')
    xlabel('Time[s]')
    ylabel('Velocity[m/s]')
    grid on;
    title(sprintf('\\textbf{RK4 method for h =} %.0f', h(j)),'Interpreter','latex');
    % Save the RK4 results for future use
    v_RK4(1:length(v),j) = v;
end
sgtitle('\textbf{Analytical vs Numerical for RK4}','Interpreter','latex');
hold off;
t_RK4 = toc;                % Save the RK4 computational Time in t_RK4
% Compare numerical vs analytical solution
% Analytical Solution
[t1,v_ode] = ode45( @(t,v) (data.F + data.f_t - data.alpha*v)/(data.m0 - data.cm*t),[0 160],0,options); 
clear v 
% Calculate for each step sie the analytical solution
for j = 1:length(h)
    % Initialize the time interval for each step size
    t_ode = 0:h(j):160;
    % Calculate the analytical solution using ode for each time interval 
    [~,v] = ode45( @(t,v) (data.F + data.f_t - data.alpha*v)/(data.m0 - data.cm*t),t_ode,0,options); 
    % Plot the results
    figure(3)
    hold on;
    subplot(2,2,j)
    plot(t_ode,v,'Color',"g",'LineStyle','--','LineWidth',2,'DisplayName','Analytical for Tsteps')
    % Plot the solution over entire interval
    plot(t1,v_ode,'Color',"#D95319",'LineStyle',':','LineWidth',2,'DisplayName','Analytical for T = 0:160')
    legend('Location', 'northwest');
    % Calculate the velocity using the analytical formula
    v_sol = 1/data.alpha - (1/data.alpha)*(ones(1,length(t_ode)) - data.cm .* t_ode ./ data.m0) .^ (data.alpha / data.cm); 
    % Calculate the error between the analytical solution and RK4
    err = abs(v_RK4(1:length(v_sol),j) - v_sol');
    % Plot the results
    figure(4)
    hold on;
    subplot(2,2,j)
    hold on;
    semilogy(t_ode,err,'LineWidth',2,'Color',"r",'DisplayName','Error RK4','LineStyle','-')
    title(sprintf('\\textbf{Error for h =} %.0f', h(j)),'Interpreter','latex');
    grid on;
    xlabel('Time[s]')
    ylabel('Velocity [m/s]')
    legend('Location','best')
end
hold off;

for j = 1:length(h)
    t_err = 0:h(j):160;
    % Calculate the difference between the performance of RK2 and RK4
    err =  abs(v_RK4(1:length(t_err),j) -  v_RK2(1:length(t_err),j));
    figure(5)
    set(gcf,'Position',[100 100 800 500])
    sgtitle('\textbf{Difference Between RK2 and RK4}','Interpreter','latex');
    subplot(2,2,j)
    plot(t_err,err,'LineWidth',2)
    ylabel('$\Delta V$ [m/s]', 'Interpreter','latex')
    xlabel('Time [s]')
    xlim([0 160])
    grid on;
    hold on;
end
hold off;


%% Part II

clear dvdt f v t h m0 f_t alpha

% Initialize Values
data.rho = 900;              % Density [kg/m^3]
data.Cd = 2.05;              % Drag Coefficient
data.Am = 1;                 % Object Area [m^2]

% Define the Drag force
F = @(v) -0.5 * data.rho * data.Cd * data.Am * v^2;
m_t = @(t) data.m0 - data.cm * t;
f = @(t,v) (F(v) + data.f_t - data.alpha * v)/m_t(t);
% Choose two step sizes
h_vec = [1 0.1];
% Calculate for both step sizes
for j = 1:2
    h = h_vec(j);
    % Define the time interval for each step size
    t_HEUN = 0:h:160;
    % Initialize the velocity vector
    v = zeros(length(h),1);
    % Set the initial condition
    v(1) = 0;
    % Begin RK2
    for i = 1:(length(t_HEUN)-1)
        v_pred = v(i) + 1 * h * f(t_HEUN(i),v(i)); % Predictor
        v(i+1) = v(i) + h * (1/2 * f(t_HEUN(i),v(i)) + 1/2 * f(t_HEUN(i + 1),v_pred)); % Corrector
    end
    % Calculate the analytical solution using ode for the time interval
    [t1,v_ODE] = ode45(f,t_HEUN,0,options); 
    % Plot the results
    figure(15)
    set(gcf,'Position',[100 100 800 500])
    subplot(1,2,j)
    hold on;
    plot(t_HEUN,v,'Color',"#77AC30",'DisplayName','RK2 Method $\rho = 900$','LineWidth',2.5)
    grid on; 
    plot(t1,v_ODE,'k','DisplayName','Analytical Method $\rho = 900$','LineWidth',2.5,'LineStyle',':')
    legend('Location','best','Interpreter','latex')
    title(sprintf('\\textbf{RK2 method for h = }%.1f',h),'Interpreter','latex')
    ylabel('$\Delta V$ [m/s]', 'Interpreter','latex')
    xlabel('Time [s]')
    hold off;
    % Clear the vector for second iteration
    v = [];
end

%% Exercise 3

clearvars;
close all;
clc;
%% Point 1 RK2
% Define range of alpha and initial guesses for h
alpha0 = linspace(pi, 0, 160);
h0 = 2;

% Preallocate results matrix
h_sol = NaN(length(alpha0), 1);  % Use NaN to ignore invalid solutions directly

% Loop over initial guesses for h

    % Loop over all alpha values
    for i = 1:length(alpha0)
        % Define the matrix A
        A = [0, 1; -1, 2 * cos(alpha0(i))];
        
        % Compute the matrix FRK2 as a function of h
        FRK2_func = @(h) eye(2) + h * A + (h^2 / 2) * (A^2);
        
        % Define the objective function to satisfy |eig(FRK2)| = 1 condition
        stability_func = @(h) max(abs(eig(FRK2_func(h)))) - 1;
        
        % Solve for h numerically using fzero (faster than fsolve for scalar functions)
        try
            h_solution = fsolve(stability_func, h0);
            if h_solution > 0  % Ensure h is positive
                h_sol(i) = h_solution;
            end
        catch
            h_sol(i) = NaN;  % If no solution, keep NaN
        end
        h0 = h_sol(i);
    end

% Pre-allocate the error check vector
f_check = NaN(length(alpha0), 1);

% Loop through each alpha to check if solutions meet the error tolerance
for i = 1:length(alpha0)
    % Define the matrix A for the current alpha
    A = [0 1; -1 2 * cos(alpha0(i))];
    
    % Get the last column of h_sol for the current alpha, which we want to check
    h_val = h_sol(i, end);
    
    % Check only if h_val is not NaN
    if ~isnan(h_val)
        % Compute FRK2 directly with h_val (no symbolic operations)
        FRK2 = eye(2) + h_val * A + (h_val^2 / 2) * A^2;
        
        % Compute eigenvalues and check maximum absolute value
        max_eigenvalue = max(abs(eig(FRK2)));
        
        % Calculate the error relative to 1 (for stability boundary)
        f_check(i) = max_eigenvalue - 1;
    else
        f_check(i) = NaN;  % If h is NaN, keep error as NaN
    end
end


% Plot the Eigenvalues in the hlambda plane
alpha_lambda1 = zeros(1,length(alpha0));
alpha_lambda2 = zeros(1,length(alpha0));
for i = 1:length(alpha0)
    A = [0 1; -1 2 * cos(alpha0(i))];
    alpha_lambda = h_sol(i,end) * eig(A);
    alpha_lambda1(i) = alpha_lambda(1);
    alpha_lambda2(i) = alpha_lambda(2); 
end
figure(1);
hold on;
plot(real(alpha_lambda1),imag(alpha_lambda1),'-','Color',"#0072BD",'DisplayName','RK2 stability region margin','LineWidth',1.5);
plrk2 = plot(real(alpha_lambda2),imag(alpha_lambda2),'-','Color',"#0072BD","LineWidth",1.5);
plrk2.Annotation.LegendInformation.IconDisplayStyle = 'off';
grid on;
title('Stability Regions')
legend('Location','best')
xlabel('Re($\lambda$h)','Interpreter','latex');
ylabel('Im($\lambda$h)','Interpreter','latex');
axis([-3 0.5 -4 4])


%% Point 1 RK4

% RK4 stability region that works this is gold
close all;
clearvars;

% Define range of alpha and initial guess for h
alpha0 = linspace(pi, 0, 1400);
h0 = 2; % Initial guess for h
options = optimset('TolX', 1e-6, 'TolFun', 1e-6, 'Display', 'off'); % fsolve options

A = @(aa) [0, 1; -1, 2 * cos(aa)]; % The matrix A(aa)
% RK4 stability function
FRK4_func = @(h, aa) eye(2) + h * A(aa) + h ^ 2  / 2 * A(aa)^2 + h ^3 / 6 * A(aa) ^ 3+ h ^ 4 / 24 * A(aa) ^ 4;
h_solution = zeros(1,length(alpha0));
% Loop over all alpha values to solve for h using fsolve
for i = 1:length(alpha0)
    h_solution(i) = fsolve(@(h) max(abs(eig(FRK4_func(h, alpha0(i))))) - 1, h0, options);
    h0 = h_solution(i);  % Update h0 for next iteration
    a(i) = alpha0(i);
    if h0 < 2.2  % Stop when h becomes very small
        break;
    end
end

% Reverse the loop to find stability for smaller alpha values
ii = i - 1;

j = 0;
h0 = h_solution(i);
for i = ii:-1:1
    h_solution(ii + j) = fsolve(@(h) max(abs(eig(FRK4_func(h, alpha0(i))))) - 1, h0, options);
    h0 = h_solution(ii + j);   
    a(ii + j) = alpha0(i);
    j = j + 1;
    if h0 < 0 
        break;
    end
end
a = a(a~=0);
h_solution = h_solution(h_solution~=0);
% Calculate the lambda*h for each alpha
for i = 1:length(a)
    alpha_lambda = h_solution(i) * eig(A(a(i))); % Compute the eigenvalue product
    alpha_lambda1(1,i) = alpha_lambda(1);  % First eigenvalue
    alpha_lambda2(1,i) = alpha_lambda(2);  % Second eigenvalue
end

% Plot the results: RK4 stability region
figure(1);
hold on;
plot(real(alpha_lambda1),imag(alpha_lambda1),'o-','Color',"#D95319",'DisplayName','RK4 stability region margin')
plrk4 = plot(real(alpha_lambda2),imag(alpha_lambda2),'o-','Color',"#D95319");
plrk4.Annotation.LegendInformation.IconDisplayStyle = 'off';
xlabel('Re($\lambda h$)','Interpreter','latex');
ylabel('Im($\lambda h)$','Interpreter','latex');
title('RK4 Stability Region');

grid on;

f_check = NaN(length(a), 1);

% Loop through each alpha to check if solutions meet the error tolerance
    for i = 1:length(a)
        
        % Get the last column of h_sol for the current alpha, which we want to check
        h_val = h_solution(i);
        
        % Check only if h_val is not NaN
        if ~isnan(h_val)
          
            % Compute eigenvalues and check maximum absolute value
            max_eigenvalue = max(abs(eig(FRK4_func(h_val,a(i)))));
            
            % Calculate the error relative to 1 (for stability boundary)
            f_check(i) = max_eigenvalue - 1;
        else
            f_check(i) = NaN;  % If h is NaN, keep error as NaN
        end
    end

% Optional: Evaluate if all errors are below 10^-6
is_within_tolerance = all(abs(f_check(~isnan(f_check))) < 1e-4);

% Point 1) RK4 for pi
h_sol_pi = h_solution(1);
% Point 2) RK4 for [pi 0]
h_sol_pi_0 = h_solution;

%% Point 2 BI20.3

clc;
clearvars;

%implement also here numerical integration
tic
options = optimset('TolX',1e-6,'TolFun',1e-6,'Display','off');
alpha0 = linspace(0,pi,160);
h0_0 = [5 5 10 -10 -5];
theta = [0.2, 0.3, 0.4, 0.6, 0.8];
col = ["#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30"];	

h = zeros(1,length(alpha0));
for k = 1:length(theta)
     h0 = h0_0(k);
     clear hsol;
     hsol = zeros(length(alpha0));
for i = 1:length(alpha0)        

        Abi = [0, 1; -1, 2 * cos(alpha0(i))];
        % Compute the matrix BI2cl as a function of h
        BI2_func = @(h) (eye(2) - (1 - theta(k)) * h * Abi + ( Abi * h * (1 - theta(k)))^2/2) \ (eye(2) + theta(k) * h * Abi + (Abi * theta(k) * h) ^ 2 / 2);
        
        % Define the objective function to satisfy |eig(FRK2)| = 1 condition
        stability_func = @(h) max(abs(eig(BI2_func(h)))) - 1;
        hsol(i) = fsolve(stability_func, h0,options);
        h0 = hsol(i);
end
h = hsol;
lambalpha = zeros(2,length(alpha0));
for i = 1:length(alpha0)
    Abi = [0 1; -1 2 * cos(alpha0(i))];
    alpha_lambda = h(i) * eig(Abi);
    lambalpha(:,i) = alpha_lambda;
    alpha_lambda1(i) = alpha_lambda(1);
    alpha_lambda2(i) = alpha_lambda(2);
end

figure(3)
hold on;
plot(real(alpha_lambda1),imag(alpha_lambda1),'-','Color',col(k),'LineWidth',1.5,'DisplayName',sprintf('\\theta = %.1f',theta(k)));
n1 = plot(real(alpha_lambda2),imag(alpha_lambda2),'-','Color',col(k),'LineWidth',1.5);
n1.Annotation.LegendInformation.IconDisplayStyle = 'off';
text(real(alpha_lambda1(1)) + 0.1,imag(alpha_lambda1(1)),string(real(alpha_lambda1(1))))
n2 = plot(real(alpha_lambda1(1)),imag(alpha_lambda1(1)),'o','MarkerFaceColor','b','DisplayName','');
n2.Annotation.LegendInformation.IconDisplayStyle = 'off';
grid on;
end
legend('Location','best')
title('Stability Regions for $BI2_{0.3}$','Interpreter','latex');
xlabel('Re($\lambda$h)','Interpreter','latex');
ylabel('Im($\lambda$h)','Interpreter','latex');
axis([-12 12 -8 8])
hold off;
time = toc;

%% Exercise 4 
clearvars;
close all;
clc;

%% 4.2) RK2

Kc = 0.0042;        % Convective heat loss coeff [J/(sK)]
Kr = 6.15*10^-11;   % Radiation Heat loss coeff [J/(sK^4)]
T_alpha = 277;      % Surrounding air temperature [K]
C = 45;             % Mass Thermal capacity [J/K]
T_0 = 555;          % Initial Temperature [K]
h.RK2 = 720;        % RK2 step size
h.RK4 = 1440;     % RK4 step size

% 1 - lumped capacitance modeling
% Mathematical model of the system
dTdt = @(t,T) -Kc/C * (T - T_alpha) - Kr/C * (T^4 - T_alpha^4);

% RK2
% Initialize the temperature vector and time vector
T = zeros(1,100);
tt_rk2 = zeros(1,length(T));
% Initialize the initial condtion
T(1) = T_0;             % Initial temperature of the body
i = 1;                  
t = 0;                  % Initial time
err = 1;                % Set the initial error, which shall reach below limit
                        % (10^-5) in order to account for the equilibrium 

% Use RK2 until the equilibrium is reached
while T(i) > T_alpha && abs(err) > 10^-5
    
    T_k1_p = T(i) + h.RK2 * dTdt(t,T(i));           % Predictor
    T (i + 1) = T(i) + h.RK2/2 * (dTdt(t,T(i)) + dTdt(t + h.RK2, T_k1_p)); % Corrector
    % Update the time
    t = t + h.RK2;
    % Update the error
    err = T(i) - T(i+1);

    i = i + 1;
    % Save the time in a vector
    tt_rk2(i) = t;
end
% Eliminate the zero elements from the time vector 
tt_rk2 = [0 tt_rk2(tt_rk2 ~= 0)];
% Eliminate the zero elements from the temperature vector
T_rk2 = T(T ~= 0);

% RK4
% Initialize the vectors
T = zeros(1, 100);              % Temperature vector
tt_rk4 = zeros(1,100);          % Time vector
% Set the initial conditions
T(1) = T_0;                     % Initial Temperature
i = 1;
t = 0;                          % Initial Time
err = 1;                        % Set the initial error, which shall reach below limit
                                % (10^-5) in order to account for the equilibrium 
% Use RK4 until you reach equillibrium
while T(i) > T_alpha && abs(err) > 10^-5

    k1 = dTdt(t,T(i));          % First predictor
    k2 = dTdt(t + h.RK4/2, T(i) + h.RK4/2 * k1);        % Second predictor
    k3 = dTdt(t + h.RK4/2, T(i) + h.RK4/2 * k2);        % Third predictor
    k4 = dTdt(t + h.RK4, T(i) + h.RK4 * k3);            % Forth Predictor
    T(i + 1) = T(i) + h.RK4/6 * (k1 + 2 * k2 + 2 * k3 + k4);    % Corrector
       
    t = t + h.RK4;              % Update the time
    err = T(i) - T(i+1);        % Update the error
    i = i + 1;
    tt_rk4(i) = t;              % Save time step
end
% Eliminate the zero elements from the time and temperature vectors
tt_rk4 = [0 tt_rk4(tt_rk4 ~= 0)];
T_rk4 = T(T ~= 0);
% Set the Tolerances for the ode integrator
options = odeset('AbsTol',10^-5,'RelTol',10^-5,'Events',@(t,T) EventFct4(t,T,1));
% Calculate the analytical solution using ode45
[tt, TT] = ode45(dTdt,[0 t],T_0,options);
% Plot the results for RK2, RK4 and analytical solution
figure(1)
hold on;
plot(tt,TT,'LineWidth',3,'LineStyle','-','DisplayName','Analytical solution')
plot(tt_rk2,T_rk2,'LineWidth',3,'LineStyle',':','DisplayName','RK2 solution');
plot(tt_rk4,T_rk4,'LineWidth',3,'LineStyle','--','DisplayName','RK4 solution')
xlabel('Time[s]')
ylabel('Temperature[K]')
title('Temperature evolution with time')
legend('Location','best')
grid on;
hold off;
% zoom on the initial part of the graph where the error is larger
axes('position',[.65 .475 .25 .25])
box on
plot(tt,TT,'LineWidth',3)
hold on;
plot(tt_rk2,T_rk2,'LineWidth',3,'LineStyle',':');
plot(tt_rk4,T_rk4,'LineWidth',3,'LineStyle','--')
axis([0 1500  420 555]);
grid on;

%% Plot the errors
% Calculate the errors
% As the RK2 and RK4 have different time steps they need to be interpolated
T_rk4_interp = interp1(tt_rk4, T_rk4, tt_rk2, 'linear');
% Calculate the error between RK2 and RK4
err_rk2_rk4 = abs(T_rk4_interp - T_rk2);
% Calculate the analytical solution for the RK2 time interval
options = odeset('AbsTol',10^-5,'RelTol',10^-5,'Events',@(t,T) EventFct4(t,T,1));
[tt_1, TT_rk2] = ode45(dTdt,tt_rk2,T_0,options);
% Calculate the error betweent RK2 and analytical
err_rk2 = abs(TT_rk2 - T_rk2');
% Calculate the analytical solution for the RK4 time interval
options = odeset('AbsTol',10^-5,'RelTol',10^-5,'Events',@(t,T) EventFct4(t,T,1));
[tt_2, TT_rk4] = ode45(dTdt,tt_rk4,T_0,options);
% Calculate the error betweent RK4 and analytical
err_rk4 = abs(TT_rk4 - T_rk4');
% Plot the results
figure()
semilogy(tt_1,err_rk2,'LineWidth',1.75)
xlabel('Time[s]', 'FontSize',14)
ylabel('Temperature[K]', 'FontSize',14)
title('Error between RK2 and analytical', 'FontSize',16)
grid on;
figure
semilogy(tt_2,err_rk4,'LineWidth',1.75)
title('Error between RK4 and analytical', 'FontSize',16)
xlabel('Time[s]', 'FontSize',14)
ylabel('Temperature[K]', 'FontSize',14)
grid on;
figure
semilogy(tt_rk2,err_rk2_rk4,'LineWidth',1.75)
title('Error between RK2 and RK4', 'FontSize',16)
xlabel('Time[s]', 'FontSize',14)
ylabel('Temperature[K]', 'FontSize',14)
grid on;

%% Exercise 5
clearvars;
close all;
clc;
%% Exercise 5 1)

%% Exercise 5 2)
% Initialize the data of the circuit

R = 25;     % resistance in ohms
L = 20e-3;  % inductance in Henry
C = 200e-3; % capacitance in Farad

% Create the state matrix of the system
A = [0 1; -1/(L * C) -R/L];
% Calculate the eignevalues
lambda = eig(A);
% Check the stifness of the system
if abs(lambda(2)) * 10^-3 >= abs(lambda(1)) && lambda(2) <= lambda(1)
    disp('The system is stiff')
else
    disp('The system is not stiff')
end
% Consider a step size that is not in stable condition the maximum for RK4
h_max = -2.8/lambda(2);

h = h_max;

% Eigenvalues in (hλ) domain
h_lambda1 = h * lambda(1);
h_lambda2 = h * lambda(2);

% Plot stability regions (RK2)
figure(1);
hold on;
% Plot RK2 stability region
real_part = -4:0.01:4;
imag_part = -4:0.01:4;
[Re, Im] = meshgrid(real_part, imag_part);
z = Re + 1i * Im;
G_RK2 = 1 + z + (z.^2) / 2;     % From the system
contourf(Re, Im, abs(G_RK2) <= 1, [1 1]);
colormap([0.99 0.99 0.99; 0.96 1 0.95]);

% Plot eigenvalues in the (hλ) plane
plot(real(h_lambda1), imag(h_lambda1), 'r', 'Marker', 'o','MarkerFaceColor','r');
plot(real(h_lambda2), imag(h_lambda2), 'b', 'Marker','o','MarkerFaceColor','b');
xlabel('Re(h\lambda)');
ylabel('Im(h\lambda)');
title('Eigenvalues in the (h\lambda)-Plane for RK2 Stability Domain');
legend({'RK2 stability domain','h\lambda_1', 'h\lambda_2'});
grid on;
axis equal;

% Implicit Euler Stability Function (as an approximation)
% For implicit Euler, R(z) = 1 / (1 - z)
Rz = 1 ./ (1 - z);      % Stability of BE B = (I - h * A)^-1;

% Plot the stability region for |R(z)| <= 1
figure(2);
hold on;
contourf(Re, Im, abs(Rz) <= 1, [1, 1],'DisplayName','IEX4 Stability boundary'); % Stability boundary
map = [0.99 0.99 0.99; 0.96 1 0.95];
colormap(map); 
title('Approximate Stability Region for IEX4');
xlabel('Re(h$\lambda$)','Interpreter','latex');
ylabel('Im(h$\lambda$)','Interpreter','latex');
plot(real(h_lambda1), imag(h_lambda1), 'r', 'Marker', 'o','MarkerFaceColor','r','DisplayName','h$\lambda _1$');
plot(real(h_lambda2), imag(h_lambda2), 'b', 'Marker','o','MarkerFaceColor','b','DisplayName','h$\lambda _2$');
legend('Interpreter','latex');
grid on;
axis equal;
hold off;

%% IEX4 

clear x t;% Initial conditions
v0 = 12; % Initial voltage in volts
x = [v0*C; 0]; % Initial state [charge (q), current (dq/dt)]
h = 0.00002; % Step size
t_0 = 0;
t_end = 20; % Final time
num_steps = floor(t_end / h);

% Define the function handle for the matrix-vector multiplication
f = @(x) A * x; % Only dependent on x, t is unused
tic
% Initialize time vector and solution storage
x_vals = zeros(num_steps, 2); % Preallocate state storage
x_vals(1,:) = [v0*C, 0]; % Start with initial state
% Time-stepping loop
t = linspace(0,t_end,num_steps);
for j = 1:num_steps - 1
       % The above equations are equivallent to the k = fsolve(k - x -
       % h*f(k,t)) but they are more computational efficient
      k1 = (eye(2) - h*A)\x;
      k2a = (eye(2) - h/2*A)\x;
      k2 = (eye(2) - h/2*A)\k2a;
      k3a = (eye(2) - h/3*A)\x;
      k3b = (eye(2) - h/3*A)\k3a;
      k3 = (eye(2) - h/3*A)\k3b;
      k4a = (eye(2) - h/4*A)\x;
      k4b = (eye(2) - h/4*A)\k4a;
      k4c = (eye(2) - h/4*A)\k4b;
      k4 = (eye(2) - h/4*A)\k4c;

    % Update state with a weighted combination of k1, k2, k3, and k4
    x = (-1/6 * k1 + 4 * k2 - 27/2 * k3 + 32/3 * k4);
    
    x_vals(j + 1, :) = x';

end
t_vals = t;
% Check the computational efficiency of the method
time_IEX4 = toc;
% Implement an ODE for checking the solution of the IEX4
options = odeset('AbsTol',1e-5,'RelTol',1e-5);
[tt,xx] = ode15s(@(t,x) A * x,t_vals,[v0*C; 0],options);
err_tran = abs(xx - x_vals);
% Plot the Results
figure(3)
plot(t_vals,x_vals,'LineWidth',1.5)
title('Transient respons through the IEX4')
xlabel('Time[s]')
ylabel('x')
legend('$x_1 = q$','$x_2 = \frac{dq}{dt}$','Interpreter','latex')
grid on;
figure(4)
semilogy(tt,err_tran,'LineWidth',1.5)
title('Error between analytical and numerical')
xlabel('Time[s]')
ylabel('x')
legend('$q$ $error$','$\frac{dq}{dt}$ $error$','Interpreter','latex','Location','best')
grid on;

figure(13)
plot(tt,xx,'LineWidth',1.5)
title('Transient respons through the IEX4')
xlabel('Time[s]')
ylabel('x')
legend('$x_1 = q$','$x_2 = \frac{dq}{dt}$','Interpreter','latex')
grid on;

%% Check RK2 stability and maximum h
clearvars -except h x_vals; clc;

clear h;
% Define parameters
R = 25;     % resistance in ohms
L = 20e-3;  % inductance in Henry
C = 200e-3; % capacitance in Farad
v0 = 12;    % initial voltage in volts

% System matrix for state-space representation
A = [0, 1; -1/(L*C), -R/L];
lambda = eig(A);
% Initial conditions
x0 = [v0*C; 0];
t_end = 20; % end time
% Implement the initial guess for step size
h0 = 0.02;

optionsfsolve = optimset('TolFun', 1e-6, 'Display', 'off'); % fsolve options
% Define the function of RK2
stability_func = @(h) max(abs(eig(eye(2) + h*(A) + (h^2/2)*A^2))) - 1;
stab = @(h) max(max(abs(eye(2) + h*eig(A) + (h^2/2)*eig(A).^2))) - 1;
% Use fsolve to find the solution
h_max =  fsolve(stability_func, h0,optionsfsolve);
h_max1 = fsolve(stab, h0,optionsfsolve);

% First approach
h_values = [0.01, 0.005,h_max + 1e-5, h_max, h_max - 1e-5, 0.001]; % try different step sizes for stability testing
tic
% Loop over step sizes
n = t_end/h_values(end);
t_vals = zeros(n,1);
x_vals = zeros (n,2); 
for i = 1:length(h_values)
    h = h_values(i);
    % Initialize variables
    t = 0;
    x = x0;
    
    t_vals(1) = t;
    x_vals(1,:) = x';

    j = 2;
    % RK2 Simulation
    while t < t_end
        % RK2 steps
        k1 = A * x;
        k2 = A * (x + h * k1 / 2);
        x = x + h * k2; % RK2 update step

        % Store values for plotting
        t = t + h;
        t_vals(j) = t;
        x_vals(j, :) = x';
        j = j + 1;
    end
    t_vals = nonzeros(t_vals);
    x_vals = x_vals(1:length(t_vals),:);
    % Plot results for this step size
    figure(5)
    subplot(2,3,i)
    sgtitle('Stability Analysis RK2')
    plot(t_vals, x_vals(:,1), 'DisplayName', '$q$','LineWidth',1.5);
    hold on;
    plot(t_vals, x_vals(:,2), 'DisplayName', '$\frac{dq}{dt}$','LineWidth',1.5);
    xlabel('Time [s]');
    ylabel('State Variables');
    title([' h = ', num2str(h)]);
    legend('Interpreter','latex','Location','best');
    grid on;
    hold off;
end
time_RK2 = toc;
% Stability and Accuracy Checks
% Compute FRK2(h, alpha) and check eigenvalues
for h = h_values
    FRK2 = eye(2) + h * A + (h^2 / 2) * (A^2);
    eig_FRK2 = eig(FRK2);
    disp(['Eigenvalues of FRK2 with h = ', num2str(h), ':']);
    disp(eig_FRK2);
    disp(['Magnitude of eigenvalues: ', num2str(abs(eig_FRK2)')]);

    % Check if eigenvalues' magnitudes are <= 1
    if all(abs(eig_FRK2) <= 1)
        disp(['The solution with h = ', num2str(h), ' is stable.', newline]);
    else
        disp(['The solution with h = ', num2str(h), ' is unstable.', newline]);
    end
end

%% Plot the solution of Causal (Matlab) - already ploted above together to the solution of acausal
% model, Simscape.
% Define the Simscape model name
modelName = 'Martin10903720_Assign1_EX5'; % Replace with your model name
% Load the Simulink model
load_system(modelName);

% Run the simulation
out = sim(modelName);
% Plot the result from Simulink
figure()
plot(out.ex5{1}.Values.Time,out.ex5{1}.Values.Data,'LineWidth',2.5)
xlabel('Time[s]')
ylabel('I[A]')
ylim([-0.5 0.02])
title('Initial Circuit Current Simulation')
grid on;

%% Exercise 6 
clearvars;
close all;
clc;
%% Part 1)
% Set the initial condition of the ball
x(1) = 10;                  % Height [m]
v(1) = 0;                   % Velocity [m/s]
% Initialize the data
data.k = 0.9;               % Attenuation factor
data.rho = 1.225;           % Density of medium [kg/m^3]
data.mb = 1;                % Ball mass [kg]
data.A_b = 0.7;             % Ball Area [m^2]
data.V_b = 0.014;           % Ball Volume [m^3]
data.Cd = 1.17;             % Ball Drag Coefficient
% Time interval for the calculation
time_vec = [0 10];
% Set the initial state of the ball
x0 = [x(1) v(1)];
% Set the ode tolerance high to be as close to the ground as possible
options = odeset('RelTol',1e-9, 'AbsTol',1e-9,'Events',@(t,s) EventFct(t,s,1));
tp = zeros(10000,1);
xp = zeros(10000,2);
% Use ode to propagate the motion of the body for the first bounce
[tt,xx] = ode78(@(t,x) fun(t,x,data),time_vec, x0, options);
iter = 1;
tp(1:length(tt)) = tt;
xp(1:length(tt),:) = xx;
while tt(end) < 10
    % Implement the velocity decrease by the impact with the body
    v_after_contact = -data.k*xx(end,2); 
    % Update the ball's state
    x0 = [0 v_after_contact];
    % Reduce the time interval, eliminating the already propagate part (time interval)
    time_vec = [tt(end) 10];
    % Repeat the process for the rest of the time interval
    [tt2,xx2] = ode78(@(t,x) fun(t,x,data),time_vec, x0, options);
    % Check if the where is the maximum of the propagation
    [~, id] = max(abs(xx2));
    % If the maximum of the propagation is below zero, then the ball
    % reacheed equilibrium and it cannot be propagated anymore
    if xx2(id) < 0
        % Set the rest of the interval with equilibrium conditions
        xx2 = zeros(length(tt2),2);
    end
    % Update the entire time interval
    tt = tt2;
    n = length(nonzeros(tp)) + 1;
    tp(n+1:n+length(tt)) = tt;
    % Save the states in a single vector
    xx = xx2;
    xp(n+1:n+length(tt),:) = xx;
    iter = iter + 1;
end
l = length([0;nonzeros(tp)]);
% Plot the results
figure(1)
plot([0;nonzeros(tp)],xp(1:l,:),'LineWidth',1.25)
xlabel('Time[s]')
ylabel('x(t)')
title('\textbf{Falling ball for $\rho$ = 1.225}','Interpreter','latex');
legend('Height[m]','Velocity[m/s]')
hold on;
grid on;
hold off;

%% Exercise 6.2) rho = 15; 
clearvars -except data x0  options;
time_vec = [0 10];
x0 = [10 0];
% Update the value of density
data.rho = 15;
tp = zeros(10000,1);
xp = zeros(10000,2);
% Use ode to propagate the motion of the body for the first bounce
[tt,xx] = ode78(@(t,x) fun(t,x,data),time_vec, x0, options);
tp(1:length(tt)) = tt;
xp(1:length(tt),:) = xx;
iter = 1;
while tt(end) < 10
    % Implement the velocity decrease by the impact with the body
    v_after_contact = - data.k*xx(end,2); 
    % Update the ball's state
    x0 = [0 v_after_contact];
    % Reduce the time interval, eliminating the already propagate part (time interval)
    time_vec = [tt(end) 10];
    % Repeat the process for the rest of the time interval
    [tt2,xx2] = ode78(@(t,x) fun(t,x,data),time_vec, x0, options);
        % Update the entire time interval
    tt = tt2;
    n = length(nonzeros(tp)) + 1;
    tp(n+1:n+length(tt)) = tt;
    % Save the states in a single vector
    xx = xx2;
    xp(n+1:n+length(tt),:) = xx;
    iter = iter + 1;

    iter = iter + 1;
end
l = length([0;nonzeros(tp)]);
% Plot the results
figure(2)
plot([0;nonzeros(tp)],xp(1:l,:),'LineWidth',1.25)
xlabel('Time[s]')
ylabel('x(t)')
title('\textbf{Falling ball for $\rho$ = 15}','Interpreter','latex');
legend('Height[m]','Velocity[m/s]')
grid on;

%% Exercise 6.3) rho = 60; 
clearvars -except data options;
time_vec = [0 10];
x0 = [10 0];
data.rho = 60;
% Use ode to propagate the motion of the body for the first bounce as it
% doesn't reach the ground
[tt,xx] = ode78(@(t,x) fun(t,x,data),time_vec, x0, options);
% Plot the results
figure(3)
plot(tt,xx,'LineWidth',1.25)
xlabel('Time[s]')
ylabel('x(t)')
title('\textbf{Falling ball for $\rho$ = 60}','Interpreter','latex');
legend({'Height[m]','Velocity[m/s]'})
grid on;

%% Exercise 6.4) rho = 60 h = 3 RK4
%% Point 3 RK4
% Begin time counter
tic
% Initialize
clearvars -except tt xx data x0 time_vec options;

% Initialize the data
rho = data.rho;
V_b = data.V_b;
A_b = data.A_b;
Cd = data.Cd;


% RK4 method 
% Define the system as a function handle
dvdt = { @(t,d) -9.81 - 0.5 * rho * Cd * A_b * d(1) * norm(d(1))/mb + rho * V_b * 9.81; @(t,d) d(1)};

% Define time span  
h = 3;                
% Use RK4 function to calculate the results with h = 3
[values] = rk4(time_vec, x0, h, data);
% Save the results of RK4 into the vector of state and time vector
t_values = values(:,1);
x_values = values(:,2:3);

% Plot the results
figure(4)
plot(t_values, x_values, '-o','LineWidth',1.5);
xlabel('Time [s]');
ylabel('x(t)');
title('\textbf{Runge-Kutta 4th Order Solution for $\rho$ = 60}','Interpreter','latex');
legend('Height[m]','Velocity[m/s]')
grid on;

data.rho = 60;
% Use the rk5 function to update the step size and obtain the RK4 solution
% for the new step size
h = 3;
[T, Y , h1,err,er] = rk5( time_vec, x0, h, data);
% Plot the results from the RK5 function
figure(6)
plot(T, Y , '-', 'LineWidth',1.5);
xlabel('Time [s]');
ylabel('x(t)');
title('\textbf{Runge-Kutta 4th Order Solution for $\rho$ = 60}','Interpreter','latex');
legend('Height[m]','Velocity[m/s]')
grid on;
% Calculate he analytical solutionn with ode78
[tt,xx] = ode78(@(t,x) fun(t,x,data),T, x0, options);
axis([0 10 -2 10])

err_rk4 = abs(xx - Y');

figure(7)
sgtitle('\textbf{Runge-Kutta 4th Order Error for $\rho$ = 60}','Interpreter','latex','FontSize',11);
subplot(1,2,1)
semilogy(T, err_rk4(:,1),'LineWidth',1.5);
title('Height Error')
xlabel('Time [s]');
ylabel('$\Delta$ h[m]','Interpreter','latex');
grid on;
subplot(1,2,2)
semilogy(T, err_rk4(:,2),'LineWidth',1.5);
title('Velocity Error')
xlabel('Time [s]');
ylabel('$\Delta$ v[m/s]','Interpreter','latex');
grid on;

%% END
%% Functions

function [dxdt] = fun(~,vec_in,data)
% Function that represents the dynamic of the ball
% Input: vec_in = state vector with: vec_in(1) = position(height) [m]
%                                    vec_in(2) = velocity         [m/s]   
%       data = structure that contains the environment and ball characteristics
%              constants such as density, Drag Coefficient, ball's mass, area, volume 
% Output: dxdt = containing dxdt(1) = velocity     [m/s]
%                           dxdt(2) = acceleration [m/s^2]
    rho = data.rho;                 
    m_b = data.mb;
    A_b = data.A_b;
    V_b = data.V_b;
    Cd = data.Cd;
    v = vec_in(2);
    % Initialize the vector
    dxdt = zeros(2,1);
    dxdt(1) = v;
    dxdt(2) = -9.81 + (-0.5 * rho * Cd * A_b * v * norm(v))/m_b + (rho * V_b * 9.81)/m_b; 

end

function [value,isterminal,direction] = EventFct(~,xx,direction)
% Event Function that stops the integration when the ball reaches the x =
% 0, hits the ground
% Input:
% xx = the ball's state - xx(1) = height (m)
%                       - xx(2) = velocity (m/s)
% Direction: the ball is falling on the ground

     x_floor(2,1) = 0;
     value = dot(x_floor(1)-xx(1), x_floor(2) - xx(2));
     isterminal = 1;

end

function [values] = rk4(tspan, y0, h, data)
    % f      - Function handle f(t, y) representing the system
    % tspan  - Time range [t0, tf] where t0 is the initial time, and tf is the final time
    % y0     - Initial condition (y(t0) = y0)
    % h      - Step size (time increment)
    
    f = @(t,x) fun(t,x, data);

    t0 = tspan(1);  % Initial time
    tf = tspan(2);  % Final time

    n_steps = round((tf - t0) / h);
    
    % Initialize time and solution arrays
    t_values = zeros(1, n_steps + 1);
    y_values = zeros(2, n_steps + 1);
    
    % Initial conditions
    t_values(1) = t0;
    y_values(:,1) = y0;
    
    % Perform RK4 integration
    for i = 1:n_steps
        t = t_values(i);
        y = y_values(:,i);
        
        % Compute RK4 intermediate steps (k1, k2, k3, k4)
        k1 =  f(t, y);
        k2 =  f(t + h/2, y + k1 * h/2);
        k3 =  f(t + h/2, y + k2 * h/2);
        k4 =  f(t + h, y + h * k3);
        
        % Update the next value of y
        y_values(:,i + 1) = y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        
        % Update the next time value
        t_values(i + 1) = t + h;
        
    end

    values = [t_values;y_values]';

end

function [T, Y , h1,err,er] = rk5( tspan, y0, h, data)
    % f      - Function handle f(t, y) representing the system
    % tspan  - Time range [t0, tf] where t0 is the initial time, and tf is the final time
    % y0     - Initial condition (y(t0) = y0)
    % h      - Step size (time increment)

    f = @(t,x) fun(t,x, data);

    t0 = tspan(1);  % Initial time
    tf = tspan(2);  % Final time
    
    tol = 1e-6;          % Desired tolerance
    h_min = 1e-6;        % Minimum step size
    h_max = 1.0;         % Maximum step size
    %Initialize the data
    T(1) = t0;                  % Initial time 
    y_values4 = zeros(2,100);   % RK4 vector
    y_values5 = zeros(2,100);   % RK5 vector
    err_r = 1;                  % Position error
    err_v = 1;                  % Velocity error
    j = 1;
    h1 = zeros(1,tf/1e-4);
    % Continue the process until the error is under the desired tolerance
    while max(err_r) > tol || max(err_v) > tol
      % Calculate the number of steps for the new step size
      n_steps = round((tf - t0) / h);
      % Erase the values of errors
      err_r = zeros(1,n_steps);
      err_v = zeros(1,n_steps);
      % Clear the values of the state
      Y = zeros(2,n_steps);
      T = zeros(1,n_steps);
      % Set the inital conditions
      Y(:,1) = y0';
      T(1) = t0;
      y = y0';
      t = t0;
      for i = 1:n_steps 
        
        % Compute RK5 intermediate steps (k1, k2, k3, k4, k5, k6)_5
        k15 = h * f(t, y);                                                                               % 1st predictor
        k25 = h * f(t + h * 1/4, y + 1/4  * k15);                                                        % 2nd predictor
        k35 = h * f(t + h * 1/4, y + 1/8  * k15 + 1/8 * k25);                                            % 3rd predictor
        k45 = h * f(t + h * 1/2, y - 0    * k15 - 1/2 * k25 + 1 * k35);                                  % 4th predictor
        k55 = h * f(t + h * 3/4, y + 3/16 * k15 - 0   * k25 + 0 * k35 + 9/16 * k45);                     % 5th predictor
        k65 = h * f(t + h * 1/2, y - 3/7  * k15 + 2/7 * k25 + 12/7 * k35 - 12/7 * k45 + 8/7 * k55);      % 6th predictor
        % Update the next value of y
        y_values5(:,i + 1) = y + 1/90 * (7 * k15 + 32 * k35 + 12 * k45 + 32 * k55 + 7 * k65);               % Corrector
        % Compute RK4 intermediate steps (k1, k2, k3, k4)
        k1 =  f(t, y);                                                          % 1st predictor
        k2 =  f(t + h/2, y + k1 * h/2);                                         % 2nd predictor
        k3 =  f(t + h/2, y + k2 * h/2);                                         % 3rd predictor
        k4 =  f(t + h, y + h * k3);                                             % 4th predictor
        
        % Update the next value of y
        y_values4(:,i + 1) = y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        
        % Error 
        err_r(i + 1) = abs(y_values5(1,i + 1) - y_values4(1,i + 1));
        err_v(i + 1) = abs(y_values5(2,i + 1) - y_values4(2,i + 1));
       % Udate the time step
        t = t + h;
        % Update the y _values
        y = y_values4(:,i+1);
        % Save the state in a vector
        Y(:,i + 1) = y';
        % Save the time step
        T(i+1) = t;
      end
        err_ru = abs(y_values5(1, i + 1) - y_values4(1, i + 1));  % Error in position
        err_vu = abs(y_values5(2, i + 1) - y_values4(2, i + 1));  % Error in velocity
        % Compute the overall error (norm of position and velocity errors)
        err = max([err_ru, err_vu]);
    
        % Calculate the scaling factor
        scaling_factor = min(5, max(0.1, (tol / err)^(1/5)));  % RK5 has order 5
    
        % Update the step size
        h_new = h * scaling_factor;
    
        % Ensure h stays within bounds
        h = min(max(h_new, h_min), h_max);
        h1(j) = h;
        j = j + 1;
    end
    er = [err_r;err_v];
    h1 = nonzeros(h1);
end

function [alpha_new, iter] = NS(alpha,f,fprime, tol)
 % This function simulates the Newton Solver
 % Input :
 % alpha: is the initial value [deg]
 % f: is the function handle
 % fprime: is a fct handle containing the derivative of the function f'
 % tol : the desired tolerance for the Newton Solver
 % Output:
 % alpha_new: the update value obtained from Newton Solver [deg]
 % iter: the number of iteration need by NS to reach the desired tolerance
     iter = 1;
     err = 1;
     max_iter = 50;         % set maximum nr of iterations for computation efficiency
     
     % Repeat the process until the tolerance is satisfied or the maximum
     % nr of iteration is reached
     while iter <= max_iter && err > tol
        % Update the alpha value
        alpha_new = alpha - f(alpha)/fprime(alpha);
    
        % Check for convergence
        err = abs(alpha_new - alpha);
        % Save the new alpha 
        alpha = alpha_new;
        iter = iter + 1;
     end   
end
function [value,isterminal,direction] = EventFct4(~,T_in,direction)
% Function meant to stop the ode45 when the T_in reaches the outside
% temperature
% T_in [K] = Temperature of the body
% T_out[K] = Surrounding  temperature, can be considered as an input, in
% this case was defined in the function
% direction - set the ode to stop when the temperature is decreasing to
% outside temperature
     T_out = 277;
     value = T_in - T_out;
     isterminal = 1;
     
end