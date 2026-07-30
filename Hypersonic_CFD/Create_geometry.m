clc;
close all;
clearvars;
%cd ("VKI\NLAB\")
%% Create the circle

p_2pi_1 = linspace(0,      pi/2,   100);
p_2pi_2 = linspace(pi/2,   pi,     100);
p_2pi_3 = linspace(pi,     3*pi/2, 100);
p_2pi_4 = linspace(3*pi/2, 2*pi,   100);

p_2pi = [p_2pi_1, p_2pi_2, p_2pi_3, p_2pi_4];
p_2pi = unique(p_2pi);

x = 0.2 + 0.1 * cos(p_2pi);
y = 0.2 + 0.1 * sin(p_2pi);
circle1 = [x',y',zeros(size(x))'];
x = 0.2 + 0.125 * sin(p_2pi);
y = 0.2 + 0.125 * cos(p_2pi);
circle2 = [x',y',zeros(size(x))'];
%%
fileID = fopen('points.dat','w');
fprintf(fileID, '%i\n', length(p_2pi));
for i = 1:length(p_2pi)
    fprintf(fileID, '%.6e %.6e %.6e \n', circle1(i,:));
end
fprintf(fileID, '%i\n', length(p_2pi));
for i = 1:length(p_2pi)
    fprintf(fileID, '%.6e %.6e %.6e \n', circle2(i,:));
end
%% Connecting the circles
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.2 + 0.1 * cos(p_2pi(1)),0.2 + 0.1 * sin(p_2pi(1)),0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.2 + 0.125 * cos(p_2pi(1)),0.2 + 0.1 * sin(p_2pi(1)),0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.2 + 0.1 * cos(p_2pi_1(end)),0.2 + 0.1 * sin(p_2pi_1(end)),0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.2 + 0.1 * cos(p_2pi_1(end)),0.2 + 0.125 * sin(p_2pi_1(end)),0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.2 + 0.1 * cos(p_2pi_2(end)),0.2 + 0.1 * sin(p_2pi_2(end)),0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.2 + 0.125 * cos(p_2pi_2(end)),0.2 + 0.1 * sin(p_2pi_2(end)),0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.2 + 0.1 * cos(p_2pi_3(end)),0.2 + 0.1 * sin(p_2pi_3(end)),0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.2 + 0.1 * cos(p_2pi_3(end)),0.2 + 0.125 * sin(p_2pi_3(end)),0]);


%% Outer bottom line — split into segments
fprintf(fileID, '%i\n',2);  % corner to first divider
fprintf(fileID, '%.6e %.6e %.6e \n', [0.0, 0.0, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.5, 0.0, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.5, 0.0, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.0, 0.0, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.0, 0.0, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.5, 0.0, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.5, 0.0, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.0, 0.0, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.0, 0.0, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.2, 0.0, 0]);

%% Outer right line — no dividers, but split at corner offsets
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.2, 0.0, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.2, 0.025, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.2, 0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.2, 0.41-0.025, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.2, 0.41-0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.2, 0.41, 0]);

%% Outer top line — split into segments (right to left)
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.2, 0.41, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.0, 0.41, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.0, 0.41, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.5, 0.41, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.5, 0.41, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.0, 0.41, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.0, 0.41, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.5, 0.41, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.5, 0.41, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.0, 0.41, 0]);

%% Outer left line — split at corner offsets
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.0, 0.41, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.0, 0.41-0.025, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.0, 0.41-0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.0, 0.025, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.0, 0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.0, 0.0, 0]);

%% Inner bottom line — split into segments
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.025, 0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.5, 0.025, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.5, 0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.0, 0.025, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.0, 0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.5, 0.025, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.5, 0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.0, 0.025, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.0, 0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.2-0.025, 0.025, 0]);

%% Inner right line
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.2-0.025, 0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.2-0.025, 0.41-0.025, 0]);

%% Inner top line — split into segments
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.2-0.025, 0.41-0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.0, 0.41-0.025, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.0, 0.41-0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.5, 0.41-0.025, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.5, 0.41-0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.0, 0.41-0.025, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.0, 0.41-0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.5, 0.41-0.025, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.5, 0.41-0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.025, 0.41-0.025, 0]);

%% Inner left line
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.025, 0.41-0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.025, 0.025, 0]);
%% Squares - Bottom Left (correct)
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.025, 0.0, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.025, 0.025, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.0, 0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.025, 0.025, 0]);

%% Squares - Bottom Right (added missing base line)
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.2-0.025, 0.0, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.2-0.025, 0.025, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.2-0.025, 0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.2, 0.025, 0]);

%% Squares - Top Right (added missing top line)
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.2, 0.41-0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.2-0.025, 0.41-0.025, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.2-0.025, 0.41-0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.2-0.025, 0.41, 0]);

fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.025, 0.41-0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.025, 0.41, 0]);
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.025, 0.41-0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.0, 0.41-0.025, 0]);
%%
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.5, 0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.5, 0.41-0.025, 0]);

fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.0, 0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.0, 0.41-0.025, 0]);

fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.5, 0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.5, 0.41-0.025, 0]);

fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.0, 0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.0, 0.41-0.025, 0]);

%% Lower Boundary
%%
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.5, 0.0, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.5, 0.025, 0]);

fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.0, 0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.0, 0.0, 0]);

fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.5, 0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.5, 0.0, 0]);

fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.0, 0.025, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.0, 0.0, 0]);

%% Upper Boundary
%%
fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.5, 0.41, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [0.5, 0.41-0.025, 0]);

fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.0, 0.41, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.0, 0.41-0.025, 0]);

fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.5, 0.41, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [1.5, 0.41-0.025, 0]);

fprintf(fileID, '%i\n',2);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.0, 0.41, 0]);
fprintf(fileID, '%.6e %.6e %.6e \n', [2.0, 0.41-0.025, 0]);




% fprintf(fileID, '%i\n',2);
% fprintf(fileID, '%.6e %.6e %.6e \n', [0.2,0,0]);
% fprintf(fileID, '%.6e %.6e %.6e \n', [0.2,0.1,0]);


fclose(fileID);


%%

%% Read and plot the dat file
fileID = fopen('points.dat', 'r');

figure; hold on; axis equal;
xlabel('x'); ylabel('y');
title('Points from points.dat');
colors = lines(150); % enough colors for all groups
groupIdx = 1;

while ~feof(fileID)
    % Read number of points in this group
    n = fscanf(fileID, '%i', 1);
    if isempty(n), break; end

    % Read the n points (x y z)
    data = fscanf(fileID, '%e %e %e', [3, n])';

    if size(data, 1) == n
        x = data(:, 1);
        y = data(:, 2);
        c = colors(groupIdx, :);

        if n == 1
            plot(x, y, 'o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 5);
        else
            plot(x, y, '-o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 3, 'LineWidth', 1.2);
        end
        groupIdx = groupIdx + 1;
    end
end

fclose(fileID);
grid on;