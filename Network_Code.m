% This code combines both the contacts and degree + Clustering
clc; clearvars;
format shortG;
close all;

%% Load DATA
tic
addpath('C:\Users\mihne\Documents\GitHub\Chrono_Projects\files\spherical_large_no_core_fixed')
format_ast = {"Kleopatra.obj\", "Geographos.obj\","Bennu_v20_200k.obj\"};%"Kleopatra.obj\", "Geographos.obj\","Bennu_v20_200k.obj\"
core = {"spherical_large_no_core_fixed\"}; 
densitiesall  = {"rho 1200\","rho 2000\","rho 2400\","rho 1200\","rho 2000\"}; 
rotationsall =  {"Rotation Period 0.720\", "Rotation Period 0.720\","Rotation Period 1.980\","Rotation Period 3.240\", "Rotation Period 1.260\"};

% Possible vector solutions
% densitiesall = {"rho 3000\","rho 1600\","rho 2800\","rho 2400\","rho 2800\","rho 2000\"};
% rotationsall = {"Rotation Period 0.720\", "Rotation Period 1.260\", "Rotation Period 3.240\","Rotation Period 3.060\", "Rotation Period 1.980\", "Rotation Period 2.520\"};

 linestyle = {"-", ":", "--", "-.", "-"};

 %% Initialize DATA

for doc_i = 1:length(format_ast)
    for core_i = 1:length(core)
        % Put r_i = 1 if the data must go in parallel in r_i and d_i
        for r_i = 1 %:length(rotationssall)
            for d_i = 1:length(densitiesall)


                density_folder = densitiesall{d_i};
                rotation_folder = rotationsall{d_i};
                r_i = d_i;

                % Load a folder to save data
                photo_folder = append('C:\Users\mihne\Documents\GitHub\Chrono_Projects\files\spherical_large_no_core_fixed\Photos_Matlab\',format_ast{doc_i},density_folder,rotation_folder);
                mkdir(photo_folder);

                load_directory =append("C:\Users\mihne\Documents\GitHub\Chrono_Projects\files\",core{core_i},"Chrono_Reintroduction\CutOut_simulation_build\Release\Results_New_Background\",format_ast{doc_i},rotation_folder,density_folder);


                addpath(load_directory);

                % Load the data saved from Chrono in txt file
                timeAll.positions = readmatrix("Positions.txt");
                timeAll.velocities = readmatrix("Velocities.txt");
                timeAll.mass = readmatrix("Mass.txt");
                timeAll.radius = readmatrix("Radius.txt");
                timeAll.angmomentum = readmatrix("AngularMomentum.txt");   % This is the angular velocity
                timeAll.inertia = readmatrix("Inertia.txt");               % Moment of inertia of each body around itself in the form [x,y,z] - > diag([x,y,z])
                timeAll.accforces = readmatrix("AccumulatedForce.txt");    
                timeAll.contforces = readmatrix("ContactForces.txt");      % The sum of all contact forces on each body 
                timeAll.contactpairs = readlines("ContactPairs.txt");
                timeAll.eachforce = readlines("EachForce.txt");            % The contact forces for each contact pair
                timeAll.eachforcesave = timeAll.eachforce;
                timeAll.contactpairssave = timeAll.contactpairs;
                sim_inputs_path = append("C:\Users\mihne\Documents\GitHub\Chrono_Projects\files\",core{core_i},"results\", format_ast{doc_i},density_folder,rotation_folder,"simInputs.txt");
                field = "Universal gravity constant G = ";
                if (format_ast{doc_i} == "Geographos.obj\")
                    sim_inputs_path  = append("C:\Users\mihne\Documents\GitHub\Chrono_Projects\files\spherical_large_no_core_fixed\results\Geographos Radar-based, mid-res.obj\",density_folder,rotation_folder,"simInputs.txt");
                end

                %% 
                Grav = readdata(sim_inputs_path, field);
                % Create colors for different plots
                hex_colors = { "#FF0000";  "#00FF00";   "#0000FF";   "#00FFFF";  "#FF00FF";   "#FFFF00";  "#000000";   "#FFFFFF"; "#0072BD";   "#D95319";  "#EDB120";   "#7E2F8E"; "#77AC30";   "#4DBEEE";  "#A2142F"; "#1A2B3C"; "#4F5E6D"; "#A1B2C3"; "#D4E5F6";"#FF5733"; "#33FF57";"#5733FF"; "#C0C0C0";"#800000"; "#008000";"#000080"; "#FFA500";"#4B0082"; "#EE82EE";"#4682B4"; "#20B2AA";"#DC143C"; "#8B0000";"#556B2F2"; "#2F4F4F"};
                % Colors = { Red, Green, Blue, Cyan, Magenta, Yellow, Black, White, Dark, Blue, Dark Orange, Dark Yellow, Dark Purple, Medium Green, Light Blue, Dark Red};
                %%

                % Check to see other distance multipliers
                threshold.distance = 10 * (timeAll.radius(1) + timeAll.radius(2));
                sizes.time = length(timeAll.mass(:,1));
                sizes.bodies = length(timeAll.mass(1,:));
                A.distance = zeros(sizes.bodies,sizes.bodies);
                %Initialize data
                centralities_distance.closeness = zeros(sizes.time, sizes.bodies);
                centralities_distance.closesum = zeros(sizes.time, sizes.bodies);

                % Check the number of instances that can be created

                % Set to maximum 96 only for Kleopatra and Geographos, Bennu becomes too large (~ 60RAM)
                pairs.time_i = zeros(sizes.bodies,sizes.bodies, 96);  % [Matrix; time instance]
                pairs.forces_time_i = zeros(sizes.bodies,sizes.bodies);
                pairs.forces_time_All = zeros(sizes.bodies, sizes.bodies,96);  % [Matrix; time instance]
                pairs.time_All = zeros(sizes.bodies,sizes.bodies, 96);         % [Matrix; time instance]
                % pairs.time_All2 = zeros(sizes.bodies, sizes.bodies,15);
                sumtimeAll = zeros(sizes.bodies,sizes.time);
                time_vec = round(linspace(1,sizes.time,96));                   % [Matrix; time instance]
                time_ii = 1;

                % Create a system tolerance
                tol = 1e-6;



                % Transform Adjacency List into Adjacency Matrix
                for time_i = 1:size(timeAll.positions,1)

                    % Save the lines that have more values in common in pairs
                    savelines = [];

                    if(any(time_vec == time_i))
                    length_max_local = 0;
                    for pairs_i = 1:sizes.bodies
                        length_max_local =  max(length_max_local, length(str2num(timeAll.contactpairs(pairs_i))));
                    end
                    
                    length_max = length_max_local;

                    if (length_max > 2)
                        for i = 1:sizes.bodies
                            line = str2num(timeAll.contactpairs(i)) + 1;
                            timell(time_ii) = line(1);
                            line(1:2) = [];
                            line2 = str2num(timeAll.eachforce(i));
                            id_line_init = line2(2) + 1;
                            line2(1:2) = [];        
                            if(length(line) >= 1)
                                pairs.time_All(i,line(line > i),time_ii) = 1;
                                pairs.time_All2(i,line,time_ii) = 1;
                                % Put in also the contact forces
                                for line_i = 1:length(line)
                                    line3 = str2num(timeAll.eachforce(line(line_i)));
                                    line4 = str2num(timeAll.contactpairs(line(line_i))) + 1;
                                    id_line_comp = line4(2);
                                    line3(1:2) = [];
                                    line4(1:2) = [];

                                    count = 0;
                                    for line_j = 1:length(line2)
                                        if (any ( abs(line2(line_j) - line3) < tol))     
                                            count = count + 1;
                                            save_j = line_j;
                                        end
                                    end
                                    if (count > 1)
                                        % If more than to values are identical save the pairs in a vector and move on
                                        savelines = [savelines; [id_line_init id_line_comp time_ii]];
                                        continue;
                                    else
                                        pairs.forces_time_All(i,line(line_i),time_ii) = line2(save_j);
                                    end

                                end
                            end
                        end
                    end

                % check the pairs between the combinations of pairs that have not been
                % imposed yet
                for k = 1:size(savelines,1)

                    i = savelines(k,1);
                    j = savelines(k,2);
                    t = savelines(k,3);

                    % Skip if already assigned
                    if pairs.forces_time_All(i,j,t) ~= 0
                        continue;
                    end

                    % Companion already knows → inherit
                    if pairs.forces_time_All(j,i,t) ~= 0
                        pairs.forces_time_All(i,j,t) = pairs.forces_time_All(j,i,t);
                        continue;
                    end

                    % Extract force lists
                    f_i = str2num(timeAll.eachforce(i)); f_i(1:2) = [];
                    f_j = str2num(timeAll.eachforce(j)); f_j(1:2) = [];

                    used_i = nonzeros(pairs.forces_time_All(i,:,t));
                    used_j = nonzeros(pairs.forces_time_All(j,:,t));

                    rem_i = setdiff(f_i, used_i);
                    rem_j = setdiff(f_j, used_j);

                    % Find common force (tolerance-based)
                    common = rem_i(ismembertol(rem_i, rem_j, tol));

                    if isempty(common)
                        continue;   % NOTHING written → safe
                    end

                    % If multiple identical forces, sum or take one (your choice)
                    force_value = sum(common);

                    pairs.forces_time_All(i,j,t) = force_value;
                    pairs.forces_time_All(j,i,t) = force_value;
                end


                    sumtimeAll(1:sizes.bodies,time_ii) = sum(pairs.time_All(:,:,time_ii),2);
                    time_ii = time_ii + 1;
                    pairs.time_i = timeAll.contactpairs(1:sizes.bodies);
                    pairs.forces_time_i = timeAll;
                    end
                timeAll.contactpairs(1:sizes.bodies) = [];
                timeAll.eachforce(1:sizes.bodies) = [];
                end
                % Check symmetry - Adjacency must be symmetric
                
                for i = 1:length(time_vec)
                    pairs.time_All(:,:,i) = pairs.time_All(:,:,i) + transpose(pairs.time_All(:,:,i));
                    pairs.time_All(:,:,i) = abs(pairs.time_All(:,:,i));
                    if (issymmetric(pairs.forces_time_All(:,:,i)) && issymmetric(pairs.time_All(:,:,i)))
                        disp(true);
                    else
                        disp(i);
                        disp(false);
                    end
                end
                toc


                %% Create an initialization for smoother running
                E_time = zeros(sizes.bodies,numel(time_vec));
                Ltotal = zeros(1,numel(time_vec));
                L_total_vec = zeros(sizes.bodies,3);
                L_total = zeros(numel(time_vec),3);
                Adj = cell(1,numel(time_vec));

                
                for time_i = 1:length(time_vec)
                    % Clear some values
                clear E.c E.pe E.r E.t E.escape 
                tic
                    % Initialize each time step
                    time_ii = time_vec(time_i);
                    disp(time_ii)
                    positions = transpose(reshape(timeAll.positions(time_ii,:),[3,sizes.bodies]));  
                    velocities = transpose(reshape(timeAll.velocities(time_ii,:),[3,sizes.bodies]));
                    radius = transpose(timeAll.radius(time_ii,:));
                    Mass.all = transpose(timeAll.mass(time_ii,:));
                    acc_forces = transpose(reshape(timeAll.accforces(time_ii,:),[3,sizes.bodies]));
                    contact_forces = transpose(reshape(timeAll.contforces(time_ii,:),[3,sizes.bodies]));
                    contact_forces_norm = vecnorm(contact_forces,2,2);
                    acc_forces_norm = vecnorm(acc_forces,2,2);
                    omega = transpose(reshape(timeAll.angmomentum(time_ii,:),[3,sizes.bodies]));
                    inertia = transpose(reshape(timeAll.inertia(time_ii,:),[3,sizes.bodies]));

                %% Calculate the limit for simulation time-step
                R = zeros(sizes.bodies,sizes.bodies);
                V = zeros(sizes.bodies,sizes.bodies);
                for i = 1:sizes.bodies
                    for j = i+1:sizes.bodies
                        R(i,j) = norm(positions(i,:) - positions(j,:));
                        R(j,i) = R(i,j);
                        V(i,j) = norm(velocities(i,:) - velocities(j,:));
                        V(j,i) = V(i,j);
                    end
                end

                %%
                % Create the adjacency matrices : contacts, dists, forces, vels, omega,
                % inertia

                    % Contact pairs

                    A.pairs = pairs.time_All(:,:,time_i);
                    G.pairs = graph(A.pairs, string(1:sizes.bodies));
                    density.pairs(time_i) = 2 * size(G.pairs.Edges,1)/(size(G.pairs.Nodes,1) *(size(G.pairs.Nodes,1) - 1));

                    % Distances - considered the one between each pair of bodies
                    dist = zeros(sizes.bodies,sizes.bodies);
                    for bodies_i = 1:sizes.bodies
                        % + 1 to avoid self loops
                        for bodies_j = bodies_i + 1:sizes.bodies   
                            dist(bodies_i,bodies_j) = norm(positions(bodies_i,:) - positions(bodies_j,:));
                                A.distance(bodies_i, bodies_j) = dist(bodies_i,bodies_j);               % 1
                                A.distance(bodies_j, bodies_i) = dist(bodies_i,bodies_j);               % 1
                        end
                    end

                    G.distances = graph(A.distance, string(1:sizes.bodies));

                    density.distances(time_i) = 2 * size(G.distances.Edges,1)/(size(G.distances.Nodes,1) *(size(G.distances.Nodes,1) - 1));

                    % Forces 
                    % Adjacency matrix that takes the pairs of forces
                    A.forces = abs(pairs.forces_time_All(:,:,time_i));
                    G.forces = graph(A.forces, string(1:sizes.bodies));


                    density.forces(time_i) = 2 * size(G.forces.Edges,1)/(size(G.forces.Nodes,1) *(size(G.forces.Nodes,1) - 1));

                    % Gravity Adjacency
                    Contact = zeros(sizes.bodies,sizes.bodies);
                    for i = 1:sizes.bodies
                        for j = i+1:sizes.bodies
                            dist_ij = norm(positions(i,:) - positions(j,:));
                            if dist_ij <= radius(i) + radius(j)
                                Contact(i,j) = 1;
                                Contact(j,i) = 1;
                            end
                            A.gravity(i,j) = Grav * Mass.all(i) * Mass.all(j)/dist_ij^2;
                            A.gravity(j,i) = Grav * Mass.all(i) * Mass.all(j)/dist_ij^2;
                        end
                    end
                    G.gravity = graph(A.gravity, string(1:sizes.bodies));

                    % This density approaches 1 which is normal as each vertex is in
                    % contact with all the others
                    density.gravity(time_i) = 2 * size(G.gravity.Edges,1)/(size(G.gravity.Nodes,1) *(size(G.gravity.Nodes,1) - 1));
                    toc
                    %% Calculate the adjacency for omega
                    for i = 1:sizes.bodies
                        for j = i+1:sizes.bodies
                            A.omega(i,j) = norm(omega(i,:) - omega(j,:));
                            A.omega(j,i) = A.omega(i,j);
                        end
                    end
                    G.omega = graph(A.omega, string(1:sizes.bodies));

                    %% Calculate the adjacency for angular momentum
                   toc
                    Mass.total = sum(Mass.all);
                    r_com = sum(Mass.all .* positions) / Mass.total ;
                    velocity_com = sum(Mass.all .* velocities) /Mass.total;

                    L = zeros(sizes.bodies,3);

                    % Compute the angular momentum of each particle
                    for i = 1:sizes.bodies
                        L(i,:) = cross((positions(i,:) - r_com),Mass.all(i)*(velocities(i,:) - velocity_com));
                    end
                    % So the Angular momentum and angular velocity can be calculated as the
                    % correlation dot(Li,Lj)/|Li|*|Lj|
                    L_norm = transpose(vecnorm(L'));
                    for i = 1: sizes.bodies
                        for j = i+1:sizes.bodies
                            A.ang_momentum(i,j) = norm(L(i,:) - L(j,:));
                            A.ang_momentum(j,i) = norm(L(i,:) - L(j,:));
                        end
                    end
                    G.ang_momentum = graph(A.ang_momentum,string(1:sizes.bodies));

                    %% 
                    % A.acc_forces = zeros(sizes.bodies,sizes.bodies);
                    % % Accumulated forces
                    for i = 1: sizes.bodies
                        for j = 1:sizes.bodies
                            if i == j
                                A.acc_forces(i,i) =  0;
                            else
                                A.acc_forces(i,j) = Grav * Mass.all(i) * Mass.all(j)/norm(positions(i,:) - positions(j,:));
                                A.acc_forces(j,i) = A.acc_forces(i,j);
                            end
                        end
                    end
                    % % % Keep it like this Acc - Forces
                    if (issymmetric(A.acc_forces))
                        A.forces_total = abs(A.acc_forces - A.forces);
                    end
                    G.forces_total = graph(A.forces_total, string(1:sizes.bodies));

                    % Calculate the difference of forces and the energies
                    E.c = zeros(sizes.bodies,1);
                    E.pe = zeros(sizes.bodies,1);
                    E.r = zeros(sizes.bodies,1);
                    
                    % Formulas for the energies
                    % E.c = 1/2mv2
                    % E.p_self = -3/5*G*m^2/R - R = radius
                    % E.p_unit = -G*mi*mj/rij - rij = distance btw the particles 
                    % E.r = 1/2 omega' I omega

                    for i = 1:sizes.bodies
                        E.c(i) = Mass.all(i) * norm(velocities(i,:))^2/2;
                        for j = 1:sizes.bodies
                            if i == j
                                continue;
                            end
                            E.pe(i) = E.pe(i) - Grav * Mass.all(j) * Mass.all(i) / norm(positions(i,:) - positions(j,:));
                        end
                        % Self Potential energy
                        U_ii = -3/5*Grav * Mass.all(i)/radius(i); 
                        % Rotation energy
                        E.r(i) = 0.5 * omega(i,:) * diag(inertia(i,:)) * omega(i,:)';
                        E.t(i) =  E.c(i) + E.r(i) + U_ii + E.pe(i); %% + E.pe + E.pc which one??????????
                        % Escape Energy
                        E.escape(i) = Grav * (Mass.total)/norm(positions(i,:)); 
                    
                        % Save E.r + E.c
                        E.r_k (i,time_i) = E.r(i) + E.c(i) + U_ii;
                    end
                    delta_E = E.t(:) - E.escape(:);
                    E.total(time_i) = sum(E.t);
                    E_time(:,time_i) = E.t';


                    E.mags = abs(E.t(1:end));   % Energy magnitude (absolute value)

                    E.normalized = (E.mags - min(E.t(1:end))) / (max(E.t(1:end)) - min(E.t(1:end)));  % Scale to [0,1]
                    % Normalize marker size

                    toc
                    %% Calculate the relative energies and the adjacency matrix
                    
                    A.energies = zeros(sizes.bodies,sizes.bodies);
                    D = pdist(E.t');          % E is a column vector of energies
                    sigma = median(pdist(E.t')); 
                    A.energies = exp(-(E.t' - E.t).^2/(2 * sigma^2));
                    for i = 1:sizes.bodies , A.energies(i,i) = 0; end
                    A.energies = pdist2(E.t', E.t');

                    G.energies = graph(abs(A.energies), string(1:sizes.bodies));

                    A.potentials = zeros(sizes.bodies, sizes.bodies);

                    A.potentials = abs(E.pe' - E.pe);

                    G.potentials = graph(A.potentials, string(1:sizes.bodies));

                    A.kinetics = abs(E.c' - E.c);
                    G.kinetics = graph(A.kinetics, string(1:sizes.bodies));

                    

                    %% This section is used to create the clusters
                    %     clustervalue = clusterpairs(A.pairs,zeros(10));
                    % 
                    %     figure()
                    %     hold on;
                    %     grid on;
                    %     num_clusters = length(clustervalue);
                    %     colors = hsv(num_clusters);
                    %     for i = 1:length(clustervalue)
                    %          if length(clustervalue{i}) > 1
                    %             scatter3(positions(clustervalue{i},1),positions(clustervalue{i},2),positions(clustervalue{i},3),'filled','MarkerFaceColor',colors(i,:))
                    %             % Plot the edges
                    %             for j = 1:length(clustervalue{i})
                    %                 pairs_j = find(A.pairs(clustervalue{i}(j),:));
                    %                 for pairs_i = 1:length(pairs_j)
                    %                     x = [positions(clustervalue{i}(j),1), positions(pairs_j(pairs_i),1)];
                    %                     y = [positions(clustervalue{i}(j),2), positions(pairs_j(pairs_i),2)];
                    %                     z = [positions(clustervalue{i}(j),3), positions(pairs_j(pairs_i),3)];
                    % 
                    %                     plot3(x, y, z, 'Color',[0 0 0]+0.1);  % 'k' = black line
                    %                 end
                    %             end
                    %         else
                    %             body = clustervalue{i};
                    %             scatter3(positions(body,1),positions(body,2),positions(body,3),'filled','MarkerFaceColor',hex_colors{end})
                    %         end
                    %     end
                    %     xlabel('x [m]','Interpreter','latex','FontSize',16); ylabel('y [m]','Interpreter','latex','FontSize',16); zlabel('z [m]','Interpreter','latex','FontSize',16);
                    %     view([-15 15 5]);
                    %     title(sprintf('T = %.2f hrs', round(time_vec(time_i) * 300/ 3600,2)));  
                    %      saveas(gcf,append("C:/Users/mihne/Documents/GitHub/Chrono_Projects/files/spherical_large_no_core_fixed/Phtots_thesis/Energies/Bennu/",...
                    %      append("Cluster Mine"," ",string(time_i)," "  ,".svg")));
                    % hold off;
                    % 
                    % 
                    % %    Calculate the Energies with respect to the clusters
                    %     % Each value must be recalculated wrt to the cluster - R = position
                    %     % relative to cluster's CoM, V = velocity relative to cluster's
                    %     % velocity, rotation the same and moment of inertia
                    % 
                    %      sizes.cluster = numel(clustervalue);
                    %     CoMc = zeros(sizes.cluster,3);
                    %     min_dist_cluster = ones(sizes.cluster,1);
                    % 
                    %     E.t_each_cluster = {};
                    % 
                    %     cluster_velocities = zeros(sizes.cluster,3);
                    %     % Calculate the CoM of each cluster
                    %     for i = 1:sizes.cluster
                    %         Mass.cluster(i) = sum(Mass.all(clustervalue{i}));
                    %         CoMc(i,:) = sum(positions(clustervalue{i},:).* Mass.all(clustervalue{i}),1)/Mass.cluster(i);
                    %         [~,id] = min(vecnorm(ones(numel(clustervalue{i}),3) .* CoMc(i,:) - positions(clustervalue{i},:),2,2)); 
                    %         min_dist.cluster(i) =  clustervalue{i}(id); % vector that stores the points closest to the CoM
                    %         min_dist.all{i} = ones(numel(clustervalue{i}),3) .* CoMc(i,:) - positions(clustervalue{i},:);
                    %         cluster_velocities(i,:) = sum(velocities(clustervalue{i},:).* Mass.all(clustervalue{i}),1)/Mass.cluster(i);
                    %         relative_positions{i} = positions(clustervalue{i},:) - CoMc(i,:);
                    %         relative_velocities{i} = velocities(clustervalue{i},:) - cluster_velocities(i,:); 
                    % 
                    %     end
                    % 
                    %      cluster_velo = velocities(min_dist_cluster,:);
                    % 
                    %      E.escape_each_cluster = cell(1,sizes.cluster);
                    %      E.c_each_cluster = cell(1,sizes.cluster);
                    %      E.pe_each_cluster = cell(1,sizes.cluster);
                    % 
                    %      figure(5999+time_ii)
                    %      colormap(cool)
                    %      grid on;
                    % 
                    %      for i = 1:sizes.cluster     
                    %         % Initialize
                    %         E.c_each_cluster{i} = zeros(numel(clustervalue{i}),1);
                    %         E.pe_each_cluster{i} = zeros(numel(clustervalue{i}),1);
                    %         E.r_each_cluster{i} = zeros(numel(clustervalue{i}),1);
                    %         % Calculate the escape energy        
                    %         for bod_i = 1:numel(clustervalue{i})
                    % 
                    %             if (numel(clustervalue{i}) ~= 1)
                    %                 E.escape_each_cluster{i}(bod_i) = Grav * (Mass.cluster(i))/norm(relative_positions{i}(bod_i));
                    %                 % Calculate the entire energy of each body
                    %                 E.c_each_cluster{i}(bod_i)= 0.5 * Mass.all(clustervalue{i}(bod_i)) * norm(relative_velocities{i}(bod_i,:))^2;
                    % 
                    %                 for bod_j = 1:numel(clustervalue{i})
                    %                     if bod_i == bod_j
                    %                         continue;
                    %                     end
                    %                     % Potential btw two sphere
                    %                     E.pe_each_cluster{i}(bod_i) = E.pe_each_cluster{i}(bod_i) - Grav * Mass.all(clustervalue{i}(bod_j)) * Mass.all(clustervalue{i}(bod_i)) / norm(relative_positions{i}(bod_i,:) - relative_positions{i}(bod_j,:));
                    %                 end
                    %                 % Self Potential energy
                    %                 U_ii = -3/5*Grav * Mass.all(clustervalue{i}(bod_i))/radius(clustervalue{i}(bod_i)); 
                    %                 % Rotation energy
                    %                 E.r_each_cluster{i}(bod_i) = 0.5 * omega(clustervalue{i}(bod_i),:) * diag(inertia(clustervalue{i}(bod_i),:)) * omega(clustervalue{i}(bod_i),:)';
                    %                 E.t_each_cluster{i}(bod_i) =  E.c_each_cluster{i}(bod_i) + E.r_each_cluster{i}(bod_i) + U_ii + E.pe_each_cluster{i}(bod_i);
                    %             else
                    %                 E.t_each_cluster{i}(bod_i) = 0;  
                    %             end
                    %         end
                    % 
                    %         colormap(cool)
                    %         pointsizes = (E.t_each_cluster{i} - min(E.t_each_cluster{i}))/(max(E.t_each_cluster{i} - min(E.t_each_cluster{i})));
                    %         scatter3(positions(clustervalue{i},1),positions(clustervalue{i},2),positions(clustervalue{i},3),36,E.t_each_cluster{i}','filled')
                    %         c = colorbar;
                    %         c.FontSize = 14;
                    %         title(sprintf('T = %.2f hrs', round(time_vec(time_i) * 300/ 3600,2)));
                    %         xlabel('x [m]','Interpreter','latex','FontSize',16)
                    %         ylabel('y [m]','Interpreter','latex','FontSize',16)
                    %         zlabel('z [m]','Interpreter','latex','FontSize',16)
                    %         hold on;
                    %      end
                    %      saveas(gcf,append("C:/Users/mihne/Documents/GitHub/Chrono_Projects/files/spherical_large_no_core_fixed/Phtots_thesis/Energies/Bennu/",...
                    %      append("Cluster Energies"," ",string(time_i)," "  ,".svg")));
                    %      %title("Cluster Energies")
                    %     hold off;

                     %% Hierarchical Clustering

                     rng("default")
              %      YY = squareform(pdist(positions,"euclidean"));
              %      Y = YY;
              %      Z = linkage(Y,'single','euclidean');
              %      c = cophenet(Z,Y);
              %      incos = inconsistent(Z);
              %      T1 = cluster(Z,'maxclust',numel(clustervalue));
              %      figure() 
              %      colormap("lines");
              %      scatter3(positions(:,1),positions(:,2),positions(:,3),30,T1,'filled');
              %      set(gca,"FontSize",14)
              %      hold on;
              %      grid on;
              %      title(sprintf('T = %.2f hrs', round(time_vec(time_i) * 300/ 3600,2)));
              %      xlabel('x [m]','Interpreter','latex','FontSize',16)
              %      ylabel('y [m]','Interpreter','latex','FontSize',16)
              %      zlabel('z [m]','Interpreter','latex','FontSize',16)
              %      view([90 90 90]);
              %      grid on; hold off;
              %      saveas(gcf,append("C:\Users\mihne\Documents\GitHub\Chrono_Projects\files\spherical_large_no_core_fixed\Phtots_thesis\Clusters\",append('Spectral_clustering_position',string(time_vec(time_i))),'.png'))

              %      clus = cell(1,max(T1));
              %      for i = 1:length(T1)
              %          val_clus = T1(i);
              %          clus{val_clus} = [clus{val_clus}, i];
              %      end

                    %% Energy
              %      rng("default")
              %      YY = squareform(pdist(E.t',"euclidean"));
              %      Y = YY;
              %      Z = linkage(Y,'single','euclidean');
              %      c = cophenet(Z,Y);
              %      incos = inconsistent(Z);
              %      T1 = cluster(Z,'maxclust',numel(clustervalue) + 10);
              %      figure() 
              %      colormap("lines");

              %     scatter3(positions(:,1),positions(:,2),positions(:,3),30,T1,'filled');
              %      set(gca,"FontSize",14)
              %      hold on;
              %      grid on;
              %      title(sprintf('T = %.2f hrs', round(time_vec(time_i) * 300/ 3600,2)));
              %      xlabel('x [m]','Interpreter','latex','FontSize',16)
              %      ylabel('y [m]','Interpreter','latex','FontSize',16)
              %      zlabel('z [m]','Interpreter','latex','FontSize',16)
              %      view([90 90 90]);     grid on; hold off;
              %      clus = cell(1,max(T1));
              %      for i = 1:length(T1)
              %          val_clus = T1(i);
              %          clus{val_clus} = [clus{val_clus}, i];
              %      end



                    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    % Check sigma and epsillon
                %    sigma = 0.4;
                %    epsilon = 0.75;
                %    S.energies = exp( - (A.energies.^2) / (2 * sigma^2) );
                %    S.energies = (S.energies > epsilon) .* S.energies;
                %    for i = 1:sizes.bodies
                %        S.energies(i,i) = 0.1;
                %    end
                % 
                %    sigma = 1;
                %    epsilon = 0.75;
                %    S.ang_momentum = exp( - ((A.ang_momentum/1e4).^2) / (2 * sigma^2) );
                %    S.ang_momentum = (S.ang_momentum > epsilon) .* S.ang_momentum;
                %    for i = 1:sizes.bodies
                %        S.ang_momentum(i,i) = 1e-4;
                %    end
                % 
                    % Add the self loops
                %    [clusters_markov.positions,points2clusters.positions] = markov_clustering(A.distance +  0.25 * diag(max(A.distance)),100, 2, 10, 1e-9);
                    % % For A.pairs the best solution of exp power coefficient is between 1.1 and
                    % % 1.3
                %    [clusters_markov.pairs,points2clusters.pairs] = markov_clustering(A.pairs + 0.25 * eye(size(A.pairs)),50, 2, 1.1, 1e-6);
                %    [clusters_markov.forces_total,points2clusters.forces_total] = markov_clustering(A.forces_total + 0.25 * diag(max(A.forces_total)),100, 2, 1.25, 1e-6);
                %    [clusters_markov.energies,points2clusters.energies] = markov_clustering(A.energies + 0.1 * diag(max(A.energies)),100, 2, 3, 1e-3);
                %    [clusters_markov.ang_momentum,points2clusters.ang_momentum] = markov_clustering(A.ang_momentum + 0.1 * diag(max(A.ang_momentum)),100, 2, 4, 1e-3);

                    % Plot markov clusters
                %   plotmarkov(clusters_markov.positions, ' Positions',positions,points2clusters.positions,'Positions Markov', time_vec,time_i,photo_folder);
                %   plotmarkov(clusters_markov.pairs, ' Pairs',positions,points2clusters.pairs,'Pairs Markov',time_vec,time_i,photo_folder);
                %   plotmarkov(clusters_markov.energies, ' Energies',positions,points2clusters.energies,'Energies Markov',time_vec, time_i,photo_folder);
                %   plotmarkov(clusters_markov.ang_momentum, ' Ang Momentum',positions,points2clusters.ang_momentum,'Angular momentum Markov',time_vec, time_i,photo_folder);
                %   plotmarkov(clusters_markov.forces_total,' Forces Total',positions,points2clusters.forces_total,'Forces Total Markov', time_vec, time_i,photo_folder);

                % Markov clustering is good only for pairs and distance not for the others, try dbscan, hdbscan, hierarchical, and others
                % Save the cluster value of each particle in time

                %   clusterassign.pairs(:,time_i) = points2clusters.pairs;
                %   clusterassign.positions(:,time_i) = points2clusters.positions;
                %c  lusterassign.energies(:,time_i) = points2clusters.energies;
                %   clusterassign.ang_momentum(:,time_i) = points2clusters.ang_momentum;
                %   clusterassign.forces_total(:,time_i) = points2clusters.forces_total;


                    %% Graph properties

                %% Degree Distribution 

                 degree.distance_common = sum(A.distance,2);
                 [unique_degrees, ~, idx] = unique(degree.distance_common);
                 degree_counts = accumarray(idx, 1);
                 degree_distribution.distance_common = degree_counts;
                 mean_degree = sum(unique_degrees .* degree_distribution.distance_common);
                 mean_deviation = sum(abs(unique_degrees - mean_degree) .* degree_distribution.distance_common);


                clear degree_counts uniques_degrees idx;

                
                tic
                degree.pairs_common = sum(A.pairs,2);
             
                figure(9090 + time_i)
                subplot(1, numel(time_vec),time_i)
                history =  histogram(degree.pairs_common,'FaceColor',cool(1));
                 
                 % Get bin center positions and heights
                binCenters{time_i} = history.BinEdges(1:end-1) + history.BinWidthbinCenters/2;
                binHeights{time_i} = history.Values;
                 
                xlabel('Node Degree [-]','Interpreter','latex','FontSize',14)
                ylabel('Count [-]','Interpreter','latex','FontSize',14)
                if (time_ii == time_vec(end))
                    saveas(gcf,append(photo_folder,append('Degree Distribution ',string(time_vec(time_i)),'.png')))
                 end
                avg_degree.pairs(time_i,d_i, r_i) = mean(degree.pairs_common);
             
             
                figure(898989)
                hold on;
             for i = 1:length(binCenters)
                if i >= 5
                    markerstyle = 'o';
                else
                    markerstyle = 'none';
                end
                plot(binCenters{i},binHeights{i},'LineWidth',2,'Marker',markerstyle,'DisplayName', append("Time = ", string(round(time_vec(i) * 300 / 3600,1))))
             end
             grid on;
                legend;
                xlabel('Node Degree [-]','Interpreter','latex','FontSize',14)
                ylabel('Count [-]','Interpreter','latex','FontSize',14)
                hold off;
                toc

             avg_degree.pairs_sorted(:,1:2) = avg_degree.pairs(:,2:3);
             avg_degree.pairs_sorted(:,3) = avg_degree.pairs(:,1);
             avg_degree.pairs_sorted(:,4:5) = avg_degree.pairs(:,4:5);
                 figure()
                 plot(time_vec,avg_degree.pairs_sorted,'LineWidth',2)
                 xlabel('Time [hrs]','Interpreter','latex','FontSize',16)
                 ylabel('Average Degree [-]','Interpreter','latex','FontSize',16)
                 grid on;
                rotationsall = {"Rotation Period 2.520\","Rotation Period 0.288\","Rotation Period 1.260\","Rotation Period 1.980\","Rotation Period 3.060\","Rotation Period 3.240\"};
                legend(["Rotation Period 21.8 hrs","Rotation Period 5 hrs","Rotation Period 3.2 hrs","Rotation Period 2.5 hrs","Rotation Period 2.1 hrs","Rotation Period 1.9 hrs"],'FontSize',12)
                 legend(["$\rho = 1200$", "$\rho = 1600$", "$\rho = 2000$","$\rho = 2400$", "$\rho = 2800$", "$\rho = 3000$"],'Interpreter','latex','FontSize',12);
                 xticks(time_vec);
                 xticklabels(string(round(time_vec * 300 / 3600,2)));
            
             for i = 1:size(clustertime,1)
                 for j = 1:size(clustertime,2)
                     c = clustertime{i,j};
                     clusterlengths = cellfun(@length,c);
                     ccc(i,j) = sum(nnz(clusterlengths == 1));
                     ccc4(i,j) = sum(nnz(clusterlengths < 4));
                 end
             end
                figure()
                hold on;
                for i = 1:size(ccc,2)
                    hold on;
                    plot(time_vec, ccc(i),'LineWidth',2,'Color',hex_colors{i})
                end
                plot(time_vec,ccc)
                grid on;
                xlabel('Time [hrs]','Interpreter','latex','FontSize',16)
                xticks(time_vec);
                xticklabels(string(round(time_vec * 300 / 3600,2)));



            %%

                degree.forces_common = sum(A.forces,2);
                degree.forces_total_common = sum(A.forces_total,2);
                degree.energies_common = sum(A.energies, 2);
                degree.potentials_common = sum(A.potentials, 2);
                degree.omega_common = sum(A.omega, 2);
                degree.ang_momentum_common = sum(A.ang_momentum, 2);
                degree.omega_common = sum(A.omega, 2);
                degree.kinetics_common = sum(A.kinetics, 2);


                %% Average degree

            % This calculates the average degree
                degree.pairs = sum(A.pairs,2) / max(sum(A.pairs,2));
                % figure()
                % histogram(degree.pairs)
            % This calculates the average degree
                degree.distance = sum(A.distance,2) / max(sum(A.distance,2));
                degree.forces = sum(A.forces,2) / max(sum(A.forces,2));
                degree.forces_total = sum(A.forces_total,2) / max(sum(A.forces_total,2));
                degree.energies = sum(A.energies,2) / max(sum(abs(A.energies),2));
                degree.potential = sum(A.potentials,2) / max(sum(A.potentials,2));
                degree.omega = sum(A.omega,2) / max(sum(A.omega,2));
                degree.ang_momentum = sum(A.ang_momentum,2)/ max(sum(A.ang_momentum,2));
                degree.kinetics = sum(A.kinetics,2)/ max(sum(A.kinetics,2));

            % Calculate the centralization index 
                maxi(1) = max(degree.distance);
                maxi(2) = max(degree.forces);
                maxi(3) = max(degree.forces_total);
                maxi(4) = max(degree.pairs);
                maxi(5) = max(degree.energies);
                sum_diff(:,1) = sum(maxi(1) - degree.distance);
                sum_diff(:,2) = sum(maxi(2) - degree.forces);
                sum_diff(:,3) = sum(maxi(3) - degree.forces_total);
                sum_diff(:,4) = sum(maxi(4) - degree.pairs);
                sum_diff(:,5) = sum(maxi(5) - degree.energies);
                max_possible_diff = (sizes.bodies - 1) * (sizes.bodies - 2);

                centralization_index.forces(time_i,r_i,d_i) = sum_diff(:,2) / max_possible_diff;
                centralization_index.forces_total(time_i,r_i,d_i) = sum_diff(:,3) / max_possible_diff;
                centralization_index.distances(time_i,r_i,d_i) = sum_diff(:,1) / max_possible_diff;
                centralization_index.pairs(time_i,r_i,d_i) = sum_diff(:,4) / max_possible_diff; 
                centralization_index.energies(time_i,r_i,d_i) = sum_diff(:,5) / max_possible_diff;

                %Forces Common Correlation
                correlation_index.forces(time_i, d_i) = degreeCorrelation(A.forces, degree.forces_common);

                %Forces Total Common Correlation
                correlation_index.forces_total_common = degreeCorrelation(A.forces_total, degree.forces_total_common);

                %Pairs Correlation
                correlation_index.pairs(time_i, d_i) = degreeCorrelation(A.pairs, degree.pairs);

                % Distance Correlation
                correlation_index.distance(time_i, d_i) = degreeCorrelation(A.distance, degree.distance);

                % Energies Correlation
                correlation_index.energy(time_i, d_i) = degreeCorrelation(A.energies, degree.energies);

            % plot the graphs,of pairs, and on it highlight the points where the degree, importance etc properties of the graph is out of common
            % update the markov cluster, and the one from the book -> maybe for the
            % distances graph and clustering, a change in the reference system must
            % be imposed: meaning that I identify two large clusters(let's say) and
            % from here I change the reference systems of all the particles wrt to
            % the centre of mass of each cluster such that m the positions are ow
            % wirt to the cluster, and each new position, is introduced again into
            % the clustering method - the thing is, if you watch carefully the
            % distance and the contact graphs, in boths some of the points are a bit
            % secluded 

                % % Forces Common Correlation
                correlation_index.forces_common = degreeCorrelation(A.forces, degree.forces_common);

                % % Forces Total Common Correlation
                correlation_index.forces_total_common = degreeCorrelation(A.forces_total, degree.forces_total_common);

                % % Pairs Correlation
                correlation_index.pairs = degreeCorrelation(A.pairs, degree.pairs);

                % % Distance Correlation
                correlation_index.distance = degreeCorrelation(A.distance, degree.distance);

            %% Do self loops to avoid division by zero

                 % % Normalise everyhting - Method 1
                 for i = 1:sizes.bodies
                     for j = 1:sizes.bodies
                         A.energies_norm1(i,j) = norm(E.t(i)/sum(abs(E.t)) - E.t(j)/sum(abs(E.t)));
                         A.ang_momentum_norm1(i,j) = norm(L_norm(i)/sum(L_norm) - L_norm(j)/sum(L_norm));
                         A.energies_norm2(i,j) = norm(E.t(i)/max(abs(E.t)) - E.t(j)/max(abs(E.t)));
                         A.ang_momentum_norm2(i,j) = norm(L_norm(i)/max(L_norm) - L_norm(j)/max(L_norm));
                     end
                 end
                 A.everything1 =  A.energies_norm1 + A.ang_momentum_norm1 + A.pairs;
                 A.everything2 = A.energies_norm2 + A.ang_momentum_norm2 + A.pairs;

                 degree.everything1 = sum(A.everything1,2)/sum(sum(A.everything1));
                 G.everything1 = graph(A.everything1);
                 centralities.everything1.closeness(time_i,:) = centrality(G.everything1,"closeness","Cost",G.everything1.Edges.Weight);
                 centralities.everything1.eigenvector(time_i,:) = centrality(G.everything1,"eigenvector","Importance",G.everything1.Edges.Weight, "MaxIterations",100,'Tolerance',1e-4);
                 % centralities.everything1.pagerank(time_i,:) = centrality(G.everything1,"pagerank","Importance",G.everything1.Edges.Weight);
                 centralities.everything1.degree(time_i,:) = centrality(G.everything1,"degree","Importance",G.everything1.Edges.Weight);
                 centralities.everything1.betweenness(time_i,:) = centrality(G.everything1,"betweenness","Cost",G.everything1.Edges.Weight);
           
                % % plotcentralities(positions,centralities.everything1.pagerank(time_i,:),1,3,time_vec(time_i) * 300/3600,'Evr pagerank',doc_i)
                % plotcentralities(positions,centralities.everything1.closeness(time_i,:),2,3,time_vec(time_i) * 300/3600,'Evr Closeness',doc_i)
                % plotcentralities(positions,centralities.everything1.eigenvector(time_i,:),3,3,time_vec(time_i) * 300/3600,'Evr Eigenvector',doc_i)
                % plotcentralities(positions,centralities.everything1.degree(time_i,:),4,3,time_vec(time_i) * 300/3600,'Evr Degree',doc_i)
                % plotcentralities(positions,centralities.everything1.betweenness(time_i,:),5,3,time_vec(time_i) * 300/3600,'Evr Betweeness 2',doc_i)
                % % 
                % plotdegree(positions,degree.everything1*1e3,1,3,time_vec(time_i),'Degree Everything 1',doc_i)
                % 
                % degree.everything2 = sum(A.everything2,2);
                % G.everything2 = graph(A.everything2);
                % centralities.everything2.closeness(time_i,:) = centrality(G.everything2,"closeness","Cost",G.everything2.Edges.Weight);
                % centralities.everything2.eigenvector(time_i,:) = centrality(G.everything2,"eigenvector","Importance",G.everything2.Edges.Weight, "MaxIterations",100,'Tolerance',1e-4);
                % % centralities.everything2.pagerank(time_i,:) = centrality(G.everything2,"pagerank","Importance",G.everything2.Edges.Weight);
                % centralities.everything2.degree(time_i,:) = centrality(G.everything2,"degree","Importance",G.everything2.Edges.Weight);
                % centralities.everything2.betweenness(time_i,:) = centrality(G.everything2,"betweenness","Cost",G.everything1.Edges.Weight);

            %     % plotcentralities(positions,centralities.everything2.pagerank(time_i,:),6,3,time_vec(time_i),'Evr pagerank 2',doc_i)
            %     plotcentralities(positions,centralities.everything2.closeness(time_i,:),7,3,time_vec(time_i) * 300/3600,'Evr Closeness 2',doc_i)
            %     plotcentralities(positions,centralities.everything2.eigenvector(time_i,:),8,3,time_vec(time_i) * 300/3600,'Evr Eigenvector 2',doc_i)
            %     plotcentralities(positions,centralities.everything2.degree(time_i,:),9,3,time_vec(time_i) * 300/3600,'Evr Degree 2',doc_i)
            %     plotcentralities(positions,centralities.everything2.betweenness(time_i,:),10,3,time_vec(time_i) * 300/3600,'Evr Betweeness 2',doc_i)
            % 
            %     plotdegree(positions,degree.everything2*1e3,2,3,time_ii,'Degree Everything 2',doc_i)
            %    % 
            %    % % Centralities
            %    % MATLAB  implementation
            % 
            %     centralities.distance.closeness(time_i,:) = centrality(G.distances,"closeness",'Cost',G.distances.Edges.Weight);
            %     centralities.pairs.closeness(time_i,:) = centrality(G.pairs,"closeness");
            %     centralities.forces.closeness(time_i,:) = centrality(G.forces,"closeness");
            %     centralities.forces_total.closeness(time_i,:) = centrality(G.forces_total,"closeness",'Cost',G.forces_total.Edges.Weight);
            %     centralities.energies.closeness(time_i,:) = centrality(G.energies,"closeness",'Cost',G.energies.Edges.Weight);
            %     centralities.omega.closeness(time_i,:) = centrality(G.omega,"closeness",'Cost',G.omega.Edges.Weight);
            %     % centralities.potentials.closeness(time_i,:) = centrality(G.potentials,"closeness",'Cost',G.potentials.Edges.Weight);
            %     % centralities.ang_momentum.closeness(time_i,:) = centrality(G.ang_momentum,"closeness",'Cost',G.ang_momentum.Edges.Weight);
            %     % centralities.kinetics.closeness(time_i,:) = centrality(G.kinetics,"closeness",'Cost',G.kinetics.Edges.Weight);
            %     % centralities.inertia.closeness(time_i,:) = centrality(G,inertia,'closeness');
            %     % 2)Betweenness centralities - no useful in the disconnected graphs
            % 
            %     centralities.distance.betweenness(time_i,:) = centrality(G.distances, 'betweenness','Cost',G.distances.Edges.Weight);
            %     centralities.pairs.betweenness(time_i,:) = centrality(G.pairs, 'betweenness');
            %     centralities.forces.betweenness(time_i,:) = centrality(G.forces, 'betweenness','Cost',G.forces.Edges.Weight);
            %     centralities.forces_total.betweenness(time_i,:) = centrality(G.forces_total, 'betweenness','Cost',G.forces_total.Edges.Weight);
            %     centralities.energies.betweenness(time_i,:) = centrality(graph(A.energies.*[A.distance <= 20]), 'betweenness','Cost',graph(A.energies.*[A.distance <= 20]).Edges.Weight);
            %     centralities.omega.betweenness(time_i,:) = centrality(G.omega,'betweenness','Cost',G.omega.Edges.Weight);
            %     centralities.potentials.betweenness(time_i,:) = centrality(G.potentials,'betweenness','Cost',G.potentials.Edges.Weight); 
            %     centralities.ang_momentum.betweenness(time_i,:) = centrality(G.ang_momentum,'betweenness','Cost',G.ang_momentum.Edges.Weight); 
            % 
            % %    3) Eigenvector centralities    
            % 
            %     centralities.distance.eigenvector(time_i,:) = centrality(G.distances,'eigenvector','Importance',G.distances.Edges.Weight, "MaxIterations",100,'Tolerance',1e-4);
            %     centralities.pairs.eigenvector(time_i,:) = centrality(G.pairs,'eigenvector', "MaxIterations",100,'Tolerance',1e-4);
            %     centralities.forces.eigenvector(time_i,:) = centrality(G.forces,'eigenvector','Importance',G.forces.Edges.Weight, "MaxIterations",100,'Tolerance',1e-4);
            %     centralities.forces_total.eigenvector(time_i,:) = centrality(G.forces_total,'eigenvector','Importance',G.forces_total.Edges.Weight, "MaxIterations",100,'Tolerance',1e-4);
            %     centralities.energies.eigenvector(time_i,:) = centrality(G.energies,'eigenvector','Importance',G.energies.Edges.Weight, "MaxIterations",100,'Tolerance',1e-4);
            %     centralities.omega.eigenvector(time_i,:) = centrality(G.omega,'eigenvector','Importance',G.omega.Edges.Weight, "MaxIterations",100,'Tolerance',1e-4);
            %     centralities.potentials.eigenvector(time_i,:) = centrality(G.potentials,'eigenvector','Importance',G.potentials.Edges.Weight, "MaxIterations",100,'Tolerance',1e-4);
            %     centralities.ang_momentum.eigenvector(time_i,:) = centrality(G.ang_momentum,'eigenvector','Importance',G.ang_momentum.Edges.Weight, "MaxIterations",100,'Tolerance',1e-4);
            %     centralities.kinetics.eigenvector(time_i,:) = centrality(G.kinetics,'eigenvector','Importance',G.kinetics.Edges.Weight, "MaxIterations",100,'Tolerance',1e-4);
            %     % 4 ) Pagerank centralities
            %     %%
            %     centralities.distance.pagerank(time_i,:) = centrality(G.distances,'pagerank', 'Importance', G.distances.Edges.Weight);
            %     centralities.pairs.pagerank(time_i,:) = centrality(G.pairs,'pagerank');
            %     centralities.forces.pagerank(time_i,:) = centrality(G.forces,'pagerank', 'Importance', G.forces.Edges.Weight);
            %     centralities.forces_total.pagerank(time_i,:) = centrality(G.forces_total,'pagerank', 'Importance', G.forces_total.Edges.Weight);
            %     centralities.energies.pagerank(time_i,:) = centrality(G.energies,'pagerank', 'Importance', G.energies.Edges.Weight);
            %     centralities.omega.pagerank(time_i,:) = centrality(G.omega,'pagerank', 'Importance', G.omega.Edges.Weight);
            %     centralities.potentials.pagerank(time_i,:) = centrality(G.potentials,'pagerank', 'Importance', G.potentials.Edges.Weight); 
            %     centralities.ang_momentum.pagerank(time_i,:) = centrality(G.ang_momentum,'pagerank', 'Importance', G.ang_momentum.Edges.Weight); 
            %     centralities.kinetics.pagerank(time_i,:) = centrality(G.kinetics,'pagerank', 'Importance', G.kinetics.Edges.Weight); 
            %     % 5) Degree centralities 
            %     %%
            %     centralities.distance.degree(time_i,:) = centrality(G.distances,'degree','Importance',G.distances.Edges.Weight);
            %     centralities.pairs.degree(time_i,:) = centrality(G.pairs,'degree');
            %     centralities.forces.degree(time_i,:) = centrality(G.forces,'degree','Importance',G.forces.Edges.Weight);
            %     centralities.forces_total.degree(time_i,:) = centrality(G.forces_total,'degree','Importance',G.forces_total.Edges.Weight);
            %     centralities.energies.degree(time_i,:) = centrality(G.energies,'degree','Importance',G.energies.Edges.Weight);
            %     centralities.omega.degree(time_i,:) = centrality(G.omega,'degree','Importance',G.omega.Edges.Weight);
            %     centralities.potentials.degree(time_i,:) = centrality(G.potentials,'degree','Importance',G.potentials.Edges.Weight); 
            %     centralities.ang_momentum.degree(time_i,:) = centrality(G.ang_momentum,'degree','Importance',G.ang_momentum.Edges.Weight); 
            %     centralities.kinetics.degree(time_i,:) = centrality(G.kinetics,'degree','Importance',G.kinetics.Edges.Weight); 
            % %%    % Plot the data
            % 
            %     plotcentralities(positions,centralities.distance.pagerank(time_i,:),1,1,time_vec(time_i) * 300/3600,'Distance pagerank',doc_i)
            %     plotcentralities(positions,centralities.distance.closeness(time_i,:),2,1,time_vec(time_i) * 300/3600,'Distance Closeness',doc_i)
            %     plotcentralities(positions,centralities.distance.eigenvector(time_i,:),3,1,time_vec(time_i) * 300/3600,'Distance Eigenvector',doc_i)
            %     plotcentralities(positions,centralities.distance.degree(time_i,:),4,1,time_vec(time_i) * 300/3600,'Distance Degree',doc_i)
            %     plotcentralities(positions,centralities.distance.betweenness(time_i,:),5,1,time_vec(time_i) * 300/3600,'Distance Betweenness',doc_i)
            % %
            %     plotcentralities(positions,centralities.pairs.pagerank(time_i,:),6,1,time_vec(time_i) * 300/3600,'Pairs pagerank',doc_i)
            %     plotcentralities(positions,centralities.pairs.closeness(time_i,:),7,1,time_vec(time_i) * 300/3600,'Pairs Closeness',doc_i)
            %     plotcentralities(positions,centralities.pairs.eigenvector(time_i,:),8,1,time_vec(time_i) * 300/3600,'Pairs Eigenvector',doc_i)
            %     plotcentralities(positions,centralities.pairs.degree(time_i,:),9,1,time_vec(time_i) * 300/3600,'Pairs Degree',doc_i)
            %     plotcentralities(positions,centralities.pairs.betweenness(time_i,:),10,1,time_vec(time_i) * 300/3600,'Pairs Betweenness',doc_i)
            % 
            %     plotcentralities(positions,centralities.ang_momentum.pagerank(time_i,:) * 100,11,1,time_vec(time_i) * 300/3600,'Ang Momentum pagerank',doc_i)
            %     plotcentralities(positions,centralities.ang_momentum.closeness(time_i,:) * 100,12,1,time_vec(time_i) * 300/3600,'Ang Momentum Closeness',doc_i)
            %     plotcentralities(positions,centralities.ang_momentum.eigenvector(time_i,:) * 100,13,1,time_vec(time_i) * 300/3600,'Ang Momentum Eigenvector',doc_i)
            %     plotcentralities(positions,centralities.ang_momentum.degree(time_i,:) * 100,14,1,time_vec(time_i) * 300/3600,'Ang Momentum Degree',doc_i)
            %     plotcentralities(positions,centralities.ang_momentum.betweenness(time_i,:),15,1,time_vec(time_i) * 300/3600,'Ang Momentum Betweenness',doc_i)
            % 
            %     saveas(gcf,append(photo_folder,append('Centralities1 ',string(time_vec(time_i)),'.png')))
            % %%
            %     plotcentralities(positions,centralities.energies.pagerank(time_i,:),1,2,time_vec(time_i) * 300/3600,'Energies pagerank',doc_i)
            %     plotcentralities(positions,centralities.energies.closeness(time_i,:),2,2,time_vec(time_i) * 300/3600,'Energies Closeness',doc_i)
            %     plotcentralities(positions,centralities.energies.eigenvector(time_i,:),3,2,time_vec(time_i) * 300/3600,'Energies Eigenvector',doc_i)
            %     plotcentralities(positions,centralities.energies.degree(time_i,:),4,2,time_vec(time_i) * 300/3600,'Energies Degree',doc_i)
            %     plotcentralities(positions,centralities.energies.betweenness(time_i,:),5,2,time_vec(time_i) * 300/3600,'Energies Betweenness',doc_i)
            % 
            %     plotcentralities(positions,centralities.omega.pagerank(time_i,:),6,2,time_vec(time_i) * 300/3600,'Angular Velocity pagerank',doc_i)
            %     plotcentralities(positions,centralities.omega.closeness(time_i,:),7,2,time_vec(time_i) * 300/3600,'Angular Velocity Closeness',doc_i)
            %     plotcentralities(positions,centralities.omega.eigenvector(time_i,:),8,2,time_vec(time_i) * 300/3600,'Angular Velocity Eigenvector',doc_i)
            %     plotcentralities(positions,centralities.omega.degree(time_i,:),9,2,time_vec(time_i) * 300/3600,'Angular Velocity Degree',doc_i)
            %     plotcentralities(positions,centralities.omega.betweenness(time_i,:),10,2,time_vec(time_i) * 300/3600,'Ang Velocity Betweenness',doc_i)
            % %
            %     plotcentralities(positions,centralities.potentials.pagerank(time_i,:) * 100,11,2,time_vec(time_i) * 300/3600,'Potentials pagerank',doc_i)
            %     plotcentralities(positions,centralities.potentials.closeness(time_i,:) * 100,12,2,time_vec(time_i) * 300/3600,'Potentials Closeness',doc_i)
            %     plotcentralities(positions,centralities.potentials.eigenvector(time_i,:) * 100,13,2,time_vec(time_i) * 300/3600,'Potentials Eigenvector',doc_i)
            %     plotcentralities(positions,centralities.potentials.degree(time_i,:) * 100,14,2,time_vec(time_i) * 300/3600,'Potentials Degree',doc_i)
            %     plotcentralities(positions,centralities.potentials.betweenness(time_i,:),15,2,time_vec(time_i) * 300/3600,'Potentials Betweenness',doc_i)
            %     saveas(gcf,append(photo_folder,append('Centralities2 ',string(time_vec(time_i) * 300/3600),'.png')))
            %     hold off;
            % %%
            %     plotcentralities(positions,centralities.forces.pagerank(time_i,:) * 100,1,3,time_vec(time_i) * 300/3600,'Forces pagerank',doc_i)
            %     plotcentralities(positions,centralities.forces.closeness(time_i,:) * 100,2,3,time_vec(time_i) * 300/3600,'Forces Closeness',doc_i)
            %     plotcentralities(positions,centralities.forces.eigenvector(time_i,:) * 100,3,3,time_vec(time_i) * 300/3600,'Forces Eigenvector',doc_i)
            %     plotcentralities(positions,centralities.forces.degree(time_i,:) * 100,4,3,time_vec(time_i) * 300/3600,'Forces Degree',doc_i)
            %     plotcentralities(positions,centralities.forces.betweenness(time_i,:),5,3,time_vec(time_i) * 300/3600,'Forces Betweenness',doc_i)
            %     saveas(gcf,append(photo_folder,append('Centralities2 ',string(time_vec(time_i) * 300/3600),'.png')))

                %
                % Calcualte this data also for each cluster
                %%
                % plotdegree1(positions, degree.pairs, time_vec(time_i) * 300/3600, doc_i,1)
                % plotdegree1(positions, degree.energies, time_ii, doc_i,2)
                % plotdegree1(positions, degree.distance, time_ii, doc_i,3)


            
                % plotdegree(positions,degree.pairs,1,1,time_ii,'Pairs',doc_i)
                % plotdegree(positions,degree.distance,2,1,time_ii,'Distance',doc_i)
                % plotdegree(positions,degree.forces_total,3,1,time_ii,'Forces Total',doc_i)
                % plotdegree(positions,degree.potential,4,1,time_ii,'Potential',doc_i)
                % plotdegree(positions,degree.omega,5,1,time_ii,"Angular Velocity",doc_i)
                % plotdegree(positions,degree.energies,6,1,time_ii,'Energies',doc_i)
                % plotdegree(positions,degree.ang_momentum,7,1,time_ii,'Angular Momentum',doc_i)
                % plotdegree(positions,degree.forces,8,1,time_ii,'Forces',doc_i)
                % plotdegree(positions,degree.kinetics,9,1,time_ii,'Kinetic Energy',doc_i)
                % saveas(gcf,append(photo_folder,append('Degrees ',string(time_vec(time_i) * 300/3600),'.png')))
                % 
                % plotdegree(positions,degree.distance_common,1,2,time_ii,'Distance Common',doc_i)
                % plotdegree(positions,degree.pairs_common,2,2,time_ii,'Pairs Common',doc_i)
                % plotdegree(positions,degree.forces_total_common,6/2,2,time_ii,'Forces Total Common',doc_i)
                % plotdegree(positions,degree.energies_common,8/2,2,time_ii,'Energies Common',doc_i)
                % plotdegree(positions,degree.omega_common,10/2,2,time_ii,"Angular Velocity Common",doc_i)
                % plotdegree(positions,degree.potentials_common,12/2,2,time_ii,'Potential Common',doc_i)
                % plotdegree(positions,degree.ang_momentum_common,7,2,time_ii,'Angular Momentum Common',doc_i)
                % plotdegree(positions,degree.forces_common,8,2,time_ii,'Forces Common',doc_i)
                % plotdegree(positions,degree.kinetics_common,9,2,time_ii,'Kinetic Energy Common',doc_i)
                % saveas(gcf,append(photo_folder,append('Degrees Common ',string(time_vec(time_i)),'.png')))
                % 
                % hold off;

                % Calculate the Entropy of the system
                % Degree-Based Graph Entropy - Based on different methodollogy and equations
                
                  H.pairs(time_i,r_i,d_i) = calculate_entropy(A.pairs,"normal");
                 H.distance(time_i,r_i,d_i) = calculate_entropy(A.distance,"normal");
                 H.force(time_i,r_i,d_i) = calculate_entropy(A.forces,"normal");
                 H.energies(time_i,r_i,d_i) = calculate_entropy(A.energies,"normal");
                 H.forces_total(time_i,r_i,d_i) = calculate_entropy(A.forces_total,"normal");
                 H.pairs_von(time_i,r_i,d_i) = calculate_entropy(A.pairs,"Von Newmann");
                 H.distance_von(time_i,r_i,d_i) = calculate_entropy(A.distance,"Von Newmann");
                 H.force_von(time_i,r_i,d_i) = calculate_entropy(A.forces,"Von Newmann");
                 H.energies_von(time_i,r_i,d_i) = calculate_entropy(A.energies,"Von Newmann");
                 H.forces_total_von(time_i,r_i,d_i) = calculate_entropy(A.forces_total,"Von Newmann");
                 H.pairs_lambda(time_i,r_i,d_i) = calculate_entropy(A.pairs,"random walker");
                 H.distance_lambda(time_i,r_i,d_i) = calculate_entropy(A.distance,"random walker");
                 H.force_lambda(time_i,r_i,d_i) = calculate_entropy(A.forces,"random walker");
                 H.energies_lambda(time_i,r_i,d_i) = calculate_entropy(A.energies,"random walker");
                 H.forces_total_lambda(time_i,r_i,d_i) = calculate_entropy(A.forces_total,"random walker");
             
                 
             


            %%
                % figure(100)
                % hold on;
                % subplot(1,3,1)
                % hold on;
                % plot(1:sizes.bodies,omega(:,1));
                % subplot(1,3,2)
                % plot(1:sizes.bodies,omega(:,2));
                % subplot(1,3,3)
                % plot(1:sizes.bodies,omega(:,3));
                % grid on;

                %% Compute the centre of mass of the system

                  r_com = sum(Mass.all .* positions) / Mass.total ;
                  velocity_com = sum(Mass.all .* velocities) /Mass.total;
              
                  L = zeros(sizes.bodies,3);
                  L_each = zeros(sizes.bodies,3);
                  Pot_amd_each = zeros(sizes.bodies,1);
                  I_each = zeros(sizes.bodies,1);
            %     % Compute the angular momentum of each particle
                  for i = 1:sizes.bodies
                      L(i,1:3) = cross((positions(i,:) - r_com),Mass.all(i)*(velocities(i,:) - velocity_com));
                  end
              
                  Ltotal(time_i) = norm(sum(L,1));
                  L_total_vec(time_i,1:3) = sum(L,1);
                  % Moment of inertia of the entire body % If it is polar with respect to
                    CoM of initial system
                  I_prime(time_i,d_i) = sum(Mass.all.* vecnorm(positions,2,2).^2); 
                  I_prime_com(time_i,d_i) = sum(Mass.all.* vecnorm(positions - r_com,2,2).^2); 
            %     %%
                  H_angular_vector = zeros(sizes.bodies,3);
                  for i = 1:sizes.bodies
                      H_angular_vector(i,:) = Mass.all(i) * cross(positions(i,:),velocities(i,:)); 
                  end
              
            %     % Total angular momentum vector
              H_vec = sum(H_angular_vector, 1);
              H_hat = H_vec / norm(H_vec);  % Unit vector of H
              
            % % Construct total inertia tensor from particle data (example)
              I_tensor = zeros(3,3);
              for i = 1:sizes.bodies
                  r_rel = positions(i,:) - r_com;
                  m_i = Mass.all(i);
                  r_cross = [   0      -r_rel(3)  r_rel(2);
                             r_rel(3)     0     -r_rel(1);
                            -r_rel(2)  r_rel(1)     0    ];
                  I_tensor = I_tensor + m_i * (norm(r_rel)^2 * eye(3) - (r_rel') * r_rel);
              end
              
            % % Compute I_h = Ĥᵀ * I * Ĥ
              I_h(time_i,d_i) = H_hat * I_tensor * H_hat';
              
              % Amended potential
              U = sum(E.pe);  % Assuming negative potential already
              Pot_amd1(time_i,d_i) = norm(H_vec)^2 / (2 * I_h(time_i,d_i)) + U;
              
                  %
                  I_prime_h(time_i,d_i) = H_hat * diag(sum(inertia,1)) * H_hat';
            %     % Amanded potential 
                 Pot_amd(time_i,d_i) = Ltotal(time_i) / I_prime_h(time_i) /2 + sum(E.pe);
              tic
            %     % Calculate the above data for each particle
                  for i = 1:sizes.bodies
                      L_each(i,:) = cross((positions(i,:) - r_com),Mass.all(i)*(velocities(i,:) - velocity_com));
                      I_each(i) = 0.4 * Mass.all(i) * norm(positions(i,:))^2;
                      Pot_amd_each(i) = norm(L_each(i,:)) / I_each(i) / 2 + E.pe(i);
                  end
              
              
                  % % Angular m omentum of the entire system
                  L = [0,0,0];
                  for i = 1:sizes.bodies
                      L = L + cross(positions(i,:),Mass.all(i) * velocities(i,:));
                  end
                  L_total(time_i,1:3) = L; 
              
                  %% Percolation
              
                  A.pairs = pairs.time_All(:,:,time_i);
              
                  k = sum(A.pairs,2);
                  den = sum(0.5 * k .* (k-1));
              
                  num = trace(A.pairs^3);
            %     % Get the degree distribution
                  degree.pairs_common = sum(A.pairs,2); 
                  degree_graphs = degree.pairs_common;
                  [max_k,~] = max(degree_graphs);
                  [counts, edges] = histcounts(k, 0:max_k + 1,'Normalization','probability');
              
                  k_vals = 0:max_k;
                  k1 = sum(k_vals .* counts);
              
                  second_moment = sum(k_vals.^2 .* counts);
                  mean_degree = sum(degree_graphs)/sizes.bodies;
              
                     k = sum(A.pairs,2);
                  den = sum(0.5 * k .* (k-1));
                  num = trace(A.pairs^3);
              
                  C_global(time_i,d_i,r_i) = num /den;   % global clustering coefficient 
            % 
            %     % Critical Probability percolation
                  pc(time_i, d_i, r_i) = 1/(1-C_global(time_i,d_i,r_i)) * (mean_degree)/(second_moment - mean_degree);
                  pc_c(time_i, d_i, r_i) = (mean_degree)/(second_moment - mean_degree);
              
                  centers = edges(1:end-1) + 0.5;
                  pc_pois(time_i, d_i, r_i) = 1/mean_degree;
              
              
            %     %% Generating Function
                  p_c(time_i, d_i, r_i) =  generatingfunction(G.pairs);
              
                A.pairs = pairs.time_All(:,:,time_i);
                C_global(time_i,d_i,r_i) = num /den; % global clustering coefficient 
            %   % Critical Probability percolation 
              
                if (C_global(time_i, d_i, r_i) <= 0.1)
                    pc_active(time_i, d_i, r_i) = 1/(1-C_global(time_i,d_i,r_i)) * generatingfunction(G.pairs);
                else
                     pc_active(time_i, d_i, r_i) =  generatingfunction(G.pairs);
                end
              
                 toc
                %% Let's the evolution of pairs, angular momentum, energy

                % Adj{time_i} = A;

                %% Save the cluster values of the time i
                clustervalues_timei{time_i} = clustervalue;
                close all;


                % Potential Energy time
                E.pe_time(time_i,d_i) = sum(E.pe);
                end 

                hold off;

                %% 
                %    % Angular Momentum
                    figure(24072000 + doc_i)
                    semilogy(time_vec,nonzeros(Ltotal(:,d_i)),'LineWidth', 2.5, 'LineStyle',linestyle{d_i})
                    % title("Total Angular Momentum")
                    ylabel('$L$ $[kg \cdot m^2/s]$','Interpreter','latex','FontSize',16)
                    xlabel('Time [hrs]','Interpreter','latex','FontSize',16)
                    xticks([0,24,48,72,96])
                    xticklabels({'0','2','4','6','8'});

                    grid on;
                    hold on;
                    % plot(time_vec,Ltotal(1)*0.99 * ones(1,length(time_vec)),'LineWidth',2.5,'LineStyle','--')
                    % plot(time_vec,Ltotal(1)*1.01 * ones(1,length(time_vec)),'LineWidth',2.5,'LineStyle','--')
                    legend(["A","B","C","D","E","F"],'interpreter','latex','FontSize',16)
                    saveas(gcf,append(photo_folder,'Angular Momentum.png'))

                    figure(25072000 + doc_i)
                    semilogy(time_vec, nonzeros(I_prime(:,d_i)),'LineWidth',2.5, 'LineStyle',linestyle{d_i})
                    ylabel('$I_h$ $[kg \cdot m^2]$','Interpreter','latex','FontSize',16)
                    xlabel('Time [hrs]','Interpreter','latex','FontSize',16)
                    legend(["A","B","C","D","E","F"],'interpreter','latex','FontSize',16)
                    yl = ylim;         % Get current y-axis limits
                    ylim([5e11 1e15]); 
                    xticks([0,24,48,72,96])
                    xticklabels({'0','2','4','6','8'});

                    hold on;
                % title("Total moment of inertia")
                    saveas(gcf,append(photo_folder,'Moment of Inertia.svg'))

                    figure(26072000 + doc_i)
                    semilogy(time_vec, nonzeros(Pot_amd(:,d_i))/1e3,'LineWidth',2.5, 'LineStyle',linestyle{d_i})
                %  title("Total Amanded Potential")
                    saveas(gcf,append(photo_folder,'Amended Potential .svg'))
                    ylabel('$\varepsilon$ [kJ]','Interpreter','latex','FontSize',16)
                    xlabel('Time [hrs]','Interpreter','latex','FontSize',16)
                    legend(["A","B","C","D","E","F"],'interpreter','latex','FontSize',16)
                    xticks([0,24,48,72,96])
                    xticklabels({'0','2','4','6','8'});
                    grid off;
                    hold on;

                %% Calculate the moment of inertia of the entire system and the angular
                % momentum and the angular velocity
                    figure(27072000 + doc_i)
                    hold on;
                    for fig_i = 1:size(E.pe_time,2)
                    plot(time_vec, E.pe_time(:,fig_i)/1000, 'LineWidth', 2.5, 'LineStyle',linestyle{fig_i})
                    end
                    ylabel('$U$ [kJ]','Interpreter','latex','FontSize',16)
                    ylim([-1300 -100])
                    xlabel('Time [hrs]','Interpreter','latex','FontSize',16)
                    legend(["A","B","C","D","E","F"],'interpreter','latex','FontSize',16,'Location','northeastoutside')
                    xticks([0,24,48,72,96])
                    xticklabels({'0','2','4','6','8'});
                    grid off;
                    hold off;
%%
            end
        end
    end
end

%    save(append("C:/Users/mihne/Documents/GitHub/Chrono_Projects/files/",core{core_i},formatast(1:end-5),"_",density_folder(5:8),"_",rotation_folder(end-5:end-1),".mat"));
    %%

%%
markerstyle = {"o", "x", "square", "v", "pentagram"};

%% Plot section
figure()
hold on;
if (numel(densitiesall) > 1)
    for i = 1:size(pc_active,2)
        pc_save(:,i) = pc_active(:,i,i);
        pp(i) = plot(time_vec,smoothdata(pc_save(:,i), "gaussian", 6),'LineWidth',3, "LineStyle", linestyle{i});
    end
else
    for i = 1:size(pc_active,3)
     pc_save(:,i) = pc_active(:,:,i);
     pp(i) = plot(time_vec,smoothdata(pc_save(:,i), "gaussian", 6),'LineWidth',3, "LineStyle", linestyle{i});
    end
end
hold off;
grid on;

xlabel('Time [h]','Interpreter','latex','FontSize',18)
ylabel('Percolation Threshold $p_c$[-]','Interpreter','latex','FontSize',18)
ylim([0.15 1.65])
xticks([0,24,48,72,96])
xticklabels({'0','2','4','6','8'});
% if (numel(densitiesall) > 1)
% % legend([pp(2), pp(3), pp(1), pp(4), pp(5), pp(6)], ...
% %        {'$\rho = 1200 \ \mathrm{kg/m^3}$', '$\rho = 1600 \ \mathrm{kg/m^3}$', '$\rho = 2000 \ \mathrm{kg/m^3}$', '$\rho = 2400 \ \mathrm{kg/m^3}$', '$\rho = 2800 \ \mathrm{kg/m^3}$', '$\rho = 3000 \ \mathrm{kg/m^3}$'}, ...
% %        'Interpreter', 'latex','FontSize',18);
% legend(["Case A", "Case B", "Case C", "Case D", "Case E"], 'Interpreter', 'latex','FontSize',18);
%  % legend([pp(2), pp(3),pp(1),pp(4),pp(5),pp(6)], ...
%  %     {"D.","D.","D.","D.","D.","R."}, 'Interpreter','latex','FontSize',18,'Orientation','horizontal','Location','northoutside');
if (numel(densitiesall) > 1)
 legend([pp(1), pp(2), pp(3), pp(4), pp(5)], ...
        {"Case A", "Case B", "Case C", "Case D", "Case E"}, ...
        'Interpreter', 'latex','FontSize',18);
 %       {'$\rho = 1200 \ \mathrm{kg/m^3}$', '$\rho = 1600 \ \mathrm{kg/m^3}$', '$\rho = 2000 \ \mathrm{kg/m^3}$', '$\rho = 2400 \ \mathrm{kg/m^3}$', '$\rho = 2800 \ \mathrm{kg/m^3}$', '$\rho = 3000 \ \mathrm{kg/m^3}$'}, ...
        
else
 legend([pp(2), pp(3), pp(4), pp(1), pp(5), pp(6)],...
    {"$P = 21.8$ hrs", "$P = 5$ hrs","$P = 3.2$ hrs", "$P = 2.5$ hrs", "$P = 2.1$ hrs", "$P = 1.9$ hrs" }, 'Interpreter','latex','FontSize',18);
end

%%
cmap = cool;       % Default gives 64 colors
firstColor = cmap(1, :);
lastColor = cmap(end,:);
figure()
hold on;
scatter(1:10, sin(linspace(0,pi,10)),"filled","MarkerFaceColor",firstColor);
scatter(1:10, cos(linspace(0,pi,10)),"filled","MarkerFaceColor",lastColor);
legend({"Small Node Degree", "High Node Degree"},'Orientation','horizontal','Location','northoutside','Interpreter','latex','FontSize',18)
hold off;

%% Plot the pcs
% Rotation 
figure() 
hold on;
pc_save = zeros(numel(time_vec), numel(densitiesall));
for i = 1:size(pc_active,2)
    pc_save(:,i)  = pc_active(:,i,i);
    pp(i) = plot(time_vec,pc_save(:,i),'LineWidth',3, "LineStyle", linestyle{i});
end


hold off;
grid on;

xlabel('Time [h]','Interpreter','latex','FontSize',16)
ylabel('Percolation Threshold $p_c$[-]','Interpreter','latex','FontSize',16)
xticks([3,24,48,72,96])
xticklabels({'0.25','2','4','6','8'});
legend([pp(1), pp(2), pp(3), pp(4), pp(5)], ...
    {"Case A", "Case B", "Case C", "Case D", "Case E"}, ...
    'Interpreter', 'latex','FontSize',18);

 %% Plot cluster global

 figure() 
hold on;
for i = 1:size(pc_active,2)
    C_global_save(:,i)  = C_global(:,i,i);
    pp(i) = plot(time_vec,smooth(C_global_save(:,i)),'LineWidth',3,"LineStyle", linestyle{i});
end
hold off;
grid on;
xlabel('Time [h]','Interpreter','latex','FontSize',16)
ylabel('Global Clustering Coefficient $c$[-]','Interpreter','latex','FontSize',16)
xticks([0,24,48,72,96])
xticklabels({'0','2','4','6','8'});
legend([pp(1), pp(2), pp(3), pp(4), pp(5)], ...
    {"Case A", "Case B", "Case C", "Case D", "Case E"}, ...
    'Interpreter', 'latex','FontSize',18);

%% Plot Pcs density cases
figure()
plot(time_vec, pc,'LineWidth',3)
xlabel('Time[H]','Interpreter','latex','FontSize',16)
ylabel('Percolation Threshold $p_c$[-]','Interpreter','latex','FontSize',16)
xticks([0.25,24,48,72,96])
xticklabels({'0.25','2','4','6','8'});
legend({'$\rho = 2000 \ \mathrm{kg/m^3}$','$\rho = 1200 \ \mathrm{kg/m^3}$','$\rho = 1600 \ \mathrm{kg/m^3}$','$\rho = 2400 \ \mathrm{kg/m^3}$','$\rho = 2800 \ \mathrm{kg/m^3}$'},'Interpreter','latex');
%% Plot the averages
figure()
hold on;
for i = 1:size(avg_degree.pairs,2)
    pp(i) = plot(time_vec, smooth(avg_degree.pairs(:,i,i)),'LineWidth',3,"LineStyle", linestyle{i});
end


grid on;
hold off;
xlabel('Time [h]','Interpreter','latex','FontSize',18)
ylabel('Average Degree of the system [-]','Interpreter','latex','FontSize',18)
xticks([0,24,48,72,96])
xticklabels({'0','2','4','6','8'});
legend([pp(1), pp(2), pp(3), pp(4), pp(5)], ...
    {"Case A", "Case B", "Case C", "Case D", "Case E"}, ...
    'Interpreter', 'latex','FontSize',18);

    %% Plot Graph Entryopy - lambda
figure()
hold on;
for i = 1:size(H.pairs,2)
    H.pairs2(:,i) = H.pairs_lambda(:,i,i);
    pp(i) = plot(time_vec,smoothdata(H.pairs2(:,i), "gaussian", 6),'LineWidth',3, "LineStyle", linestyle{i});
end
hold off;
grid on;

xlabel('Time [h]','Interpreter','latex','FontSize',18)
ylabel('Entropy Pairs $h_P$[-]','Interpreter','latex','FontSize',18)
xticks([3,24,48,72,96])
xticklabels({'0.25','2','4','6','8'});
legend([pp(1), pp(2), pp(3), pp(4), pp(5)], ...
    {"Case A", "Case B", "Case C", "Case D", "Case E"}, ...
    'Interpreter', 'latex','FontSize',18);


figure()
hold on;
for i = 1:size(H.force,2)
    H.force2(:,i) = H.force_lambda(:,i,i);
    pp(i) = plot(time_vec,smoothdata(H.force2(:,i), "gaussian", 6),'LineWidth',3, "LineStyle", linestyle{i});
end


xlabel('Time [h]','Interpreter','latex','FontSize',18)
ylabel('Entropy Force $h_F$[-]','Interpreter','latex','FontSize',18)
xticks([3,24,48,72,96])
xticklabels({'0.25','2','4','6','8'});
legend([pp(1), pp(2), pp(3), pp(4), pp(5)], ...
    {"Case A", "Case B", "Case C", "Case D", "Case E"}, ...
    'Interpreter', 'latex','FontSize',18);

figure()
hold on;
if (numel(densitiesall) > 1)
    for i = 1:size(H.energies,2)
        H.energies2(:,i) = H.energies_lambda(:,i,i);
        pp(i) = plot(time_vec,smoothdata(H.energies2(:,i), "gaussian", 6),'LineWidth',3, "LineStyle", linestyle{i});
    end
else
    for i = 1:size(pc_active,3)
     pc_save(:,i) = pc_active(:,:,i);
     pp(i) = plot(time_vec,smoothdata(pc_save(:,i), "gaussian", 6),'LineWidth',3, "LineStyle", linestyle{i});
    end
end
hold off;
grid on;

xlabel('Time [h]','Interpreter','latex','FontSize',18)
ylabel('Entropy Energy $h_E$[-]','Interpreter','latex','FontSize',18)
xticks([3,24,48,72,96])
xticklabels({'0.25','2','4','6','8'});
legend([pp(1), pp(2), pp(3), pp(4), pp(5)], ...
    {"Case A", "Case B", "Case C", "Case D", "Case E"}, ...
    'Interpreter', 'latex','FontSize',18);

figure()
hold on;
for i = 1:size(H.pairs,2)
    H.pairs1(:,i) = H.pairs(:,i,i);
    pp(i) = plot(time_vec,smoothdata(H.pairs1(:,i), "gaussian", 4),'LineWidth',3, "LineStyle", linestyle{i});
end


xlabel('Time [h]','Interpreter','latex','FontSize',18)
ylabel('Entropy Pairs $h_P$[-]','Interpreter','latex','FontSize',18)
xticks([0,24,48,72,96])
xticklabels({'0','2','4','6','8'});
legend([pp(1), pp(2), pp(3), pp(4), pp(5)], ...
    {"Case A", "Case B", "Case C", "Case D", "Case E"}, ...
    'Interpreter', 'latex','FontSize',18);

figure()
hold on;
for i = 1:size(H.force,2)
    H.force1(:,i) = H.force(:,i,i);
    pp(i) = plot(time_vec,smoothdata(H.force1(:,i), "gaussian", 4),'LineWidth',3, "LineStyle", linestyle{i});
end
xlabel('Time [h]','Interpreter','latex','FontSize',18)
ylabel('Entropy Force $h_F$[-]','Interpreter','latex','FontSize',18)
xticks([0,24,48,72,96])
xticklabels({'0','2','4','6','8'});
legend([pp(1), pp(2), pp(3), pp(4), pp(5)], ...
    {"Case A", "Case B", "Case C", "Case D", "Case E"}, ...
    'Interpreter', 'latex','FontSize',18);

figure()
hold on;
if (numel(densitiesall) > 1)
    for i = 1:size(H.energies,2)
        H.energies1(:,i) = H.energies(:,i,i);
        pp(i) = plot(time_vec,smoothdata(H.energies1(:,i), "gaussian", 4),'LineWidth',3, "LineStyle", linestyle{i});
    end
else
    for i = 1:size(pc_active,3)
     pc_save(:,i) = pc_active(:,:,i);
     pp(i) = plot(time_vec,smoothdata(pc_save(:,i), "gaussian", 4),'LineWidth',3, "LineStyle", linestyle{i});
    end
end
hold off;
grid on;

xlabel('Time [h]','Interpreter','latex','FontSize',18)
ylabel('Entropy Energy $h_E$[-]','Interpreter','latex','FontSize',18)
xticks([0,24,48,72,96])
xticklabels({'0','2','4','6','8'});
legend([pp(1), pp(2), pp(3), pp(4), pp(5)], ...
    {"Case A", "Case B", "Case C", "Case D", "Case E"}, ...
    'Interpreter', 'latex','FontSize',18);


figure()
hold on;
for i = 1:size(H.energies,2)
    H.distance(:,i) = H.distance(:,i,i);
    pp(i) = plot(time_vec,smoothdata(H.distance(:,i), "gaussian", 4),'LineWidth',3, "LineStyle", linestyle{i});
end
hold off;
grid on;

xlabel('Time [h]','Interpreter','latex','FontSize',18)
ylabel('Entropy Distance $h_D$[-]','Interpreter','latex','FontSize',18)
xticks([0,24,48,72,96])
xticklabels({'0','2','4','6','8'});
legend([pp(1), pp(2), pp(3), pp(4), pp(5)], ...
    {"Case A", "Case B", "Case C", "Case D", "Case E"}, ...
    'Interpreter', 'latex','FontSize',18);



%%
figure()
hold on;
for i = 1:size(H.pairs,2)
    H.pairs5(:,i) = H.pairs_von(:,i,i);
    pp(i) = plot(time_vec,smoothdata(H.pairs5(:,i), "gaussian", 4),'LineWidth',3, "LineStyle", linestyle{i});
end


xlabel('Time [h]','Interpreter','latex','FontSize',18)
ylabel('Entropy Pairs $h_P$[-]','Interpreter','latex','FontSize',18)
xticks([0,24,48,72,96])
xticklabels({'0','2','4','6','8'});
if (numel(densitiesall) > 1)
legend([pp(1), pp(2), pp(3), pp(4), pp(5)], ...
    {"Case A", "Case B", "Case C", "Case D", "Case E"}, ...
    'Interpreter', 'latex','FontSize',18);


figure()
hold on;
for i = 1:size(H.force,2)
    H.force5(:,i) = H.force_von(:,i,i);
    pp(i) = plot(time_vec,smoothdata(H.force5(:,i), "gaussian", 4),'LineWidth',3, "LineStyle", linestyle{i});
end

xlabel('Time [h]','Interpreter','latex','FontSize',18)
ylabel('Entropy Force $h_F$[-]','Interpreter','latex','FontSize',18)
xticks([0,24,48,72,96])
xticklabels({'0','2','4','6','8'});
legend([pp(1), pp(2), pp(3), pp(4), pp(5)], ...
    {"Case A", "Case B", "Case C", "Case D", "Case E"}, ...
    'Interpreter', 'latex','FontSize',18);

figure()
hold on;
for i = 1:size(H.energies,2)
    H.energies5(:,i) = H.energies_von(:,i,i);
    pp(i) = plot(time_vec,smoothdata(H.energies5(:,i), "gaussian", 4),'LineWidth',3, "LineStyle", linestyle{i});
end
hold off;
grid on;

xlabel('Time [h]','Interpreter','latex','FontSize',18)
ylabel('Entropy Energy $h_E$[-]','Interpreter','latex','FontSize',18)
xticks([0,24,48,72,96])
xticklabels({'0','2','4','6','8'});
legend([pp(1), pp(2), pp(3), pp(4), pp(5)], ...
    {"Case A", "Case B", "Case C", "Case D", "Case E"}, ...
    'Interpreter', 'latex','FontSize',18);


figure()
hold on;
for i = 1:size(H.energies,2)
    H.distance(:,i) = H.distance(:,i,i);
    pp(i) = plot(time_vec,smoothdata(H.distance(:,i), "gaussian", 4),'LineWidth',3, "LineStyle", linestyle{i});
end

hold off;
grid on;

xlabel('Time [h]','Interpreter','latex','FontSize',18)
ylabel('Entropy Distance $h_D$[-]','Interpreter','latex','FontSize',18)
xticks([0,24,48,72,96])
xticklabels({'0','2','4','6','8'});
legend([pp(1), pp(2), pp(3), pp(4), pp(5)], ...
    {"Case A", "Case B", "Case C", "Case D", "Case E"}, ...
    'Interpreter', 'latex','FontSize',18);


%% Check the clusters coefficient

clusters = clusterpairs(pairs.time_All(:,:,96));

kcl1 = sum(cluster1,2);
dencl1 = sum(0.5 * kcl1 .* (kcl1-1));
numcl1 = trace(cluster1^3);

C_clusters1 = numcl1 /dencl1;

kcl2 = sum(cluster2,2);
dencl2 = sum(0.5 * kcl2 .* (kcl2-1));
numcl2 = trace(clusters2^3);

C_clusters2 = numcl2 /dencl2;


%% ================================================================================
%% ================================================================================
%% Functions

function clusters = clusterpairs(A, delta)
    % Function to find clusters based on adjacency matrix A.
    % Each row represents an entity, and a '1' indicates a connection.
    
    clusters = {};  % Initialize empty cell array for clusters
    
    for i = 1:size(A, 1)  % Loop over each row
        contacts = find(A(i, :) == 1);  % Find all connected elements in row i
        contacts = [i, contacts];  % Include the current row in the cluster
        
        % Check if any contact is already in an existing cluster
        found = false;
        for j = 1:length(clusters)
            if any(ismember(contacts, clusters{j}))  % If any contact is in an existing cluster
                clusters{j} = unique([clusters{j}, contacts]);  % Merge into existing cluster
                found = true;
                break;  % Exit loop after merging
            end
        end
        
        % If no existing cluster found, create a new one
        if ~found
            clusters{end+1} = contacts;
        end
    end

    for clu_i = length(clusters):-1:2
        for clu_j = clu_i-1 :-1:1
            if any(ismember(clusters{clu_i}, clusters{clu_j}))  % If any contact is in an existing cluster
                clusters{clu_j} = unique([clusters{clu_j}, clusters{clu_i}]); 
                clusters{clu_i} = [];
            end
        end
    end
    for i = length(clusters):-1:1
        if (numel(clusters{i}) == 0)
            clusters(i) = [];
        end
    end
end


function plotmarkov(clusters_mark,Acase,positions,points2clusters,name,time_vec,time_i,photo_folder)
    % Function to plot the results obtained from applying the method Markov cluster
            % hex_colors = { "#FF0000";  "#00FF00";   "#0000FF";   "#00FFFF";  "#FF00FF";   "#FFFF00";  "#000000";   "#FFFFFF"; "#0072BD";   "#D95319";  "#EDB120";   "#7E2F8E"; "#77AC30";   "#4DBEEE";  "#A2142F"; "#1A2B3C"; "#4F5E6D"; "#A1B2C3"; "#D4E5F6";"#FF5733"; "#33FF57";"#5733FF"; "#C0C0C0";"#800000"; "#008000";"#000080"; "#FFA500";"#4B0082"; "#EE82EE";"#4682B4"; "#20B2AA";"#DC143C"; "#8B0000";"#556B2F2"; "#2F4F4F"};
            % % Colors = { Red, Green, Blue, Cyan, Magenta, Yellow, Black, White, Dark, Blue, Dark Orange, Dark Yellow, Dark Purple, Medium Green, Light Blue, Dark Red};
            figure();
            hold on;
            colors = lines(length(clusters_mark));  % Generate distinct colors for clusters
            colormap('lines')
            scatter3(positions(:,1), positions(:,2), positions(:,3), 50,points2clusters, 'filled', 'o');
            set(gca,"FontSize",14)
            hold off;
            grid on;
            view([90 90 90]);
            xlabel('x [m]', 'Interpreter','latex'); ylabel('y [m]', 'Interpreter','latex'); zlabel('z [m]', 'Interpreter','latex');
            title(sprintf('T = %.2f hrs', round(time_vec(time_i) * 300/ 3600,2)));
            xlabel('x [m]','Interpreter','latex','FontSize',16)
            ylabel('y [m]','Interpreter','latex','FontSize',16)
            zlabel('z [m]','Interpreter','latex','FontSize',16)
          %  title(append('3D Clustered Points ', Acase)); 
            saveas(gcf,append(photo_folder,append(name,string(time_vec(time_i))),'.png'))

end

%% Markov 
function [clusters, points2clusters] = markov_clustering(A,max_iter, exp_power, inflate_power, tol)
        
    % Input : data - the consideed data
    % max_iter - the maximum nr of iterations
    % inflate_power - Inflation step  power (>1, default,2);
    % exp_power - Expansion step power ( default : 2);
    % tol - Convergence tolerance (default: 1e-6)

    adjacency_matrix = A;
    A_norm = adjacency_matrix./sum(adjacency_matrix,1);

    for iter = 1:max_iter
        A_norm_new = A_norm^exp_power;
        A_norm_new = A_norm_new.^inflate_power;
        A_norm_new = A_norm_new./sum(A_norm_new,1);

        if max(abs(A_norm_new(:) - A_norm(:))) < tol
            break;

        end
        A_norm = A_norm_new;
    end

    
       for i = 1:numel(A_norm)
           if(A_norm(i) < 1e-6)
                A_norm(i) = 0;
           end
       end

       
       cl_i = 1;
       for j = 1:size(A_norm,2)
           if any(A_norm(j,:))
            numcl(cl_i) = j;
            cl_i = cl_i + 1;
           end
       end

     % Step 2: Create Graph
    G = graph(adjacency_matrix);

    clusters = cell(1,length(numcl));
   for i = 1:length(numcl)
       clusters{i} = find(A_norm(numcl(i),:));
   end

    S = cellfun(@(x) sprintf('%g_', x), clusters, 'UniformOutput', false);
    [~, idx] = unique(S, 'stable');
    clusters = clusters(idx);

    % Points2Clusters
    [numbodies,~] = size(A); 
    points2clusters = zeros(numbodies,1);
    for i = 1:length(clusters)
        points2clusters(clusters{i}) = i; 
    end

end


function communities = recursive_partition(A)

    % Method that applies a recursive partitioning on a graph 
    communities = {};
    queue = {1:size(A,1)};
    while ~isempty(queue)
        idx = queue{end};
        queue(end) = [];

        A_sub = A(idx, idx);
        m = sum(sum(A_sub))/2;
        k = sum(A_sub, 2);
        B = A_sub - (k * k') / (2*m);
        
        [V, D] = eig(B);
        [~, maxIdx] = max(diag(D));
        v1 = V(:, maxIdx);
        s = sign(v1);
        deltaQ = (1/(4*m)) * s' * B * s;
        
        if deltaQ > 1e-5  % Accept split
            group1 = idx(s > 0);
            group2 = idx(s <= 0);
            queue{end+1} = group1;
            queue{end+1} = group2;
        else
            communities{end+1} = idx;
        end
    end
end


function lambda_k = myLaplacian(A,lim)
    % Function that takes the most important eigenvalue of the network (from the Laplacian)
    D = sum(A,2);
    L = diag(D) - A;
    lambda = sort(eig(L),'ascend');
    lambda_k = lambda(lambda < lim);
end


function value = readdata(path, field)
    % Function that reads a file and assignes the values in the file to variable
    data = readlines(path);
    lengthoffield = length(char(field)); 
    
    for i = 1:length(data)
        data_char = char(data(i));
        if (length(data_char) > lengthoffield)
            if (data_char(1:lengthoffield) == char(field))
                value = str2num(data_char(lengthoffield+1:end));
                break;
            end
        end
    end
end

%% Degree Correlation Function
function r = degreeCorrelation(A, degree)
    [i, j] = find(triu(A));  % Extract edges from the upper triangle
    ki = degree(i);          % Degree of node i
    kj = degree(j);          % Degree of node j
    k_avg = mean(degree);    % Average degree
    numerator = sum((ki - k_avg) .* (kj - k_avg));
    denominator = sum((ki - k_avg).^2);
    r = numerator / denominator;
end

function plotdegree(positions,degree,i,m,figurenr,titlefig,fig_doc)
    % Function that plots the degree distribution (not evolution) different from the other plots
    pointsizes = ((degree(1:end) - min(degree(1:end)))/(max(degree(1:end)) - min(degree(1:end))));
    figure(599 + figurenr + m + fig_doc)
    set(gcf, 'WindowState', 'maximized');
    colormap(cool)
    grid on;
    scatter3(positions(1:end,1),positions(1:end,2),positions(1:end,3),36,degree(1:end),'filled')
    view([1 0.2 1])
    xlabel('x [m]','Interpreter','latex','FontSize',18); ylabel('y [m]','Interpreter','latex','FontSize',18); zlabel('z [m]','Interpreter','latex','FontSize',18);
    %title(sprintf('T = %.2f hrs', round(titlefig,2)));
    c = colorbar;
    c.FontSize = 16;
    ax = gca;
    ax.FontSize = 16; 
    colorbar;
end

function plotdegree1(positions,degree,time_ii,fig_doc,i)
    % Function that plotsdegree using function ploting
    figurenumber = 599 + time_ii + i;
    ploting(positions, degree, figurenumber);
end
function plotcentralities1(positions,centralities,i,figurenr)
    % Function that plots centrality using function ploting
    figurenumber = figurenr+ i + 699;
    ploting(positions,centralities,figurenumber)
end

function ploting(positions, data,figurenumber,titlefig,figure_name)
    % Function that applies a general format of plotting to both centralities and degree plots
    figure();
    colormap(cool)
    grid on;
    scatter3(positions(:,1),positions(:,2),positions(:,3),50,data,'filled')
    grid on;
    xlabel('x [m]','Interpreter','latex','FontSize',16); ylabel('y [m]','Interpreter','latex','FontSize',16); zlabel('z [m]','Interpreter','latex','FontSize',16);
    %title(sprintf('T = %.2f hrs', round(titlefig,2)));
    c = colorbar;
    c.Location = 'eastoutside';
    c.FontSize = 16;
    ax = gca;
    ax.FontSize = 16; 
    view([-1.391841209537282e+02,8.382043566548514]);
end

function plotcentralities(positions,centralities,i,m,figurenr,titlefig,fig_doc)
    % Function that plots the centrality distribution (not evolution) different from the other plots
    pointsizes = ((centralities(1:end) - min(centralities(1:end)))/(max(centralities(1:end)) - min(centralities(1:end)))) + 1e-6;
    centralities = centralities / max(centralities);
    figurenumber = figurenr+ i + 699;
    ploting(positions,centralities,figurenumber,figurenr,titlefig)
    saveas(gcf,append('C:\Users\mihne\Documents\GitHub\Chrono_Projects\files\spherical_large_no_core_fixed\Phtots_thesis\Centralities\',append(titlefig,string(figurenr),'.svg')))
end


function [entropy] = calculate_entropy(A,entropy_case)
% Calculate the entropy of the system
% Calcualte the strength / degree of each edge
entropy = 0;
switch entropy_case
    case "normal"
    p_i = sum(A,2)/sum(A,'all');
    % The following part implements the limit x->0 xlog(x)
    for i = 1:numel(p_i)
        if (p_i(i) == 0)
            continue;
        else
            entropy = entropy - p_i(i) * log(p_i(i));
        end
    end
    case "Von Newmann"
    L = laplacian(graph(A));
    if (trace(L) == 0)
        entropy = 0;
        return;
    end
    rho = L/trace(L);
    lambda = eig(rho);
    for i = 1:numel(lambda)
        if(lambda(i) <= 1e-3)
            continue;
        else
            entropy = entropy - lambda(i) * log(lambda(i));   
        end
    end
    case "random walker"
    % 1. Find the largest eigenvalue of the adjacency matrix
    % For sparse matrices, eigs is much faster than eig
    if isempty(find(A, 1))
        entropy = 0; % No edges = no paths = no entropy
    else
        lambda_max = eigs(sparse(A), 1);
        
        % 2. Entropy is the log of the spectral radius
        % Use log2 for bits, log for nats
        entropy = log2(lambda_max);
    end
end
end


    % Min span depth first breadth first
    % Comment out  if needed because first breadth search is a very expensive algorithm for a network this large
    %  graphstruct_distance.bfirsts = zeros(sizes.bodies,sizes.bodies);
    %  graphstruct_pairs.bfirsts = zeros(sizes.bodies,sizes.bodies);
    %  graphstruct_distance.dfirsts = zeros(sizes.bodies,sizes.bodies);
    %  graphstruct_pairs.bdfirsts = zeros(sizes.bodies,sizes.bodies);
    % for nodes_i = 1:sizes.bodies 
    % %Breadth-first search
    %     a1 = bfsearch(G.distances,nodes_i);
    %    graphstruct_distance.bfirsts(nodes_i,1:numel(a1)) = a1;
    %     a1 = bfsearch(G.pairs,nodes_i);
    %    graphstruct_pairs.bfirsts(nodes_i,1:numel(a1)) = a1;
    % %Depth-first search
    %     a1 = dfsearch(G.distances,nodes_i);
    %    graphstruct_distance.dfirsts(nodes_i, 1:numel(a1)) = a1;
    %     a1 = dfsearch(G.pairs,nodes_i);
    %    graphstruct_pairs.dfirsts(nodes_i,1:numel(a1)) = a1;
    % end

    % % Minimum spannning tree
    % graphstruct.minspan_distances = minspantree(G.distances);
    % graphstruct.minspan_pairs = minspantree(G.pairs);


function p_c = generatingfunction(G)
    % This function creates the generating function of a graph and creates the percolation threshold
    % Input: Graph - G
    % Output: percolation threshold p_c
    d = degree(G); 
    k_values = min(d):max(d); 

    %  Degree distribution
    Pk = histcounts(d, 'Normalization', 'probability');

    syms x;
    g0 = sum(Pk .* x.^k_values);
    g0_p = diff(g0,x);
    g1 = (g0_p/subs(g0_p,x,1)); 
    
    g1_p = matlabFunction(diff(g1));
    p_c = 1/g1_p(1);

end

function [p_axis, S] = karrer_percolation(G)
    % A percolation function similar to generating function but uses a different definition for p_c
    % G is a MATLAB graph object
    Adj = adjacency(G);
    [row, col] = find(Adj); 
    num_edges = length(row); % Number of directed edges
    N = numnodes(G);
    
    p_axis = 0:0.02:1;
    S = zeros(size(p_axis));
    
    % Precompute neighbors to speed up the product
    neighbors = cell(N, 1);
    for i = 1:N
        neighbors{i} = find(Adj(i, :));
    end

    for p_idx = 1:length(p_axis)
        p = p_axis(p_idx);
        % Initialize u_ij (probability edge doesn't lead to giant component)
        u = ones(num_edges, 1) * 0.5; 
        
        % Iterative Message Passing (Fixed Point)
        for iter = 1:100
            u_old = u;
            for e = 1:num_edges
                i = row(e);
                j = col(e);
                
                % Product over neighbors k of i, excluding j
                prod_term = 1;
                nbors = neighbors{i};
                for k = nbors
                    if k ~= j
                        % Find index of directed edge k -> i
                        idx_ki = find(row == k & col == i, 1);
                        prod_term = prod_term * u(idx_ki);
                    end
                end
                u(e) = (1 - p) + p * prod_term;
            end
            if max(abs(u - u_old)) < 1e-6, break; end
        end
        
        % Calculate size of Giant Component S
        node_not_in_GC = zeros(N, 1);
        for i = 1:N
            prod_all = 1;
            for k = neighbors{i}
                idx_ki = find(row == k & col == i, 1);
                prod_all = prod_all * u(idx_ki);
            end
            node_not_in_GC(i) = (1 - p) + p * prod_all;
        end
        
        % S is the fraction of nodes in the giant component
        S(p_idx) = 1 - mean(node_not_in_GC);
    end
    
    % Plotting the result
    plot(p_axis, S, 'LineWidth', 2);
    xlabel('Occupation Probability p');
    ylabel('Giant Component Size S');
    grid on;
end