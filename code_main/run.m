% init
clear; close all; clc
t_start = tic;
warning('off','all')
addpath('functions')

%% user input

% hmax = [1, 1/2, 1/4, 1/8, 1/16, 1/32];
hmax = [1, 1/2, 1/4, 1/8, 1/16];

% blue line VW
V = [-2.0 -1.0];
W = [3.0 1.0];

n_line = [4, 8, 16, 32, 64, 128, 256, 512];

% snapshot used for the plots
hmax_snap = 1/4;
n_line_snap = 512;
i_snap = find(hmax == hmax_snap, 1);
j_snap = find(n_line == n_line_snap, 1);

% number of timing repetitions
n_reps = 100;

% number of random points to be evaluated
n_rand = 20;

% postprocessing time per (hmax, n_line)
t_post = zeros(length(hmax), length(n_line));

% line integrals (QoIs for the convergence study)
L = norm(W - V);
I_u = zeros(length(hmax), length(n_line));
mean_u_line = zeros(length(hmax), length(n_line));


%% tractor shape mesh initialization

[coords, coords_circles, big_radii, small_radii] = define_tractor_parameters;

for i = 1:length(hmax)

    % p: point coordinates
    % e: edges
    % t: elements
    
    [p, e, t, geom] = get_geometry(coords, coords_circles, big_radii, small_radii, hmax(i));
    

    %% function handles
    
    f = @(location, state) zeros(1, length(location.x));
    g = @(location, state) zeros(1, length(location.x));
    u_D = @(location, state) u_dirichlet([location.x location.y], coords_circles(1, :), big_radii(1), coords(8, :), coords(9, :), 10^(-13));
    

    %% boundary conditions
    
    [bdry_dirichlet_idx, bdry_dirichlet, bdry_neumann_idx, bdry_neumann] = define_boundary_conditions(e);
    
    
    %% solve Poisson problem
    
    model = createpde();
    
    geometryFromEdges(model, geom);
    mesh = generateMesh(model, 'hmax', hmax(i), 'GeometricOrder', 'linear');
    
    applyBoundaryCondition(model, 'dirichlet', 'Edge', bdry_dirichlet_idx, 'u', u_D);
    applyBoundaryCondition(model, 'neumann', 'Edge', bdry_neumann_idx, 'g', g);
    specifyCoefficients(model, m=0, d=0, c=1, a=0, f=f);
    
    u_sol_model = solvepde(model);
    u_sol = u_sol_model.NodalSolution;
        
    %% blue line analysis
    for j = 1:length(n_line)
        
        %% sample points along the line
        
        % direction vector (v->w)
        beta = (W - V) / norm(W - V);
        
        % creating points on VW
        s = linspace(0, 1, n_line(j))';
        line_points = (1 - s).*V + s.*W;
        
        %% evaluate u and its gradient on the line
        
        % getting mesh from pde
        mesh_points = model.Mesh.Nodes.';
        mesh_elements = model.Mesh.Elements(1:3, :).';
        
        % time the postprocessing over n_reps repetitions
        timer = tic;
        for rep = 1:n_reps

            % eval u on the line
            u_line = artificial_point_eval(mesh_points, mesh_elements, u_sol, line_points);

            % eval gradient on the line
            grad_u_line = artificial_point_grad(mesh_points, mesh_elements, u_sol, line_points);

            % directional derivative
            beta_grad_u = grad_u_line * beta.';

        end
        t_post(i, j) = toc(timer) / n_reps;

        % line integral (QoI for the convergence study)
        [I_u(i, j), mean_u_line(i, j)] = line_integral(s, u_line, L);

        % snapshot for the bottom plots
        if i == i_snap && j == j_snap
            geom_snap = geom;
            p_snap = p;
            t_snap = t;
            bdry_dirichlet_snap = bdry_dirichlet;
            bdry_neumann_snap = bdry_neumann;
            model_snap = model;
            u_sol_snap = u_sol;
            s_snap = s;
            line_points_snap = line_points;
            u_line_snap = u_line;
            beta_grad_u_snap = beta_grad_u;
        end

    end
end


%% artificial point evaluation with random points

% getting mesh data from the snapshot model
mesh_points_snap = model_snap.Mesh.Nodes.';
mesh_elements_snap = model_snap.Mesh.Elements(1:3, :).';

% creating the random points
random_points = artificial_point_rand(mesh_points_snap, ...
    mesh_elements_snap, n_rand);

% eval u for the random points
u_random = artificial_point_eval(mesh_points_snap, ...
    mesh_elements_snap, u_sol_snap, random_points);

runtime = toc(t_start);
fprintf('\nruntime: %.6f sec\n', runtime);

%% Plots

% output folder
out_dir = "test_output";

if ~exist(out_dir, "dir")
    mkdir(out_dir);
end

%% model geometry with edge labels

figure(1)
pdegplot(model_snap, 'EdgeLabels', 'on')
xlabel('x'), ylabel('y')
title('PDE Model Geometry')

exportgraphics(gcf, fullfile(out_dir, "01_model_geometry.png"), "Resolution", 300);
savefig(gcf, fullfile(out_dir, "01_model_geometry.fig"));

%% tractor mesh

figure(2)
p1 = trisurf(t_snap, p_snap(:,1), p_snap(:,2), 0.*p_snap(:,2), ...
    'edgecolor', 'k');
hold on
p2 = plot(reshape(p_snap(bdry_dirichlet_snap(:, 2:3),1),[],2)', ...
reshape(p_snap(bdry_dirichlet_snap(:, 2:3),2),[],2)','b', 'LineWidth', 3);
p3 = plot(reshape(p_snap(bdry_neumann_snap(:, 2:3),1),[],2)', ...
reshape(p_snap(bdry_neumann_snap(:, 2:3),2),[],2)','r', 'LineWidth', 3);
view(2)
xlabel('x'), ylabel('y')
axis equal
legend([p1(1), p2(1), p3(1)], 'Mesh', 'Dirichlet Boundary', 'Neumann Boundary')
title('Triangular Mesh')

exportgraphics(gcf, fullfile(out_dir, "02_triangular_mesh.png"), "Resolution", 300);
savefig(gcf, fullfile(out_dir, "02_triangular_mesh.fig"));

%% 3d FEM solution

figure(3)
pdeplot(model_snap, 'XYData', u_sol_snap, 'ZData', u_sol_snap)
colormap(turbo)
colorbar
grid on
xlabel('x'), ylabel('y'), zlabel('u')
title('FEM Solution u')

exportgraphics(gcf, fullfile(out_dir, "03_fem_solution_3d.png"), "Resolution", 300);
savefig(gcf, fullfile(out_dir, "03_fem_solution_3d.fig"));

%% line points evaluation

figure(4)
plot(s_snap, u_line_snap, 'LineWidth', 2), grid on;
xlabel('s'), ylabel('u')
title('Function Evaluation along the Line')

exportgraphics(gcf, fullfile(out_dir, "04_line_evaluation.png"), "Resolution", 300);
savefig(gcf, fullfile(out_dir, "04_line_evaluation.fig"));

%% directional derivative along the line

figure(5)
plot(s_snap, beta_grad_u_snap, 'LineWidth', 2), grid on;
xlabel('s'), ylabel('\beta \cdot \nabla u')
title('Directional Derivative along the Line')

exportgraphics(gcf, fullfile(out_dir, "05_directional_derivative.png"), "Resolution", 300);
savefig(gcf, fullfile(out_dir, "05_directional_derivative.fig"));
%% PDE solution with displayed line

figure(6)
pdeplot(model_snap, 'XYData', u_sol_snap), hold on;
plot(line_points_snap(:,1), line_points_snap(:,2), 'r-', 'LineWidth', 2)
plot(V(1), V(2), 'ko', 'MarkerFaceColor', 'k')
plot(W(1), W(2), 'ko', 'MarkerFaceColor', 'k')
text(V(1), V(2), ' V')
text(W(1), W(2), ' W')
colormap(turbo)
colorbar
axis equal; grid on;
hold off
xlabel('x'), ylabel('y')
title('Solution with the Line VW')

exportgraphics(gcf, fullfile(out_dir, "06_solution_with_line.png"), "Resolution", 300);
savefig(gcf, fullfile(out_dir, "06_solution_with_line.fig"));

%% postprocessing time vs. line resolution

figure(7)
hold on
for i = 1:length(hmax)
    plot(n_line, t_post(i, :), '-o', 'LineWidth', 1.5)
end
hold off
grid on
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('number of line points $n_{line}$', 'Interpreter', 'latex')
ylabel('postprocessing time in s')
legend(arrayfun(@(h) sprintf('h_{max} = 1/%d', round(1/h)), hmax, 'UniformOutput', false), 'Location', 'northwest')
title('Postprocessing Time over Line Resolution')

exportgraphics(gcf, fullfile(out_dir, "07_postprocessing_time.png"), "Resolution", 300);
savefig(gcf, fullfile(out_dir, "07_postprocessing_time.fig"));


%% convergence of the line integral

col = length(n_line);   % largest n_line -> negligible quadrature error
ref = length(hmax);     % finest mesh = reference
h = hmax(1:ref-1).';

err_u = abs(I_u(1:ref-1, col) - I_u(ref, col));

% observed convergence order
p_u = polyfit(log(h), log(err_u), 1);

figure(8)
loglog(h, err_u, '-o', 'LineWidth', 1.5)
grid on
xlabel('mesh size $h_{max}$', 'Interpreter', 'latex')
ylabel('error of the line integral')
legend(sprintf('\\int u ds  (order %.2f)', p_u(1)), 'Location', 'southeast')
title(['Convergence of the Line Integral with $n_{line}$ = ' num2str(n_line(end))], 'Interpreter', 'latex')

exportgraphics(gcf, fullfile(out_dir, "08_convergence_line_integral.png"), "Resolution", 300);
savefig(gcf, fullfile(out_dir, "08_convergence_line_integral.fig"));

%% random points

figure(9)
clf

pdeplot(model_snap, 'XYData', u_sol_snap, 'ZData', u_sol_snap);
hold on

colormap(turbo)
colorbar

% plotting the artificial points
plot3(random_points(:,1), ...
      random_points(:,2), ...
      u_random, ...
      'ko', ...
      'MarkerFaceColor', 'w', ...
      'MarkerSize', 6)

% plotting the corresponding values
for i = 1:n_rand
    text(random_points(i,1) + 0.08, ...
         random_points(i,2) + 0.08, ...
         u_random(i) + 0.15, ...
         num2str(u_random(i), '%.2f'), ...
         'FontSize', 8, ...
         'Color', 'k');
end

axis equal, grid on, view(3);
xlabel('x'), ylabel('y'), zlabel('u');
title('Evaluation of Random Points')
hold off

hrot = rotate3d(gcf);
hrot.Enable = 'on';

exportgraphics(gcf, fullfile(out_dir, "09_random_points.png"), "Resolution", 300);
savefig(gcf, fullfile(out_dir, "09_random_points.fig"));

%% print result tables

fprintf('\n\n================ LINE INTEGRAL VALUES ================\n');
fprintf('Format: line integral (line integral of the valid length)\n\n');

fprintf('%12s', '');
for j = 1:length(n_line)
    fprintf('%28s', "n_line_" + string(n_line(j)));
end
fprintf('\n');

fprintf('%12s', '');
for j = 1:length(n_line)
    fprintf('%28s', repmat('-', 1, 24));
end
fprintf('\n');

% rows
for i = 1:length(hmax)

    fprintf('%12s', "hmax_" + string(hmax(i)));

    for j = 1:length(n_line)

        entry = sprintf('%.3f (%.3f)', I_u(i,j), mean_u_line(i,j));
        fprintf('%28s', entry);

    end

    fprintf('\n');
end


fprintf('\n\n================ RANDOM POINT EVALUATION ================\n\n');

% table for random point coordinates and evaluated values
random_point_table = table( ...
    random_points(:,1), ...
    random_points(:,2), ...
    u_random, ...
    'VariableNames', {'x', 'y', 'u_value'});

disp(random_point_table);