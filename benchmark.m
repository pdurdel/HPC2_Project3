clear
close all
clc

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


%% tractor shape mesh initialization

[coords, coords_circles, big_radii, small_radii] = define_tractor_parameters;

for i = 1:length(hmax)

    % p: point coordinates
    % e: edges
    % t: elements
    
    [p, e, t, geom] = get_mesh_from_points(coords, coords_circles, big_radii, small_radii, hmax(i));
    

    %% function handles
    
    f = @(location, state) zeros(1, length(location.x));
    g = @(location, state) zeros(1, length(location.x));
    u_D = @(location, state) u_dirichlet([location.x location.y], coords_circles(1, :), big_radii(1), coords(8, :), coords(9, :), 10^(-13));
    

    %% boundary conditions
    
    [bdry_dirichlet_idx, bdry_neumann] = define_boundary_conditions(e);
    
    
    %% solve Poisson problem
    
    model = createpde();
    
    geometryFromEdges(model, geom);
    mesh = generateMesh(model, 'hmax', hmax(i), 'GeometricOrder', 'linear');
    
    
    applyBoundaryCondition(model, 'dirichlet', 'Edge', bdry_dirichlet_idx, 'u', u_D)
    applyBoundaryCondition(model, 'neumann', 'Edge', bdry_neumann(:, 1).', 'g', g)
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
        I_u(i, j) = line_integral(s, u_line, L);

        % snapshot for the bottom plots
        if i == i_snap && j == j_snap
            geom_snap = geom;
            p_snap = p;
            t_snap = t;
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



%% Plots

%% model geometry with edge labels

figure(1)
pdegplot(model_snap, 'EdgeLabels', 'on')
xlabel('x'), ylabel('y')
title('PDE Model Geometry')

%% tractor mesh

figure(2)
trisurf(t_snap, p_snap(:,1), p_snap(:,2), 0.*p_snap(:,2), ...
    'edgecolor', 'k'), view(2)
xlabel('x'), ylabel('y')
axis equal
title('Triangular Mesh')

%% 3d FEM solution

figure(3)
pdeplot(model_snap, 'XYData', u_sol_snap, 'ZData', u_sol_snap)
grid on
xlabel('x'), ylabel('y'), zlabel('u')
title('FEM Solution u')

%% line points evaluation

figure(4)
plot(s_snap, u_line_snap, 'LineWidth', 2), grid on;
xlabel('s'), ylabel('u')
title('Function Evaluation along the Line')

%% directional derivative along the line

figure(5)
plot(s_snap, beta_grad_u_snap, 'LineWidth', 2), grid on;
xlabel('s'), ylabel('\beta \cdot \nabla u')
title('Directional Derivative along the Line')

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