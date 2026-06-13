clc;
close all;
clear;

%% a)

coords = [
    (-13-3*sqrt(6))/2 3
    -3.5 0
    -1.5 0
    1.5 0
    2.5 0
    2.5 3
    1.5 3
    1.5 6
    0.5 6
    0.5 3
    -3.5 5
    -3.5 8
    -7.5 8
];

coords_circles = [
    -6.5 2.25
    0 0
];

big_radii = [3.75 1.5];
small_radii = [1 0.5];
hmax = 0.5;

[p,e,t,geom] = get_mesh_from_points(coords, coords_circles, big_radii, small_radii, hmax);

figure(1)
pdegplot(geom, 'EdgeLabels', 'on', 'VertexLabels', 'on');

el = t(1:3,:)';
c = p';

figure(2)
trisurf(el, c(:,1), c(:,2), 0.*c(:,2), 'edgecolor','k'), view(2)


%% b)

bdry_dirichlet_idx_u10 = [5 12];
bdry_dirichlet_idx_u20 = 14:17;
bdry_dirichlet_idx = [bdry_dirichlet_idx_u10 bdry_dirichlet_idx_u20];
bdry_neumann_idx = [1:4 6:11 13];

% Using get_boundary_3
bdry_idcs = {bdry_dirichlet_idx_u10, bdry_dirichlet_idx_u20, bdry_neumann_idx};
bdry = get_boundary_3(e', bdry_idcs, 17);

bdry_dirichlet_u10 = bdry{1};
bdry_dirichlet_u20 = bdry{2};
bdry_dirichlet = [bdry_dirichlet_u10; bdry_dirichlet_u20];
bdry_neumann = bdry{3};


%% c)

f = @(location, state) zeros(1, length(location.x));
g = @(location, state) zeros(1, length(location.x));
u_D = @(location, state) u_dirichlet([location.x location.y], coords_circles(1, :), big_radii(1), coords(8, :), coords(9, :), 10^(-13));


%% d)

model = createpde();

geometryFromEdges(model, geom);
mesh = generateMesh(model, 'hmax', hmax, 'GeometricOrder', 'linear');
%{
[p2, e2, t2] = meshToPet(mesh);
figure(3)
pdegplot(model, 'EdgeLabels', 'on')
%}

applyBoundaryCondition(model, 'dirichlet', 'Edge', bdry_dirichlet_idx, 'u', u_D)
applyBoundaryCondition(model, 'neumann', 'Edge', bdry_neumann(:, 1).', 'g', g)
specifyCoefficients(model, m=0, d=0, c=1, a=0, f=f);

u_sol_model = solvepde(model);
u_sol = u_sol_model.NodalSolution;

%figure(4)
%pdeplot(model, 'XYData', u_sol, 'ZData', u_sol)
%grid on;

%% e)

% blue line VW
V = [-2.0 -1.0];
W = [3.0 1.0];

% direction vector (v->w)
beta = (W-V) / norm(W-V);

% creating points on VW
n_line = 200;
s = linspace(0,1,n_line)';
line_points = (1-s).*V + s.*W;

% getting mesh from pde
mesh_points = model.Mesh.Nodes.';
mesh_elements = model.Mesh.Elements(1:3,:).';

% eval u on th eline
u_line = artificial_point_eval(mesh_points, mesh_elements, u_sol, line_points);

% eval gradient on the line
grad_u_line = artificial_point_grad(mesh_points, mesh_elements, u_sol, line_points);

% direc. deviation
beta_grad_u = grad_u_line * beta.';

% plots
figure(5)
plot(s, u_line, 'LineWidth', 2), grid on;
xlabel('s'), ylabel('u'); title('function evaluation');

figure(6)
plot(s, beta_grad_u, 'LineWidth', 2), grid on;
xlabel('s'), ylabel('\beta \cdot \nabla u');
title('directional derivation');

figure(7)
pdeplot(model, 'XYData', u_sol), hold on;
plot(line_points(:,1), line_points(:,2), 'b-', 'LineWidth', 2)
plot(V(1), V(2), 'ko', 'MarkerFaceColor', 'k')
plot(W(1), W(2), 'ko', 'MarkerFaceColor', 'k')
text(V(1), V(2), ' V')
text(W(1), W(2), ' W')
colormap(turbo)
colorbar
axis equal; grid on;
title('plot of the line');


%% f)

rand_points = 20;

% getting mesh data
mesh_points = model.Mesh.Nodes.';
mesh_elements = model.Mesh.Elements(1:3,:).';

% creating the random points
random_points = artificial_point_rand(mesh_points, mesh_elements, rand_points);

% eval u for the random points
u_random = artificial_point_eval(mesh_points, mesh_elements, u_sol, random_points);

figure(8)
clf

pdeplot(model,'XYData', u_sol,'ZData', u_sol);
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
for i = 1:rand_points
    text(random_points(i,1) + 0.08, ...
         random_points(i,2) + 0.08, ...
         u_random(i) + 0.15, ...
         num2str(u_random(i), '%.2f'), ...
         'FontSize', 8, ...
         'Color', 'k');
end

axis equal, grid on, view(3);
xlabel('x'), ylabel('y'), zlabel('u');
title('evaluation of random points')

hrot = rotate3d(gcf);
hrot.Enable = 'on';