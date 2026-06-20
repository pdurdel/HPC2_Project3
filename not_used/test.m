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

%figure(1)
%pdegplot(geom, 'EdgeLabels', 'on', 'VertexLabels', 'on');

el = t(1:3,:)';
c = p';

%figure(2)
%trisurf(el, c(:,1), c(:,2), 0.*c(:,2), 'edgecolor','k'), view(2)


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

figure(4)
pdeplot(model, 'XYData', u_sol, 'ZData', u_sol)


%% e)

