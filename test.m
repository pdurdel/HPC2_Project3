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


%[bdry_dirichlet_1, bdry_neumann_1] = get_boundary_1(e', bdry_dirichlet_idx, bdry_neumann_idx);
[bdry_dirichlet_u10, bdry_neumann] = get_boundary_2(e', bdry_dirichlet_idx_u10, bdry_neumann_idx);
[bdry_dirichlet_u20, ~] = get_boundary_2(e', bdry_dirichlet_idx_u20);

%bdry_dirichlet_u10
%bdry_dirichlet_u20
%bdry_neumann

%{
bdry_dirichlet_1
bdry_dirichlet_2

bdry_neumann_1
bdry_neumann_2

dir_idx_2 = [bdry_dirichlet_2(:, 1); bdry_dirichlet_2(end, 2)];
neu_idx_2 = [bdry_neumann_2(:, 1);  bdry_neumann_2(end, 2)];

c(dir_idx_2, :)
c(neu_idx_2, :)
%}

%% c)

f = @(location, state) zeros(1, length(location.x));
g = @(location, state) zeros(1, length(location.x));
u_D = @(location, state) u_dirichlet([location.x location.y], coords_circles(1, :), big_radii(1), coords(8, :), coords(9, :), 10^(-13));


u_D_2d = @(x) u_dirichlet(x, coords_circles(1, :), big_radii(1), coords(8, :), coords(9, :), 10^(-14));
%{
rad = 3.5:0.001:4;
angle = 0:0.001:2*pi;
[X, Y] = meshgrid(rad, angle);

[X, Y] = meshgrid(0.25:0.001:1.75, 5.8:0.001:6.2);

Z = zeros(size(X));
for j=1:size(X, 2)
    %Z(:, j) = u_D_2d([coords_circles(1, 1)+X(:, j).*sin(Y(:, j)), coords_circles(1, 2)+X(:, j).*cos(Y(:, j))]);
    Z(:, j) = u_D_2d([X(:, j) Y(:, j)]);
end

%Z(Z>11) = 0;

figure(5)
%surf(coords_circles(1, 1) + X.*sin(Y), coords_circles(1, 2) + X.*cos(Y), Z);
surf(X, Y, Z);
%}




%{
coords(8, :)
coords(9, :)

disp(u_D(coords(1, :)))
disp(u_D(coords(2, :)))
disp(u_D([-6.5 -1.5]))
disp(u_D([8 8]))
disp(u_D(coords(8, :)))
disp(u_D(coords(9, :)))
disp(u_D([1 6]))
%}


%% d)

model = createpde();

geometryFromMesh(model, c', el');
%geometryFromEdges(model, geom);
    
mesh = generateMesh(model, 'hmax', hmax);
[p2, e2, t2] = meshToPet(mesh);
figure(3)
pdegplot(model, 'EdgeLabels', 'on')




applyBoundaryCondition(model, 'dirichlet', 'Edge', bdry_dirichlet_idx, 'u', u_D)
%applyBoundaryCondition(model, 'dirichlet', 'Edge', bdry_dirichlet_u10(:, 1), 'u', 10)
%applyBoundaryCondition(model, 'dirichlet', 'Edge', bdry_dirichlet_u20(:, 1), 'u', 20)

applyBoundaryCondition(model, 'neumann', 'Edge', bdry_neumann(:, 1), 'g', g)
specifyCoefficients(model, m=0, d=0, c=1, a=0, f=f);

u_sol_model = solvepde(model);
u_sol = u_sol_model.NodalSolution;

figure(4)
pdeplot(model, 'XYData', u_sol, 'ZData', u_sol)

%% e)

