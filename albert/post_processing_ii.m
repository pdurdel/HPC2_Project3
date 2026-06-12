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

big_radii = [3.75; 1.5];
small_radii = [1; 0.5];
hmax = 0.5;

[p,e,t,geom] = get_geometry(coords, coords_circles, big_radii, small_radii, hmax);


figure(1)
pdegplot(geom, 'EdgeLabels', 'on', 'VertexLabels', 'on');
title('Geometry of the Shape')


c_g = p';
el_g = t(1:3,:)';


figure(2)
trisurf(el_g, c_g(:,1), c_g(:,2), 0.*c_g(:,2), 'edgecolor','k'), view(2)
title(sprintf('Triangulation of the Shape with hmax=%f', hmax))


%% b)

bdry_dirichlet_idx_u10 = [5 12];
bdry_dirichlet_idx_u20 = 14:17;
bdry_dirichlet_idx = [bdry_dirichlet_idx_u10 bdry_dirichlet_idx_u20];
bdry_neumann_idx = [1:4 6:11 13];

% Using get_boundary_3
bdry_idcs = {bdry_dirichlet_idx_u10, bdry_dirichlet_idx_u20, bdry_neumann_idx};
bdry = get_boundary(e', bdry_idcs, 17);

bdry_dirichlet_u10 = bdry{1};
bdry_dirichlet_u20 = bdry{2};
bdry_dirichlet = [bdry_dirichlet_u10; bdry_dirichlet_u20];
bdry_neumann = bdry{3}(:, 1);


%% c)

f = @(location, state) zeros(1, length(location.x));
g = @(location, state) zeros(1, length(location.x));
u_D = @(location, state) u_dirichlet([location.x location.y], coords_circles(1, :), big_radii(1), coords(8, :), coords(9, :), 10^(-13));


%% d)

model = createpde();

geometryFromEdges(model, geom);
mesh = generateMesh(model, 'hmax', hmax);
[p, e, t] = meshToPet(mesh);

applyBoundaryCondition(model, 'dirichlet', 'Edge', bdry_dirichlet_idx, 'u', u_D);
applyBoundaryCondition(model, 'neumann', 'Edge', bdry_neumann, 'g', g);
specifyCoefficients(model, m=0, d=0, c=1, a=0, f=f);

u_sol_model = solvepde(model);
u_sol = u_sol_model.NodalSolution;


figure(3)
pdeplot(model, 'XYData', u_sol, 'ZData', u_sol)
title('Solution via SolvePDE')


%% e)  und f)

c = p';
el = t(1:3,:)';
bucket_size = hmax;

[x_min, x_max, y_min, y_max] = get_limits(coords, coords_circles, big_radii);
el_buckets = get_buckets(c, el, bucket_size, x_min, x_max, y_min, y_max);

% Evaluate points on line
N_line = 100000;

V = [-2 -1];
W = [3 1];

lambda = linspace(0, 1, N_line);
L = V + lambda' * (W-V);

tic
u_L = evaluate_pointwise(u_sol, L(:, 1), L(:, 2), c, el, el_buckets, bucket_size, x_min, x_max, y_min, y_max);
toc
disp('point eval')

beta = (W-V)/norm(W-V);

tic
[Gx_L, Gy_L] = evaluate_gradient_pointwise(u_sol, L(:, 1), L(:, 2), c, el, el_buckets, bucket_size, x_min, x_max, y_min, y_max);
toc
disp('grad eval')

beta_times_G_L = beta * [Gx_L'; Gy_L'];


figure(8)
plot(lambda, u_L);


figure(9)
plot(lambda, beta_times_G_L);


% Evaluate points on Omega
N_grad = 1000;

X_lin = linspace(x_min, x_max, N_grad);
Y_lin = linspace(y_min, y_max, N_grad);
[X, Y] = meshgrid(X_lin, Y_lin);

X_arr = reshape(X, [N_grad*N_grad, 1]);
Y_arr = reshape(Y, [N_grad*N_grad, 1]);

Z = evaluate_pointwise(u_sol, X_arr, Y_arr, c, el, el_buckets, bucket_size, x_min, x_max, y_min, y_max);
[Gx, Gy] = evaluate_gradient_pointwise(u_sol, X_arr, Y_arr, c, el, el_buckets, bucket_size, x_min, x_max, y_min, y_max);

Z = reshape(Z, size(X));
Gx = reshape(Gx, size(X));
Gy = reshape(Gy, size(Y));


figure(5)
surf(X, Y, Z, 'EdgeColor','none')
title(sprintf('Pointwise Evaluation of size N=%d', N_grad))


figure(6)
contourf(X, Y, Z)
hold on
quiver(X, Y, Gx, Gy)
plot(L(:, 1), L(:, 2), 'b')
title('Pointwise Evaluation of the gradient')


[Gx_sol, Gy_sol] = evaluateGradient(u_sol_model, X_arr, Y_arr);
Gx_sol = reshape(Gx_sol, size(X));
Gy_sol = reshape(Gy_sol, size(Y));


figure(7)
quiver(X, Y, Gx_sol, Gy_sol);
title('Gradient evaluation via evaluateGradient')

