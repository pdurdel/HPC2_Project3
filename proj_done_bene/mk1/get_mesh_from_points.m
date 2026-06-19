function [p,e,t,geom] = get_mesh_from_points(coords, coords_circles, big_radii, small_radii, hmax)
%% Check for Parameters and allocate Memory for the geometry
num_coords = 13;
if size(coords, 1) ~= num_coords
    error('There are 13 point coordinates required.')
end

num_circles = 2;
if size(coords_circles, 1) ~= num_circles
    error('There are 2 circle coordinates required.')
end

if big_radii(1) <= small_radii(1)
    error('Radius m is too large')
end

if big_radii(2) <= small_radii(2)
    error('Radius n is too large')
end

geom = zeros(10, num_coords - 2 + 3 * num_circles);
geom_idx = 1;

%% Fill the geometry matrix

% Fill the straight lines from D' to L
for i=4:num_coords-1
    geom(:,geom_idx) = [
        2
        coords(i, 1)
        coords(i+1, 1)
        coords(i, 2)
        coords(i+1, 2)
        1
        0
        0
        0
        0
    ];
    geom_idx = geom_idx + 1;
end

% Fill the line BC and LA
geom(:, geom_idx:geom_idx+1) = [
    2               2
    coords(2, 1)    coords(end, 1)
    coords(3, 1)    coords(1, 1)
    coords(2, 2)    coords(end, 2)
    coords(3, 2)    coords(1, 2)
    1               1
    0               0
    0               0
    0               0
    0               0
];
geom_idx = geom_idx + 2;

% Fill the circles connecting AB and CD'
geom(:, geom_idx:geom_idx+1) = [
    1                       1
    coords(1, 1)            coords(3, 1)
    coords(2, 1)            coords(4, 1)
    coords(1, 2)            coords(3, 2)
    coords(2, 2)            coords(4, 2)
    1                       1
    0                       0
    coords_circles(1, 1)    coords_circles(2, 1)
    coords_circles(1, 2)    coords_circles(2, 2)
    big_radii(1)            big_radii(2)
];
geom_idx = geom_idx + 2;

% Remove the smaller circles
for i=1:num_circles
    M = coords_circles(i, :);
    r = small_radii(i);
    geom(:, geom_idx:geom_idx+1) = [
        1       1
        M(1)+r  M(1)-r
        M(1)-r  M(1)+r
        M(2)    M(2)
        M(2)    M(2)
        0       0
        1       1
        M(1)    M(1)
        M(2)    M(2)
        r       r
    ];
    geom_idx = geom_idx + 2;
end

%% Compute the mesh to the given maximal mesh size

[p,e,t] = initmesh(geom, 'hmax', hmax);

end