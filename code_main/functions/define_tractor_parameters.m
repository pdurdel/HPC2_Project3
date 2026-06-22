function [coords, coords_circles, big_radii, small_radii] = define_tractor_parameters()

% This function returns parameters for a tractor shape. The parameters are
% passed to get_geometry.m along hmax to initialize the mesh.

%% point coordinate definitions

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

%% wheel parameters

coords_circles = [
    -6.5 2.25
    0 0
];

big_radii = [3.75 1.5];
small_radii = [1 0.5];

end