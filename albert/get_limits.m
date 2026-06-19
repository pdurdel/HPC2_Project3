function [x_min, x_max, y_min, y_max] = get_limits(coords, coords_circles, big_radii)

x_min = min([coords(:, 1); coords_circles(:, 1)-big_radii(:, 1)]);
x_max = max([coords(:, 1); coords_circles(:, 1)+big_radii(:, 1)]);

y_min = min([coords(:, 2); coords_circles(:, 2)-big_radii(:, 1)]);
y_max = max([coords(:, 2); coords_circles(:, 2)+big_radii(:, 1)]);

end

