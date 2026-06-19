function rand_points = artificial_point_rand(points, elements, n_points)
    % func for random points on the grid

    % coords of the knot vertices
    p1 = points(elements(:,1),:);
    p2 = points(elements(:,2),:);
    p3 = points(elements(:,3),:);

    % area of the triangles
    areas = 0.5 * abs( ...
        (p2(:,1)-p1(:,1)).*(p3(:,2)-p1(:,2)) - ...
        (p3(:,1)-p1(:,1)).*(p2(:,2)-p1(:,2)));

    % cummulated sum
    cum_areas = cumsum(areas);
    total_area = cum_areas(end);

    % choosing a random trangle
    r = rand(n_points, 1) * total_area;
    tri_id = discretize(r, [0; cum_areas]);

    % getting the vertices of the corresponding triangle
    A = p1(tri_id,:);
    B = p2(tri_id,:);
    C = p3(tri_id,:);

    % choosing a random point in the triangle
    r1 = sqrt(rand(n_points,1));
    r2 = rand(n_points,1);
    rand_points = (1-r1).*A + r1.*(1-r2).*B + r1.*r2.*C;

end