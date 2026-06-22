function rand_points = artificial_point_rand(points, elements, n_points)
    % func for radnom points on the grid

    % coords of the knot vertices
    p1 = points(elements(:,1),:);
    p2 = points(elements(:,2),:);
    p3 = points(elements(:,3),:);

    % area of the triangles
    areas = 0.5 * abs( ...
        (p2(:,1)-p1(:,1)).*(p3(:,2)-p1(:,2)) - ...
        (p3(:,1)-p1(:,1)).*(p2(:,2)-p1(:,2)));

    % cummulated sum of the areas
    cum_areas = cumsum(areas);
    total_area = cum_areas(end);

    rand_points = zeros(n_points,2);

    for i = 1:n_points

        % choosing a random triangle
        r = rand * total_area;
        tri_id = find(cum_areas >= r, 1, 'first');

        % getting the vertices of the corresponding triangle
        A = points(elements(tri_id,1),:);
        B = points(elements(tri_id,2),:);
        C = points(elements(tri_id,3),:);

        % choosing a random point in the triangle
        r1 = sqrt(rand);
        r2 = rand;
        p = (1-r1)*A + r1*(1-r2)*B + r1*r2*C;
        rand_points(i,:) = p;

    end
end