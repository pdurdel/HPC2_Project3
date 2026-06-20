function vals = artificial_point_eval(points, elements, u, line_points)
    % func vor artificial point evaluation of u on a arbitrary point

    % triangulationobject gets created
    obj_tr = triangulation(elements, points);

    % gets barycentric coords for every point
    [elem_id, bary] = pointLocation(obj_tr, line_points);

    % solution vector
    vals = NaN(size(line_points,1),1);

    for i = 1:size(line_points,1)

        % checks if point is valid
        if ~isnan(elem_id(i))

            % getting the knotindices of the corresponding traiangle of the
            % point
            tri = elements(elem_id(i),:);

            % getting the u-vals of the corresponding vertices of the
            % triangle
            u_local = u(tri);

            % calcualtes u with the baryc. coords
            vals(i) = bary(i,:) * u_local;

        end
    end
end