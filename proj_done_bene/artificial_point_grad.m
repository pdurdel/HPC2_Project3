function grad_values = artificial_point_grad(points, elements, u, line_points)
    % func for calculating the gradient of u 

    % triangulationobject gets created
    obj_tr = triangulation(elements, points);

    % gets the corresponding triangle for every point
    elem_id = pointLocation(obj_tr, line_points);

    % creates gradient matrice
    grad_values = NaN(size(line_points,1),2);

    for i = 1:size(line_points,1)
        
        % checks wether its a valid point
        if ~isnan(elem_id(i))

            % getting the knotindices of the corresponding traiangle of the
            % point
            tri = elements(elem_id(i),:);

            % getting the coords of the vertices of the triangle
            x = points(tri,1); y = points(tri,2);
            u_local = u(tri);
    
            % sys of lin.equ. gets created -> u(x,y) = a + b*x + c*y
            A = [
                1 x(1) y(1)
                1 x(2) y(2)
                1 x(3) y(3)
            ];
        
            % coefficients are calculated
            coeff = A\u_local;
            
            % saving the gradients as [du/dx, du/dy] = [b, c]
            grad_values(i,:) = [coeff(2), coeff(3)];

        end
    end
end