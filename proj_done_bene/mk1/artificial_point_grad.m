function grad_values = artificial_point_grad(points, elements, u, line_points)
    % func for calculating the gradient of u 

    % triangulationobject gets created
    obj_tr = triangulation(elements, points);

    % gets the corresponding triangle for every point
    elem_id = pointLocation(obj_tr, line_points);

    % creates gradient matrice
    grad_values = NaN(size(line_points,1),2);

    valid = ~isnan(elem_id);   

    num_elem = size(elements,1);
    elem_grad = zeros(num_elem,2);

    for k = 1:num_elem

        % getting the knotindices of the corresponding traiangle of the
        % point
        tri = elements(k,:);
        
        % getting the coords of the vertices of the triangle
        x = points(tri,1); 
        y = points(tri,2);
        u_local = u(tri);
        
        % sys of lin.equ. gets created -> u(x,y) = a + b*x + c*y
        A = [ones(3,1), x, y];
        
        % coefficients are calculated
        coeff = A\u_local;

        elem_grad(k,:) = coeff(2:3).';

    end
        
    % saving the gradients as [du/dx, du/dy] = [b, c]
    grad_values(valid,:) = elem_grad(elem_id(valid),:);

end