function y = u_dirichlet(x, M, r, H, G, tol)

y = zeros(1, size(x, 1));

for i=1:size(x, 1)

    % Check distance from the circle around M to x(i)
    dist_to_circle_around_M = abs((M(1) - x(i, 1))^2 + (M(2) - x(i, 2))^2 - r^2);

    if dist_to_circle_around_M < tol
        y(i) = 10;
        continue
    end

    min_x1 = min(H(1), G(1));
    max_x1 = max(H(1), G(1));
    min_x2 = min(H(2), G(2));
    max_x2 = max(H(2), G(2));
        
    % Check whether x(i) is in the box around H and G within tolerance
    if min_x1-tol <= x(i, 1) && x(i, 1) <= max_x1+tol && min_x2-tol <= x(i, 2) && x(i, 2) <= max_x2+tol
        % Check the area of spanned parallelogram of H, G and x(i)
        area_parallelogram_HGx = det([1 1 1; H(1) G(1) x(i, 1); H(2) G(2) x(i, 2)]);
    
        if abs(area_parallelogram_HGx) < tol
            y(i) = 10;
            continue
        end

    end

    y(i)=20;
end

end