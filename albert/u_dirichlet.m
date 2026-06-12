function y = u_dirichlet(x, M, r, H, G, tol)

y = zeros(1, size(x, 1));

for i=1:size(x, 1)


dist_to_M = abs((M(1) - x(i, 1))^2 + (M(2) - x(i, 2))^2 - r^2);

if dist_to_M < tol
    y(i) = 10;
    continue
end


HG_cross_x = abs((G(1) - H(1)) * (x(i, 2) - H(2)) - (x(i, 1) - H(1)) * (G(2) - H(2)));
min_x1 = min(H(1), G(1));
max_x1 = max(H(1), G(1));
min_x2 = min(H(2), G(2));
max_x2 = max(H(2), G(2));

if HG_cross_x < tol
    
    if min_x1 <= x(i, 1) && x(i, 1) <= max_x1 && min_x2 <= x(i, 2) && x(i, 2) <= max_x2
        y(i) = 10;
        continue
    end

end



y(i)=20;


end



end