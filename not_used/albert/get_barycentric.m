function [is_inside, b] = get_barycentric(p1, p2, p3, x)

is_inside = false;
b = NaN(3, 1);
area_123 = get_area_of_triangle(p1, p2, p3);

area_x23 = get_area_of_triangle(x, p2, p3);
b(1) = area_x23 / area_123;
if b(1) < 0
    return
end

area_1x3 = get_area_of_triangle(p1, x, p3);
b(2) = area_1x3 / area_123;
if b(2) < 0
    return
end

area_12x = get_area_of_triangle(p1, p2, x);
b(3) = area_12x / area_123;
if b(3) < 0
    return
end

is_inside = true;

end

