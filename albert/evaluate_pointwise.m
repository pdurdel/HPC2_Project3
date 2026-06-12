function z = evaluate_pointwise(u, x, y, p, el, el_buckets, bucket_size, x_min, x_max, y_min, y_max)

N = size(x, 1);

z = NaN(N, 1);

for c=1:N

    if x(c) < x_min || x(c) > x_max || y(c) < y_min || y(c) > y_max
        continue
    end

    bucket_x = 1 + floor((x(c) - x_min)/bucket_size);
    bucket_y = 1 + floor((y(c) - y_min)/bucket_size);

    bucket = el_buckets{bucket_x, bucket_y};

    for i=1:length(bucket)
        e = el(bucket(i), :);
        n_1_id = e(1);
        n_2_id = e(2);
        n_3_id = e(3);

        n_1 = p(n_1_id, :);
        n_2 = p(n_2_id, :);
        n_3 = p(n_3_id, :);
    
        [is_inside, b] = get_barycentric(n_1, n_2, n_3, [x(c) y(c)]);
   
        if ~is_inside
            continue
        end

        z(c) = b(1) * u(n_1_id) + b(2) * u(n_2_id) + b(3) * u(n_3_id);
        %z(c) = u([n_1_id n_2_id n_3_id])' * b;

        break
    end

end