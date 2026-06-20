function [g_x, g_y] = evaluate_gradient_pointwise(u, x, y, p, tr, el_buckets, bucket_size, x_min, x_max, y_min, y_max)

N = size(x, 1);

g_x = NaN(N, 1);
g_y = NaN(N, 1);

for c=1:N

    if x(c) < x_min || x(c) > x_max || y(c) < y_min || y(c) > y_max
        continue
    end

    bucket_x = 1 + floor((x(c) - x_min)/bucket_size);
    bucket_y = 1 + floor((y(c) - y_min)/bucket_size);
        
    bucket = el_buckets{bucket_x, bucket_y};

    for i=1:length(bucket)
        el = tr(bucket(i), :);
        n_1_id = el(1);
        n_2_id = el(2);
        n_3_id = el(3);

        n_1 = p(n_1_id, :);
        n_2 = p(n_2_id, :);
        n_3 = p(n_3_id, :);
    
        [is_inside, ~] = get_barycentric(n_1, n_2, n_3, [x(c) y(c)]);
   
        if ~is_inside
            continue
        end

        grad_shapefct = [1 1 1; n_1(1) n_2(1) n_3(1); n_1(2) n_2(2) n_3(2)]\[0 0; 1 0; 0 1];

        grads = u([n_1_id n_2_id n_3_id])' * grad_shapefct;
        g_x(c) = grads(1);
        g_y(c) = grads(2);

        break
    end

end