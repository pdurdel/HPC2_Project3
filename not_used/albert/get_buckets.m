function el_buckets = get_buckets(p, t, bucket_size, x_min, x_max, y_min, y_max)

num_x_buckets = 1 + ceil((x_max - x_min)/bucket_size);
num_y_buckets = 1 + ceil((y_max - y_min)/bucket_size);

el_buckets = cell(num_x_buckets, num_y_buckets);
el_buckets_idx = ones(num_x_buckets, num_y_buckets);
for i=1:num_x_buckets
    for j=1:num_y_buckets
        el_buckets{i, j} = zeros(80, 1);
    end
end

num_p = size(p, 1);
num_el = size(t, 1);

% get node_to_buckets
node_to_bucket = zeros(num_p, 2);
for i=1:num_p
    node_to_bucket(i, 1) = 1 + floor((p(i, 1) - x_min)/bucket_size);
    node_to_bucket(i, 2) = 1 + floor((p(i, 2) - y_min)/bucket_size);
end

% insert elements into buckets
for i=1:num_el
    n_1 = t(i, 1);
    n_2 = t(i, 2);
    n_3 = t(i, 3);
    
    n_1_x = node_to_bucket(n_1, 1);
    n_1_y = node_to_bucket(n_1, 2);
    
    n_2_x = node_to_bucket(n_2, 1);
    n_2_y = node_to_bucket(n_2, 2);

    n_3_x = node_to_bucket(n_3, 1);
    n_3_y = node_to_bucket(n_3, 2);

    nx_min = min([n_1_x, n_2_x, n_3_x]);
    nx_max = max([n_1_x, n_2_x, n_3_x]);

    ny_min = min([n_1_y, n_2_y, n_3_y]);
    ny_max = max([n_1_y, n_2_y, n_3_y]);
    
    for nx=nx_min:nx_max
        for ny=ny_min:ny_max
            el_buckets{nx, ny}(el_buckets_idx(nx, ny)) = i;
            el_buckets_idx(nx, ny) = el_buckets_idx(nx, ny) + 1;
        end
    end

end


for i=1:num_x_buckets
    for j=1:num_y_buckets
        el_buckets{i, j} = el_buckets{i, j}(1:el_buckets_idx(i, j)-1);
    end
end


end
