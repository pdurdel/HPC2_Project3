function bdry = get_boundary_3(e, bdry_idcs, num_boundaries)
%% Set up the output

% Find all boundary edges
idx_bdry_edges = find( e(:, 6) == 0 | e(:, 7) == 0 );

bdry_len = length(bdry_idcs);

bdry = cell(1, bdry_len);
idcs = zeros(1, bdry_len);
for i=1:bdry_len
    bdry{i} = zeros(length(idx_bdry_edges), 3);
    idcs(i) = 1;
end

% EdgeID-to-Boundary maps the Boundary index which is stored for each edge
% to the corresponding Boundary index given in bdry_idcs; Only one Boundary
% index appearance in bdry_idcs is allowed for each Boundary index.
edge_id_to_bdry = zeros(num_boundaries, 1);
for i=1:bdry_len
    for j=1:length(bdry_idcs{i})
        edge_id_to_bdry(bdry_idcs{i}(j)) = i;
    end
end

%% Map all Boundary edges to the corresponding Boundary

for i=1:length(idx_bdry_edges)

    edge = e(idx_bdry_edges(i), :);
    edge_id = edge(5);
    
    bdry_idx = edge_id_to_bdry(edge_id);

    % Check the orientation and save the edge going in mathematical
    % positive direction
    if edge(6) == 1
        bdry{bdry_idx}(idcs(bdry_idx), :) = [edge_id edge(1) edge(2)];
    else
        bdry{bdry_idx}(idcs(bdry_idx), :) = [edge_id edge(2) edge(1)];
    end
    idcs(bdry_idx) = idcs(bdry_idx) + 1;

end

%% Adjust lengths for the Boundaries
for i=1:bdry_len
    bdry{i} = bdry{i}(1:idcs(i)-1, :);
end

end
