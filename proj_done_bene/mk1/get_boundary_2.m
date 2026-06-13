function [bdry_dirichlet, bdry_neumann] = get_boundary_2(e, bdry_dirichlet_idx, bdry_neumann_idx)

if nargin < 3
    skip = true;
else
    skip = false;
end

idx1 = find( e(:, 7) == 0 );
idx2 = find( e(:, 6) == 0 );
idx = [idx1; idx2];

bdry_dirichlet = zeros(length(idx), 3);
bdry_neumann = zeros(length(idx), 3);

idx_dirichlet = 1;
idx_neumann = 1;

for i=1:length(idx)

    edge_id = idx(i);
    e_edge = e(edge_id, :);
    e_edge_id = e_edge(5);
    edge_id = e_edge(5);
    found = false;

    for j=1:length(bdry_dirichlet_idx)
        if bdry_dirichlet_idx(j) == e_edge_id
            
            if e_edge(6) == 1
                bdry_dirichlet(idx_dirichlet, :) = [edge_id e_edge(1) e_edge(2)];
            else
                bdry_dirichlet(idx_dirichlet, :) = [edge_id e_edge(2) e_edge(1)];
            end
            idx_dirichlet = idx_dirichlet + 1;
            
            found = true;
            break

        end
    end

    if ~skip && ~found
        for j=1:length(bdry_neumann_idx)
            if bdry_neumann_idx(j) == e_edge_id

                if e_edge(6) == 1
                    bdry_neumann(idx_neumann, :) = [edge_id e_edge(1) e_edge(2)];
                else
                    bdry_neumann(idx_neumann, :) = [edge_id e_edge(2) e_edge(1)];
                end
                idx_neumann = idx_neumann + 1;
            
                break
            end
        end
    %{
    elseif skip

        if e_edge(6) == 1
            bdry_neumann(idx_neumann, :) = [edge_id e_edge(1) e_edge(2)];
        else
            bdry_neumann(idx_neumann, :) = [edge_id e_edge(2) e_edge(1)];
        end
        idx_neumann = idx_neumann + 1;
    %}
    end

end

bdry_dirichlet = bdry_dirichlet(1:idx_dirichlet-1, :);
bdry_neumann = bdry_neumann(1:idx_neumann-1, :);

end