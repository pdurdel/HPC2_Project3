function [bdry_dirichlet_idx, bdry_dirichlet, bdry_neumann_idx, bdry_neumann] = define_boundary_conditions(e)
%% Define dirichlet and neumann edge IDs

bdry_dirichlet_idx = [5 12 14:17];
bdry_neumann_idx = [1:4 6:11 13];
bdry_idcs = {bdry_dirichlet_idx, bdry_neumann_idx};

%% Extract boundary edges via get_boundary.m

bdry = get_boundary(e, bdry_idcs, 17);
bdry_dirichlet = bdry{1};
bdry_neumann = bdry{2};

end
