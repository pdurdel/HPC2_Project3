function [bdry_dirichlet_idx, bdry_neumann] = define_boundary_conditions(e)

%% dirichlet and neumann edge IDs

bdry_dirichlet_idx_u10 = [5 12];
bdry_dirichlet_idx_u20 = 14:17;
bdry_dirichlet_idx = [bdry_dirichlet_idx_u10 bdry_dirichlet_idx_u20];
bdry_neumann_idx = [1:4 6:11 13];

%% extract boundary edges

bdry_idcs = {bdry_dirichlet_idx_u10, bdry_dirichlet_idx_u20, bdry_neumann_idx};
bdry = get_boundary_3(e', bdry_idcs, 17);

bdry_neumann = bdry{3};

end
