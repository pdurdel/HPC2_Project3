function [line_integral, line_data] = eval_line_integral(model, u_sol, V, W, n_line)

    if nargin < 5
        n_line = 200;
    end

    % length and direction of the line
    line_length = norm(W - V);
    beta = (W - V) / line_length;

    % creating points on VW
    s = linspace(0, 1, n_line)';
    line_points = (1 - s).*V + s.*W;

    % getting mesh from PDE model
    mesh_points = model.Mesh.Nodes.';
    mesh_elements = model.Mesh.Elements(1:3,:).';

    % eval u on the line
    u_line = artificial_point_eval(mesh_points, mesh_elements, u_sol, line_points);
    u_line = u_line(:);

    % eval gradient on the line
    grad_u_line = artificial_point_grad(mesh_points, mesh_elements, u_sol, line_points);

    % directional derivative beta * grad(u)
    beta_grad_u = grad_u_line * beta.';

    % valid points inside the geometry
    valid = isfinite(u_line);
    valid_idx = find(valid);

    % line integral
    line_integral = 0;

    % split valid indices into connected blocks
    block_start = [1; find(diff(valid_idx) > 1) + 1];
    block_end   = [block_start(2:end) - 1; length(valid_idx)];

    for k = 1:length(block_start)
        idx = valid_idx(block_start(k):block_end(k));

        if length(idx) >= 2
            line_integral = line_integral + ...
                line_length * trapz(s(idx), u_line(idx));
        end
    end

    % saving the data
    line_data.s = s;
    line_data.line_points = line_points;
    line_data.u_line = u_line;
    line_data.grad_u_line = grad_u_line;
    line_data.beta_grad_u = beta_grad_u;
    line_data.valid = valid;
    line_data.beta = beta;
    line_data.line_length = line_length;
end