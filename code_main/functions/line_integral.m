function [I, mean_u_line] = line_integral(s, vals, L)

    ds = diff(s) * L;

    % midvalue of the fem solution of every segment
    seg_avg = (vals(1:end-1) + vals(2:end)) / 2;

    % valid segments
    valid_segments = ~isnan(seg_avg);

    % lineintegral over the valid segments
    I = sum(ds(valid_segments) .* seg_avg(valid_segments));

    % length of the valid line integral
    L_valid = sum(ds(valid_segments));

    % midvalue of the fem solution within the valid points
    mean_u_line = I / L_valid;

end