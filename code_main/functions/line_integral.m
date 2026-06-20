function I = line_integral(s, vals, L)

ds = diff(s) * L;
seg_avg = (vals(1:end-1) + vals(2:end)) / 2;
valid = ~isnan(seg_avg);
I = sum(ds(valid) .* seg_avg(valid));

end
