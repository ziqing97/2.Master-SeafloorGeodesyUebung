function[y1] = nan_filter(y)
% nan filter
idx_nan = isnan(y);
y1 = y(~idx_nan);
end