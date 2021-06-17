function[t1,y1] = nan_filter(t,y)
idx_nan = isnan(y);
y1 = y(~idx_nan);
t1 = t(~idx_nan);
end