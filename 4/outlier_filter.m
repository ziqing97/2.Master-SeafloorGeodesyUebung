function[y1] = outlier_filter(y)
% outlier filter
idx_outlier = isoutlier(y);
y1 = y(~idx_outlier);
end