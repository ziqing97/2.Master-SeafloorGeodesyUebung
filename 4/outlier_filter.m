function[t1,y1] = outlier_filter(t,y)
% outlier filter
idx_outlier = isoutlier(y);
y1 = y(~idx_outlier);
t1 = t(~idx_outlier);
end