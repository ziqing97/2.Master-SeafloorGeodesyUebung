%%
clc
close all
clear
%%
load('a0a30_oneday.mat')
len = length(a0a30_oneday);
p1 = zeros(len,1);
p2 = zeros(len,1);
time = cell(len,1);
datenum = zeros(len,1);

j = 1;
k = 1;
outlier1 = zeros(2,1);
outlier2 = zeros(2,1);
temp1 = zeros(len,1);
%% average
for i = 1:len
    temp1(i) = mean(a0a30_oneday{i}.temp1,'omitnan');
    
    pressure1 = a0a30_oneday{i}.pressure1;
    pressure1 = outlier_filter(pressure1);
    pressure1 = nan_filter(pressure1);
    
    pressure2 = a0a30_oneday{i}.pressure2;
    pressure2 = outlier_filter(pressure2);
    pressure2 = nan_filter(pressure2);
    
    if length(pressure1) < 1350     % 1300 not good ,1400 too much out.
        p1(i) = NaN;
        outlier1(j) = i; 
        j = j +1;
    else
        p1(i) = mean(pressure1,'omitnan');
    end
    
    if length(pressure2) < 1350     % 1300 not good ,1400 too much out.
        p2(i) = NaN;
        outlier2(j) = i; 
        k = k +1;
    else
        p2(i) = mean(pressure2,'omitnan');
    end
    
    datenum(i) = a0a30_oneday{i}.datenum(1);
    time{i} = a0a30_oneday{i}.datetime(1);
end
figure
plot(p1)

% p1
outlier_idx = isoutlier(p1);
p1(outlier_idx) = NaN;

nan_idx = isnan(p1);
F = griddedInterpolant(datenum(~nan_idx),p1(~nan_idx));
p1 = F(datenum); % comment this line, first two are nan, which leads a wrong interpolant

% p2
outlier_idx = isoutlier(p2);
p2(outlier_idx) = NaN;

nan_idx = isnan(p2);
F = griddedInterpolant(datenum(~nan_idx),p2(~nan_idx));
p2 = F(datenum); % comment this line, first two are nan, which leads a wrong interpolant

datenum = datenum(3:end);
p1 = p1(3:end);
p2 = p2(3:end);
figure
hold on
plot(datenum,p1)
plot(datenum,p2)
datetick('x')

poly = polyfit(datenum,p1,1);
y1 = polyval(poly, datenum);
poly = polyfit(datenum,p2,1);
y2 = polyval(poly, datenum);

plot(datenum,y1);
plot(datenum,y2);
legend('p1','p2','trend p1','trend p2')
ylabel('hpa')
xlabel('time')
hold off

figure
dp = p1 - p2;
plot(datenum,dp)
hold on
poly = polyfit(datenum,dp,1);
y_dp = polyval(poly, datenum);
plot(datenum,y_dp)
datetick('x')
ylabel('hpa')
xlabel('time')
hold off