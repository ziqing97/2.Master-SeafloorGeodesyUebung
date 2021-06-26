clc
close all
clearvars
load('GSSM_a0a30.mat');

t = a0a30.datenum;
time = datetime(t,'ConvertFrom','datenum');
p1 = a0a30.pressure1;
p2 = a0a30.pressure2;

p1 = movmean(p1,1440,'omitnan');
p2 = movmean(p2,1440,'omitnan');
%% example plotting

[t_p1,p1] = outlier_filter(t,p1);    % what if we do the interpolation? not good because the gap is too big
figure
subplot(2,1,1)
plot(t_p1,p1)
datetick("x")
title('pressure gauge 1')
xlabel('time')
ylabel('hpar')


[t_p2,p2] = outlier_filter(t,p2);    
subplot(2,1,2)
plot(t_p2,p2)
datetick("x")
title('pressure gauge 2')
xlabel('time')
ylabel('pressure hpar')



