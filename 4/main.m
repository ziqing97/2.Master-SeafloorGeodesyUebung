clc
close all
clearvars
load('GSSM_a0a30.mat');

t = a0a30.datenum;
time = datetime(t,'ConvertFrom','datenum');
%%  split data into individual day
day_count = length(t) / 24 / 60;
a0a30_oneday = cell(day_count,1);
for i =1:day_count
    a0a30_oneday{i} = struct;
    a0a30_oneday{i}.datenum = a0a30.datenum((i-1)*24*60 + 1: i*24*60);
    a0a30_oneday{i}.datetime = datetime(a0a30_oneday{i}.datenum,'ConvertFrom','datenum');
    a0a30_oneday{i}.pressure1 = a0a30.pressure1((i-1)*24*60 + 1: i*24*60);
    a0a30_oneday{i}.pressure2 = a0a30.pressure2((i-1)*24*60 + 1: i*24*60);
    a0a30_oneday{i}.temp1 = a0a30.temp1((i-1)*24*60 + 1: i*24*60);
    a0a30_oneday{i}.temp2 = a0a30.temp2((i-1)*24*60 + 1: i*24*60);
    a0a30_oneday{i}.barom_pressure = a0a30.barom_pressure((i-1)*24*60 + 1: i*24*60);
    a0a30_oneday{i}.barom_temp = a0a30.barom_temp((i-1)*24*60 + 1: i*24*60);
    a0a30_oneday{i}.accelX = a0a30.accelX((i-1)*24*60 + 1: i*24*60);
    a0a30_oneday{i}.accelY = a0a30.accelY((i-1)*24*60 + 1: i*24*60);
    a0a30_oneday{i}.accelZ = a0a30.accelZ((i-1)*24*60 + 1: i*24*60);
    a0a30_oneday{i}.accel_temp = a0a30.accel_temp((i-1)*24*60 + 1: i*24*60);
    a0a30_oneday{i}.accel_rss = a0a30.accel_rss((i-1)*24*60 + 1: i*24*60);
end

%% plot example
% 2018.10.31
figure
subplot(4,1,1)
plot(a0a30_oneday{501}.datenum, a0a30_oneday{501}.pressure1)
hold on
plot(a0a30_oneday{501}.datenum, a0a30_oneday{501}.pressure2)
title('pressure')
xlabel('time')
ylabel('pressure hpa')
legend('gauge 1', 'gauge 2')
datetick('x')

subplot(4,1,2)
plot(a0a30_oneday{501}.datenum, a0a30_oneday{501}.accelX)
hold on
plot(a0a30_oneday{501}.datenum, a0a30_oneday{501}.accelY)
plot(a0a30_oneday{501}.datenum, a0a30_oneday{501}.accelZ)
title('accelleration')
xlabel('time')
ylabel('acc m/s^2')
legend('X', 'Y', 'Z')
datetick('x')

subplot(4,1,3)
plot(a0a30_oneday{501}.datenum, a0a30_oneday{501}.barom_pressure)
title('barom pressure')
xlabel('time')
ylabel('pressure hpa')
datetick('x')

subplot(4,1,4)
plot(a0a30_oneday{501}.datenum, a0a30_oneday{501}.temp1)
hold on
plot(a0a30_oneday{501}.datenum, a0a30_oneday{501}.temp2)
title('temperature')
xlabel('time')
ylabel('celcius')
legend('1','2')
datetick('x')


p1 = a0a30.pressure1;
[t_p1,p1] = outlier_filter(t,p1);    % what if we do the interpolation? not good because the gap is too big
figure
subplot(2,1,1)
plot(t_p1,p1)
datetick("x")
title('pressure gauge 1')
xlabel('time')
ylabel('hpar')

p2 = a0a30.pressure2;
[t_p2,p2] = outlier_filter(t,p2);    
subplot(2,1,2)
plot(t_p2,p2)
datetick("x")
title('pressure gauge 2')
xlabel('time')
ylabel('pressure hpar')



