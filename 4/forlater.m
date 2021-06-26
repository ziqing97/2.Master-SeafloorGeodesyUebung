data = a0a30_oneday{226};

p1 = data.pressure1;
p2 = data.pressure2;
t1 = data.temp1;
t2 = data.temp2;

p_b = data.barom_pressure;
t_b = data.barom_temp;

aX = data.accelX;
aY = data.accelY;
aZ = data.accelZ;

t_a = data.accel_temp;
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