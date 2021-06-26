load('a0a30_oneday.mat')
len = 905;

p1_gauge = zeros(len,1);
p2_gauge = zeros(len,1);
p_baro = zeros(len,1);

for i = 1:len
    p1_gauge(i) = mean(a0a30_oneday{i}.pressure1,'omitnan');
    p2_gauge(i) = mean(a0a30_oneday{i}.pressure2,'omitnan');
    p_baro(i) = mean(a0a30_oneday{i}.barom_pressure,'omitnan');
end

figure
hold on
plot(p1_gauge)
plot(p2_gauge)
plot(p_baro)

figure
plot(p1_gauge - p_baro);
hold on
plot(p2_gauge - p_baro);