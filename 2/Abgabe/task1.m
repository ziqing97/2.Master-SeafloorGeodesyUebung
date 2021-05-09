clc
close all
addpath(genpath(pwd))
load('HRT2305.mat')
load('PRS2305.mat')
load('time2305.mat')
t1 = t(1);
time = datetime(t,'InputFormat','yyyy-MM-dd''T''HH:mm');
time = datenum(time);


p1 = polyfit(time,HRT(:,2),1); 
y1 = polyval(p1,time);
f1 = figure;
hold on 
plot(time, HRT(:,2))
plot(time, y1,'LineWidth',2)
title("Temperature")
ylabel("degree")
datetick('x','yyyy-mm')
set(gcf,'unit','normalized','position',[0.2,0.2,0.8,0.6]);
saveas(f1, 'TempChange.png')

p2 = polyfit(time,PRS(:,2),1); 
y2 = polyval(p2,time);
f2 = figure;
hold on 
plot(time,PRS(:,2))
plot(time, y2,'LineWidth',2)
title("Pressure")
ylabel("hpa")
datetick('x','yyyy-mm')
set(gcf,'unit','normalized','position',[0.2,0.2,0.8,0.6]);
saveas(f2, 'PressChange.png')



