clc
close all
addpath(genpath(pwd))
load('HRT2305.mat')
load('PRS2305.mat')
load('time2305.mat')
t1 = t(1);
a = datetime(t,'InputFormat','yyyy-MM-dd''T''HH:mm');


figure
plot(a,HRT(:,2))
title("Temperature")
% xlabel("from 2014 to 2017")
ylabel("degree")


figure
plot(a,PRS(:,2))
title("Pressure")
ylabel("hpa")
% xlabel("from 2014 to 2017")



