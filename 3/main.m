% Seafloor Geodesy Lab 3
% Ziqing Yu
% 3218051

%% Initial 
clc
close all

%% Load data
addpath(genpath(pwd));
gnss_data = importGNSS('gnss.txt', 5, 296283);
speed = importSpeed('sound speed.txt',2,35);
Coor_transponder = [40.745979, 40.756364, 40.8047, 40.718027;
                144.074123, 144.038315, 144.020393, 144.055476;
                -4669.618,-4696.425,-4617.816,-4635.294];
% time
t1 = gnss_data(1:1000,1);
t2 = gnss_data(1:1000,2);

% lontitude & latitude & elevation
lon = gnss_data{:,3};
lat = gnss_data{:,4};
ele = gnss_data{:,5};

% transducer coordinate calculation using 'Transducer_position'
[lat_transducer, lon_transducer, elev_transducer] = Transducer_position(lat,lon,ele,t2);

% plot transducer(black) & transponder(red) & constline(blue)
figure
worldmap([39,46],[139 148])
load coastlines
plotm(coastlat,coastlon)
hold on
linem(lat_transducer,lon_transducer,'k-')
geoshow(Coor_transponder(1,:),Coor_transponder(2,:),'Marker','o','MarkerEdgeColor','red')

%% Task 3
function [t_oneway,r_oneway,x_range,layers] = ray_trace_chadwell(speed(:,2), speed(:,1), )
