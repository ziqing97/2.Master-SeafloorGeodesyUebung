clc
close all

addpath(genpath(pwd));
gnss_data = importGNSS('gnss.txt', 5, 296283);

t1 = gnss_data(1:1000,1);
t2 = gnss_data(1:1000,2);

lon = gnss_data{:,3};
lat = gnss_data{:,4};
ele = gnss_data{:,5};

[lat_transducer, lon_transducer, elev_transducer] = Transducer_position(lat,lon,ele,t2);

figure
worldmap([39,46],[142 145])
hold on
geoshow(lat_transducer,lon_transducer)


