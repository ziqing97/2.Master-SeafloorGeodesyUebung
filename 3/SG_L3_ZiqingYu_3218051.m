% Seafloor Geodesy Lab 3
% Ziqing Yu
% 3218051

%% Initial 
clc
close all

% Load data
addpath(genpath(pwd));
gnss_data = importGNSS('gnss.txt', 5, 296283);
speed = importSpeed('sound speed.txt',2,35);
Coor_transponder = [40.745979, 40.756364, 40.728047, 40.718027;
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

% plot with outliers transducer(black) & transponder(red) & constline(blue)
figure
load coastlines
hold on
plot(coastlon,coastlat)
scatter(Coor_transponder(2,:),Coor_transponder(1,:),'Marker','o','MarkerEdgeColor','red')
grid on
plot(lon_transducer,lat_transducer,'k-')
xlim([140,146])
ylim([37,46])

% filter
index1 = find(lat_transducer<40.76 & lat_transducer>40.71);
lat_transducer = lat_transducer(index1);
lon_transducer = lon_transducer(index1);
elev_transducer = elev_transducer(index1);

index2 = find(lon_transducer < 144.1 & lon_transducer > 144);
lat_transducer = lat_transducer(index2);
lon_transducer = lon_transducer(index2);
elev_transducer = elev_transducer(index2);

% plot transducer(black) & transponder(red) & constline(blue)
figure
hold on
plot(lon_transducer,lat_transducer,'k-')
plot(coastlon,coastlat)
scatter(Coor_transponder(2,:),Coor_transponder(1,:),'Marker','o','MarkerEdgeColor','red')
grid on
xlim([144,144.1])
ylim([40.71,40.76])

%% Task 4
geodesic_range = 0; % somehow not used in function
receiver_depth = -Coor_transponder(3,:); % sign change
source_depth = elev_transducer;
theta_check = 35; % angel in degree

% modell expand using 1D interpolation
layer_depths = 0:100:5000;
layer_depths = layer_depths';
F = griddedInterpolant(speed{:,1},speed{:,2});
sound_speeds = F(layer_depths);

% here we go
t_oneway = zeros(length(source_depth),4);
for i=1:length(source_depth)
    for j=1:4
        [t_oneway(i,j),~,~,~] = ray_trace_chadwell(sound_speeds, layer_depths,...
                            receiver_depth(j),geodesic_range, source_depth(i), theta_check/180*pi);
    end
end

% plotting
figure
hold on 
plot(t_oneway(:,1))
plot(t_oneway(:,2))
plot(t_oneway(:,3))
plot(t_oneway(:,4))
legend("Transponder A","Transponder B","Transponder C","Transponder D")
ylabel('one way time (m/s)')




