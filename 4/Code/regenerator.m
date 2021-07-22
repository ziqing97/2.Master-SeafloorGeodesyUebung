function[gssm] = regenerator(file)
% Ziqing Yu
% 13/07/2021

% this function split the whole time series by day
data = load(file);
a0a30 = data.a0a30;
gssm = struct;
t = a0a30.datenum;
time = datetime(t,'ConvertFrom','datenum');
p1 = a0a30.pressure1;
p2 = a0a30.pressure2;
t1 = a0a30.temp1;
t2 = a0a30.temp2;
p_baro = a0a30.barom_pressure;
t_baro = a0a30.barom_temp;
ac_x = a0a30.accelX;
ac_y = a0a30.accelY;
ac_z = a0a30.accelZ;
ac_t = a0a30.accel_temp;
ac_rss = a0a30.accel_rss;

days = length(p1)/(24*60);
for i = 1:days
    gssm.p1{i} = p1(24*60*(i-1)+1:24*60*i);
    gssm.p2{i} = p2(24*60*(i-1)+1:24*60*i);
    gssm.t1{i} = t1(24*60*(i-1)+1:24*60*i);
    gssm.t2{i} = t2(24*60*(i-1)+1:24*60*i);
    gssm.time{i} = datenum(time(24*60*(i-1)+1:24*60*i));
    gssm.p_baro{i} = p_baro(24*60*(i-1)+1:24*60*i);
    gssm.t_baro{i} = t_baro(24*60*(i-1)+1:24*60*i);
    gssm.ac_x{i} = ac_x(24*60*(i-1)+1:24*60*i);
    gssm.ac_y{i} = ac_y(24*60*(i-1)+1:24*60*i);
    gssm.ac_z{i} = ac_z(24*60*(i-1)+1:24*60*i);
    gssm.ac_t{i} = ac_t(24*60*(i-1)+1:24*60*i);
    gssm.ac_rss{i} = ac_rss(24*60*(i-1)+1:24*60*i);
    gssm.year(i) = year(gssm.time{i}(1));
    gssm.month(i) = month(gssm.time{i}(1));
    gssm.day(i) = day(gssm.time{i}(1));
    gssm.date(i) = datenum(datetime(gssm.year(i),gssm.month(i),gssm.day(i)));
end

a0a30.calibration_datenum = datetime(a0a30.calibration_datenum,'ConvertFrom','datenum');
gssm.time_cali = a0a30.calibration_datenum;
gssm.p1_cali = a0a30.calibration_pressure1;
gssm.p2_cali = a0a30.calibration_pressure2;
end