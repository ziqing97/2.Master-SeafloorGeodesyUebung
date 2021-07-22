%% Seafloor Geodesy Lab 4-5
% Ziqing Yu 3218051
% 14/07/2021

%% Initial
clc
close all
clear

%%  load data
load('GSSM_a0a30-2.mat')

%% refomulation
gssm = regenerator('GSSM_a0a30-2.mat');

%% Calibration
cali_p1 = gssm.p1_cali;
cali_p2 = gssm.p2_cali;
cali_time = gssm.time_cali;

[l1,l2] = size(cali_time);
p1_offset = zeros(l1,l2);
p2_offset = zeros(l1,l2);
temp_cali = zeros(l1,1);
idx1 = zeros(l1,1);
for i = 1:l1
    date_cali = datetime(year(cali_time(i,1)), month(cali_time(i,1)), day(cali_time(i,1)));
    idx1(i) = find(datenum(date_cali)==gssm.date); % which day
    idx2 = find(abs(gssm.time{idx1(i)}-datenum(cali_time(i,150))) == ...
            min(abs(gssm.time{idx1(i)}-datenum(cali_time(i,150))))); % find a time point
    baro_cali = gssm.p_baro{idx1(i)}(idx2); % baro pressure at this moment
    p1_offset(i,:) = cali_p1(i,:)-baro_cali;
    p2_offset(i,:) = cali_p2(i,:)-baro_cali;
    temp_cali(i) = gssm.t_baro{idx1(i)}(idx2);
    
    % example plot
    if i == 60
        figure
        set(gcf,'outerposition',get(0,'screensize'));
        plot(cali_p1(i,:),'LineWidth',2)
        hold on 
        plot(cali_p2(i,:),'LineWidth',2)
        plot(ones(length(cali_p1(i,:)),1)*baro_cali,'LineWidth',2)
        xlabel('time [s]')
        ylabel('pressure [hpa]')
        legend('P_{1,int}','P_{2,int}','Barometer')
        set(gca,'fontsize',20)
        title('60th A-0-A calibration')
        pbaspect([3,1,1])
        saveas(gca,'calibrationExample.png')
    end
end

% remove first and last 30 s
p1_offset_remove = p1_offset(:,31:end-30);
p2_offset_remove = p2_offset(:,31:end-30);

% mean offset
p1_offset_mean = mean(p1_offset_remove,2);
p2_offset_mean = mean(p2_offset_remove,2);
figure
set(gcf,'outerposition',get(0,'screensize'));
hold on
plot(datenum(cali_time(:,1)),p1_offset_mean,'bo')
plot(datenum(cali_time(:,1)),p2_offset_mean,'ro')
xlabel('time')
ylabel('pressure offset [hpa]')
legend('P_{1}-baro','P_{2}-baro')
datetick('x')
pbaspect([3,1,1])
set(gca,'fontsize',20)
saveas(gca,'a0a.png')

%% least square
t_cali = datenum(cali_time(:,150));
t_test = gssm.date;

temp_ref = mean(a0a30.barom_temp,'omitnan');
% p1
threshold = 1e10;
for i = 1:length(t_test)
    A = [exp(-t_cali/t_test(i)), t_cali, ones(length(t_cali),1), (temp_cali - temp_ref)];
    % A = [exp(-t_cali/t_test(i)), t_cali, ones(length(t_cali),1)];
    para = (A'*A)\A'*p1_offset_mean;
    p1_new = A * para;
    e = p1_offset_mean - p1_new;
    sigma1 = e' * e / (length(t_cali)-1);
    if sigma1<threshold
        threshold = sigma1;
        p1_offset_after = p1_new;
        para1 = para;
        d1 = i;
    end
end

% p2
threshold = 1e10;
for i = 1:length(t_test)
    A = [exp(-t_cali/t_test(i)), t_cali, ones(length(t_cali),1), (temp_cali - temp_ref)];
    % A = [exp(-t_cali/t_test(i)), t_cali, ones(length(t_cali),1)]; 
    para = (A'*A)\A'*p2_offset_mean;
    p2_new = A * para;
    e = p1_offset_mean - p2_new;
    sigma2 = e' * e / (length(t_cali)-1); 
    if sigma2<threshold
        threshold = sigma2;
        p2_offset_after = p2_new;
        d2 = i;
        para2 = para;
    end
end
figure
set(gcf,'outerposition',get(0,'screensize'));
hold on
plot(datenum(cali_time(:,1)),p1_offset_after,'bo')
plot(datenum(cali_time(:,1)),p2_offset_after,'ro')
xlabel('time')
ylabel('pressure offset [hpa]')
legend('P_{1}-baro','P_{2}-baro')
datetick('x')
pbaspect([3,1,1])
set(gca,'fontsize',20)
saveas(gca,'smoothing.png')

% 
% A_P1 = [exp(-a0a30.datenum'/t_test(d1)), a0a30.datenum', ones(length(a0a30.datenum),1)];
A_P1 = [exp(-a0a30.datenum'/t_test(d1)), a0a30.datenum', ones(length(a0a30.datenum),1), a0a30.temp1'-temp_ref];
p1_offset_all = A_P1 * para1;
%
% A_P2 = [exp(-a0a30.datenum'/t_test(d2)), a0a30.datenum', ones(length(a0a30.datenum),1)];
A_P2 = [exp(-a0a30.datenum'/t_test(d2)), a0a30.datenum', ones(length(a0a30.datenum),1), a0a30.temp2'-temp_ref];
p2_offset_all = A_P2 * para2;

% correction
p1_correct = a0a30.pressure1' - p1_offset_all;
p2_correct = a0a30.pressure2' - p2_offset_all;
days = length(p1_correct)/(24*60);
p1 = zeros(days,1);
p2 = zeros(days,1);
j = 1;
k = 1;
outlier1 = zeros(2,1);
outlier2 = zeros(2,1);
for i = 1:days
    pressure1 = p1_correct(24*60*(i-1)+1:24*60*i);
    pressure2 = p2_correct(24*60*(i-1)+1:24*60*i);
    pressure1 = outlier_filter(pressure1);
    pressure2 = outlier_filter(pressure2);
    if length(pressure1) < 1350     % 1300 not good ,1400 too much out.
        p1(i) = NaN;
        outlier1(j) = i; 
        j = j +1;
    else
        p1(i) = mean(pressure1,'omitnan');
    end
    
    if length(pressure2) < 1350     % 1300 not good ,1400 too much out.
        p2(i) = NaN;
        outlier2(j) = i; 
        k = k +1;
    else
        p2(i) = mean(pressure2,'omitnan');
    end
    
    if i == 60
        figure
        t_exp = gssm.time{i};
        set(gcf,'outerposition',get(0,'screensize'));
        hold on
        plot(t_exp,pressure1)
        plot(t_exp,pressure2)
        xlabel('time')
        ylabel('pressure [hpa]')
        legend('P_{1}','P_{2}')
        datetick('x')
        title(datestr(datetime(year(t_exp(1)),month(t_exp(1)),day(t_exp(1)))))
        pbaspect([3,1,1])
        set(gca,'fontsize',20)
        saveas(gca,'examplePressure.png')
    end
end
%% Outlier Filter & simple interpolant
% p1
outlier_idx = isoutlier(p1);
p1(outlier_idx) = NaN;

nan_idx = isnan(p1);
F = griddedInterpolant(t_test(~nan_idx),p1(~nan_idx));
p1 = F(t_test); 

% p2
outlier_idx = isoutlier(p2);
p2(outlier_idx) = NaN;

nan_idx = isnan(p2);
F = griddedInterpolant(t_test(~nan_idx),p2(~nan_idx));
p2 = F(t_test); 

%% fitting
poly1 = polyfit(t_test,p1,1);
y1 = polyval(poly1, t_test);
poly2 = polyfit(t_test,p2,1);
y2 = polyval(poly2, t_test);

figure
set(gcf,'outerposition',get(0,'screensize'));
hold on
plot(t_test,p1-90012,'b')
plot(t_test,p2-90012,'r')
plot(t_test,y1-90012,'b','LineWidth',2)
plot(t_test,y2-90012,'r','LineWidth',2)
legend('P_{1}','P_{2}','Trend P_{1}','Trend P_{2}')
xlabel('time')
ylabel('pressure - 90012 [hpa]')
datetick('x')
pbaspect([3,1,1])
set(gca,'fontsize',20)
saveas(gca,'pressure.png')

trend1 = poly1(1) * 365; % hpa/year
trend2 = poly2(1) * 365; % hpa/year

%% Accelerator
aX = zeros(length(t_test),1);
aY = zeros(length(t_test),1);
aZ = zeros(length(t_test),1);
for i = 1:length(t_test)
    aX(i) = mean(gssm.ac_x{i},'omitnan'); % vertical
    aY(i) = mean(gssm.ac_y{i},'omitnan'); % horizontal
    aZ(i) = mean(gssm.ac_z{i},'omitnan'); % horizontal
end
figure
set(gcf,'outerposition',get(0,'screensize'));
subplot(3,1,1)
plot(t_test,aX,'LineWidth',2)
datetick('x')
xlabel('time')
ylabel('acceleration [m/s^2]')
set(gca,'fontsize',20)

subplot(3,1,2)
plot(t_test,aY,'LineWidth',2)
datetick('x')
xlabel('time')
ylabel('acceleration [m/s^2]')
set(gca,'fontsize',20)
subplot(3,1,3)
plot(t_test,aZ,'LineWidth',2)
datetick('x')
xlabel('time')
ylabel('acceleration [m/s^2]')
set(gca,'fontsize',20)
saveas(gca,'acceleration.png')

% tilt
tilt = sqrt(aY.^2 + aZ.^2)./aX * 1000; % rad to mrad
figure
set(gcf,'outerposition',get(0,'screensize'));
plot(t_test,tilt,'LineWidth',2)
datetick('x')
set(gca,'fontsize',20)
ylabel('tilt [mrad]')
xlabel('time')
title("tilt")
pbaspect([3 1 1])
saveas(gca,'tilt.png')

% 
figure
plot(t_test,p1-p2,'b')


figure
set(gcf,'outerposition',get(0,'screensize'));
plot(t_cali,temp_cali,'ko')
pbaspect([3 1 1])
title('tempature')
xlabel('time')
datetick('x')
ylabel('C^{\circ}')
set(gca,'fontsize',20)
saveas(gca,'temp.png')