clc
close all
load('loihi_oceandata.mat')
sound_speeds = double(mean(loihi.sound_speed,1,'omitnan'));
% sound_speeds = loihi.sound_speed(1,:);
sound_speeds = sound_speeds';

layer_depths = double(loihi.depth);

for i = 1:5
    figure
    transponder_depth = i*200;
    transponder_range = 600;
    [ranges,depths,sound_speed,sound_speed_gradient] = ray_trace_test(sound_speeds,layer_depths,transponder_depth,transponder_range);
    plot(ranges, depths)
    xlabel('ranges')
    ylabel('depths')
    set(gca,'YDir','reverse')
end


% 
% transponder_depth = 0:1:1400;
% for i = 1:1401
%     [ranges,depths,sound_speed(i),sound_speed_gradient(i)] = ray_trace_test(sound_speeds,layer_depths,transponder_depth(i),transponder_range);
% end
% 
% figure
% plot(transponder_depth, sound_speed);
% title('sound speed')
% xlabel('depth m')
% ylabel('sound speed m/s')
% 
% figure
% plot(transponder_depth, sound_speed_gradient);
% title('sound speed gradient')
% xlabel('depth m')
% ylabel('sound speed gradient ')
