function [ranges,depths,sound_speed,sound_speed_gradient] = ray_trace_test(sound_speeds,layer_depths,transponder_depth,transponder_range)
% ray_trace_test takes a vertical sound speed profile with the depths of
% the speeds and uses Snell's law for calculate the ray path between two
% transponders separated by transponder_range at the same transponder_depth
% It returns a vector of ranges and depths for plotting and the sound_speed
% and sound_speed_gradient at the transponder depth

c0 = interp1(layer_depths,sound_speeds,transponder_depth);
ss_grad = diff(interp1(layer_depths,sound_speeds,transponder_depth+[-.5 .5]))

tantheta0 = transponder_range*ss_grad./(2.*c0);
theta0 = atan(tantheta0);
launch_angle = theta0.*180./pi;
a = cos(theta0)./c0;    % snells law: this value is conserved (cos(theta0)./c0 = cos(thetai)./ci);
cturn = 1./a;   % the sound speed at which the ray-path becomes horizontal

dz = (cturn-c0)./ss_grad;   % the depth change to reach c_turn
zturn = transponder_depth+dz; % the absolute depth for c_cturn
dr = 2.*(sqrt(1-a.^2.*c0.^2))./(a.*ss_grad); % check that the range to z_turn is half the full transponder range

% now turn that into something we can plot:
thetai = linspace(theta0,0,101);
ci = cos(thetai)./a;
dzi = (ci-c0)./ss_grad;
if ss_grad > 0
    dri = (sqrt(1-a.^2.*c0.^2)-sqrt(1-a.^2.*ci.^2))./(a.*ss_grad);
else
    dri = (sqrt(1-a.^2.*ci.^2)-sqrt(1-a.^2.*c0.^2))./(a.*ss_grad);
end

zi = [dzi fliplr(dzi(1:end-1))];
ri = [dri 2.*dri(end)-fliplr(dri(1:end-1))];

ranges = ri;
depths = transponder_depth+zi;
sound_speed = c0;
sound_speed_gradient = ss_grad;




