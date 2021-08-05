function [t_oneway,r_oneway,x_range] = ray_trace_hovem(sound_speeds,layer_depths,receiver_depth,geodesic_range,source_depth,theta_check)
%
% [t_oneway,r_oneway,x_range] = ray_trace_hovem(sound_speeds,layer_depths,receiver_depth,geodesic_range,source_depth,theta_check)
%
% ray_trace_hovem takes a profile of sound speeds and their depths and
% calculates the one-way travel time and ray-path range from the source
% offset a distance geodesic_range from the receiver.
%
% N.B. for highest accuracy results the sound-speeds and layer depths
% supplied should have had the (local) flat earth transformation applied.
%

if nargin < 5
    source_depth = 0;   % if not given assume we're pinging from the surface
end

%% start crunching:

% sanity check:
if layer_depths(end)<receiver_depth
    disp(['Problem: sound speed model does not extend below the receiver depth']);
end

% now trim the profile to start at the source depth and stop at the
% receiver depth:
z1 = source_depth; z2 = receiver_depth;

% find the sound speeds at the source and receiver depths:
c1 = interp1(layer_depths,sound_speeds,z1,'linear');
c2 = interp1(layer_depths,sound_speeds,z2,'linear');

% strip off any layers above and below the source and receiver respectively
sound_speeds(layer_depths<z1 | layer_depths > z2) = [];
layer_depths(layer_depths<z1 | layer_depths > z2) = [];

% wrap the profile with the source & receiver depths and interpolated speeds
% (checking for degenerate case where layer boundaries are exactly equal to
% either source of receiver depths)
if layer_depths(1)>z1
    sound_speeds = [c1;sound_speeds];
    layer_depths = [z1;layer_depths];
end
if layer_depths(end)< z2
    sound_speeds = [sound_speeds;c2];
    layer_depths = [layer_depths;z2];
end


% now, how many layers do we have:
nlayers = length(layer_depths)-1;   % n+1 layer interfaces
layer_thicknesses = diff(layer_depths);

% define the sound-speed gradients within each layer:
bvect = diff(sound_speeds)./layer_thicknesses; % dims; nlayers x 1

mean_layer_speeds = sound_speeds(1:end-1)+bvect.*layer_thicknesses/2;

% check: if b==0 then 1./b == NaN! Solution: make b(b==0) very small
% where it is nominally 0:
% izero = find(bvect==0);
% if (izero>0)
%     bfill = 1e-3*randn(length(izero),1);    % minimize any bias by using random perturbation?
%     bfill = 1e-4;    % minimize any bias by using random perturbation
%     bvect(izero)=bfill;
% end

% guesstimate of the launch angle: assume simple trig:
% theta_init = asin(receiver_depth./R);

%
% dtheta = sweep_half_range./nsweeps; % in dwegrees
% % Now define a range of thetas to search over to find the "eigen"ray(s):
% theta_lims = theta_init - pi.*[-sweep_half_range sweep_half_range]./180;
% if (theta_lims(1)>pi/2)
%     theta_lims(1)=pi/2;
% end
% if (theta_lims(2)<=0)
%     theta_lims(2)=0.1;
% end
%
% theta_check = linspace(theta_lims(1),theta_lims(2),nsweeps+1);
% dtheta = abs(diff(theta_lims))./nsweeps;
% %theta_check = theta_init - pi.*[-sweep_half_range:dtheta:sweep_half_range]./180; % work with a +/- 10-deg sweep (starting from high angles!)
% %theta_check(theta_check<=0 | theta_check>=pi/2)=[]; % can't be steeper than 90-deg or shallower than 0!
% disp(['Iteration ',int2str(nloop+1),': Searching between ',num2str(theta_check(1)*180./pi,'%14.10f'),' and ',num2str(theta_check(end)*180./pi,'%14.10f'),' for best launch angle'])

ntheta = length(theta_check);

k = repmat(cos(theta_check)./sound_speeds(1),nlayers,1);  % dims: nlayers x ntheta. Constant for each sample launch_angle

b = repmat(bvect,1,ntheta);

ci = repmat(sound_speeds(1:end-1),1,ntheta); % sound speeds at top of layers
ci_1 = repmat(sound_speeds(2:end),1,ntheta); % sound speeds at bottom of layers:

A0 = 1./(b.*k);
A1 = sqrt(1 - k.^2.*ci.^2);
A2 = sqrt(1 - k.^2.*ci_1.^2);

xdiff = A0.*(A1-A2);

% Hovem write "abs(b)" which for decreasing SS makes tdiff<0 so I added 
% the leading "-". This appears to be wrong - should not be "abs(b)":
%tdiff = -log((ci_1./ci).*(1+A1)./(1+A2))./abs(b);   
tdiff = log((ci_1./ci).*(1+A1)./(1+A2))./b;

turning_point_check = k.^2.*ci_1.^2;
iturn = find(turning_point_check>=1);
if ~isempty(iturn)
    % turning_point_check:
%     xdiff(iturn) = 2*A0(iturn).*(A1(iturn);
%     tdiff(iturn) = 2*log((1+A1(iturn))./(k(iturn).*ci(iturn))))./abs(b(iturn));
    % we don't care really - want these paths to be "NaN"?
    xdiff(iturn) = NaN;
    tdiff(iturn) = NaN;
end

x_range = sum(xdiff);
t_oneway = sum(tdiff);

rdiff = diag(mean_layer_speeds)*tdiff;
r_oneway = sum(rdiff);
