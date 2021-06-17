function [t_oneway,r_oneway,x_range,layers] = ray_trace_chadwell(sound_speeds,layer_depths,receiver_depth,geodesic_range,source_depth,theta_check)
%
% [t_oneway,r_oneway,x_range,layers] = ray_trace_chadwell(sound_speeds,layer_depths,receiver_depth,geodesic_range,source_depth,theta_check)
%
% ray_trace_chadwell takes a profile of sound speeds and their depths and
% calculates the one-way travel time and ray-path range from the source
% offset a distance geodesic_range from the receiver.
%
% ("layers" is a structure with the layered model resulting info for use in
% bug-check, plotting, or forming a linear perturbation model)
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
sound_speeds(layer_depths < min(z1) | layer_depths > max(z2)) = [];
layer_depths(layer_depths < min(z1) | layer_depths > max(z2)) = [];

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


ntheta = length(theta_check);

k = repmat(cos(theta_check)./sound_speeds(1),nlayers,1);  % dims: 1 x ntheta. Constant for each sample launch_angle

%     w=NaN.*ones(nlayers+1,ntheta);
%     w(1:nlayers,:) = repmat(sound_speeds(1:nlayers,:)./bvect,1,ntheta);  % only have nlayers of b - need to specially define the final interface value:
%     w(end,:) = w(nlayers,:)+layer_thicknesses(end);
% Bugger - this approach for "w" does not work as "w" is not continuous
% across layer interfaces. So need two "w"s one for the top of each layer,
% the other for the bottom. Generate these based on the mid-point of each
% layer using the gradient and mean sound speed for that layer and then
% use the layers thickness to project to the later interfaces:

w1 = repmat(mean_layer_speeds./bvect - 0.5.*layer_thicknesses,1,ntheta);
w2 = repmat(mean_layer_speeds./bvect + 0.5.*layer_thicknesses,1,ntheta);

b = repmat(bvect,1,ntheta);

% from Chadwell & Sweeney:
A1 = 1./(b .* k);                % dims: nlayers -1  x ntheta
A2 = real(sqrt(1 - b.^2.*w2.^2.*k.^2));    % dims: nlayers -1  x ntheta
A3 = real(sqrt(1 - b.^2.*w1.^2.*k.^2));    % dims: nlayers -1  x ntheta
xdiff = A1.*(-A2 + A3);
tdiff = (atanh(A3)-atanh(A2))./b;   % based on Matlab integral
%    tdiff = log((w2./w1).*(1+A3)./(1+A2))./b;   % based on Chadwell's integral tables (equivalent to the atanh version above
%    tdiff = log((w2./w1).*(1+A2)./(1+A3))./b;   % based on Chadwell's eqn 33: *wrong*

turning_point_check = b.^2.*w2.^2.*k.^2;
iturn = find(turning_point_check>=1);
if ~isempty(iturn)
    % turning_point_check:
%     xdiff(iturn) = 2*A0(iturn).*(A1(iturn);
%     tdiff(iturn) = 2*log((1+A1(iturn))./(k(iturn).*ci(iturn))))./abs(b(iturn));
    % we don't care really - want these paths to be "NaN"?
    xdiff(iturn) = NaN;
    tdiff(iturn) = NaN;
end

ibzero = find(abs(b)==0);        %  want to set this to a threshold equivalent to zero?
if ~isempty(ibzero)
    dz = repmat(layer_thicknesses,1,ntheta);
    c = repmat(mean_layer_speeds,1,ntheta);
    
    tdiff(ibzero) = dz(ibzero)./(c(ibzero).*sqrt(1-k(ibzero).^2.*c(ibzero).^2));
    xdiff(ibzero) = c(ibzero).^2.*k(ibzero).*tdiff(ibzero);
    
end

x_range = sum(xdiff);
t_oneway = sum(tdiff);

rdiff = diag(mean_layer_speeds)*tdiff;
r_oneway = sum(rdiff);

% setup the layers structure:
layers = [];
layers.algorithm = 'Chadwell & Sweeney';
layers.note = 'sea-surface sensor at top, sea-floor sensor at bottom';
layers.depths = layer_depths;
layers.mean_sound_speeds = mean_layer_speeds;
layers.delta_time = tdiff;
layers.effective_raypath_length = rdiff;  
layers.effective_horizontal_distance = xdiff;   

