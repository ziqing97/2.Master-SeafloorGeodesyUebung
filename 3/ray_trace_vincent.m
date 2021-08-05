function [t_oneway,r_oneway,launch_angle,x_range] = ray_trace_profile(sound_speeds,layer_depths,receiver_depth,geodesic_range,source_depth)
%
% [t_oneway,r_oneway,launch_angle,x_range] = ray_trace_profile(sound_speeds,layer_depths,receiver_depth,geodesic_range,source_depth)
%
% ray_trace_profile takes a profile of sound speeds and their depths and
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

% geometric range:
R = sqrt(receiver_depth.^2+geodesic_range.^2);

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
izero = find(bvect==0);
if (izero>0)
    bfill = 1e-3*randn(length(izero),1);    % minimize any bias by using random perturbation?
    bfill = 1e-4;    % minimize any bias by using random perturbation
    bvect(izero)=bfill;
end

% guesstimate of the launch angle: assume simple trig:
theta_init = asin(receiver_depth./R);

% define some initial parameters for the interation
sweep_half_range = 30;    % start with a +/- 30-deg sweep (starting from high!)
nsweeps = 100;
convergence_check = 0.001;     % 1 mm (final step should make this even smaller
x_off = 9999;           % initial misfit
x_off_last = x_off;
%% now want to run a loop that iterates in with finer and finer angular
% resolution until we have a ray-path that reaches the bottom of the
% sound-speed stack with a horizontal range within "convergence_check" of
% the provided geodesic distance between source and receiver

nloop = 0;
while (x_off > convergence_check)
    %% and this is the matrix version... maybe?
    
    clear w k theta_check w1 w2
    dtheta = sweep_half_range./nsweeps; % in dwegrees
    % Now define a range of thetas to search over to find the "eigen"ray(s):
    theta_lims = theta_init - pi.*[-sweep_half_range sweep_half_range]./180;
    if (theta_lims(1)>pi/2)
        theta_lims(1)=pi/2;
    end
    if (theta_lims(2)<=0)
        theta_lims(2)=0.1;
    end
    
    theta_check = linspace(theta_lims(1),theta_lims(2),nsweeps+1);
    dtheta = abs(diff(theta_lims))./nsweeps;
    %theta_check = theta_init - pi.*[-sweep_half_range:dtheta:sweep_half_range]./180; % work with a +/- 10-deg sweep (starting from high angles!)
    %theta_check(theta_check<=0 | theta_check>=pi/2)=[]; % can't be steeper than 90-deg or shallower than 0!
    disp(['Iteration ',int2str(nloop+1),': Searching between ',num2str(theta_check(1)*180./pi,'%14.10f'),' and ',num2str(theta_check(end)*180./pi,'%14.10f'),' for best launch angle'])
    
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
    %     w1 = w(1:end-1,:);
    %     w2 = w(2:end,:);
    
    b = repmat(bvect,1,ntheta);
    
    % from Chadwell & Sweeney:
    A1 = 1./(b.*k);                % dims: nlayers -1  x ntheta
    A2 = real(sqrt(1 - b.^2.*w2.^2.*k.^2));    % dims: nlayers -1  x ntheta
    A3 = real(sqrt(1 - b.^2.*w1.^2.*k.^2));    % dims: nlayers -1  x ntheta
    xdiff = A1.*(-A2 + A3);
    tdiff = (atanh(A3)-atanh(A2))./b;   % based on Matlab integral
    %    tdiff = log((w2./w1).*(1+A3)./(1+A2))./b;   % based on Chadwell's integral tables
    %    tdiff = log((w2./w1).*(1+A2)./(1+A3))./b;   % based on Chadwell's eqn 33: *wrong*
    ibzero = find(abs(b)==0);        %  want to set this to a threshold equivalent to zero?
    if ~isempty(ibzero)
        dz = repmat(layer_thicknesses,1,ntheta);
        c = repmat(mean_layer_speeds,1,ntheta);
        
        tdiff(ibzero) = dz(ibzero). / ( c(ibzero). * sqrt( 1 - k(ibzero).^2 * c(ibzero).^2 ) );
        xdiff(ibzero) = c(ibsero.^2).*k(ibzero).*tdiff(ibzero);
        
    end
    
    x_range = sum(xdiff);
    t_oneway = sum(tdiff);
    
    rdiff = diag(mean_layer_speeds)*tdiff;
    r_oneway = sum(rdiff);
    
    % first need to throw out any luanch angles that leads to ray paths never
    % making it to the bottom - if xdiff == 0 then ray never enters the layer,
    % or is perfectly vertical (so check for that case as well);
    check_diffs = xdiff;
    %     check_diffs(xdiff == 0) = NaN; % first cut assume 0 means no ray:
    %     check_diffs(:,cos(theta_check)==0) = 0;    % put zeros back for that column if exactly vertical (not sure this is a practical concern!)
    
    check_ranges = sum(check_diffs); % should now have NaNs for any launch angle that doesn't get to the bottom
    itrapped = find(isnan(check_ranges));
    check_thetas = theta_check;
    check_thetas(itrapped)=[];
    check_ranges(itrapped)=[];
    
    % we've swept through angles from high to low, expect that x_range will initially
    % be shorter than geodesic_range then will step past the "best" angle and
    % become longer. At longer ranges may then shorten again as ray never
    % emerges from a layer...
    ilo = find(check_ranges>geodesic_range,1); % first range that is longer than the desired range
    ihi = ilo-1;   % this should be true! ie previous angle produces a range that is too short...
    
    theta_lo = check_thetas(ilo);
    theta_hi = check_thetas(ihi);
    xdist_lo = check_ranges(ilo);
    xdist_hi = check_ranges(ihi);
    [x_off,imin]=min(abs([xdist_lo xdist_hi] - geodesic_range));
    %x_off = min(abs([xdist_lo xdist_hi] - geodesic_range));
    
    if (x_off >= x_off_last)
        % break
    end
    
    % use newton-raphson method to calculate new best guess for the launch
    % angle:
    theta_init=theta_lo+((geodesic_range)/((xdist_lo-xdist_hi)/(theta_lo-theta_hi)))-((xdist_lo)/((xdist_lo-xdist_hi)/(theta_lo-theta_hi)));
    
    % define tighter limits on the range on launch_angles to search through
    sweep_half_range = 2.5*dtheta;  % in degrees
    
    nloop = nloop+1;
    
    disp(['Iteration ',int2str(nloop),': Found horizontal range ',num2str(x_off,'%8.6f m'),' from target range'])
    disp(['Iteration ',int2str(nloop),': Guessing  ',num2str(theta_init(1)*180./pi,'%14.10f'),' for launch angle'])
    x_off_last = x_off;
    
    
    %% loop interior subsection
    
end

%% final re-calc:
% rerun this final launch angle to get the "best" ray trace:
k = repmat(cos(theta_init)./sound_speeds(1),nlayers,1);  % dims: 1 x ntheta. Constant for each sample launch_angle
% w=NaN.*ones(nlayers+1,1);
% w(1:nlayers,1) = sound_speeds(1:nlayers,:)./bvect;  % only have nlayers of b - need to specially define the final interface value:
% w(end,:) = w(nlayers,1)+layer_thicknesses(end);
%
% recreate b, w etc as vectors

w1 = mean_layer_speeds./bvect - 0.5.*layer_thicknesses; % top of each layer
w2 = mean_layer_speeds./bvect + 0.5.*layer_thicknesses; % bottom of each layer
% w1 = w(1:end-1,1);
% w2 = w(2:end,1);

b = bvect;

% from Chadwell & Sweeney:
A1 = 1./(b.*k);                % dims: nlayers -1  x ntheta
A2 = real(sqrt(1 - b.^2.*w2.^2.*k.^2));    % dims: nlayers -1  x ntheta
A3 = real(sqrt(1 - b.^2.*w1.^2.*k.^2));    % dims: nlayers -1  x ntheta
xdiff = A1.*(-A2 + A3);
tdiff = (atanh(A3)-atanh(A2))./b;
x_range = sum(xdiff);
t_oneway = sum(tdiff);

rdiff = diag(mean_layer_speeds)*tdiff;
r_oneway = sum(rdiff);

launch_angle = theta_init*180./pi;

%%
return

%% OK - loop this to make sure we're setting it all up correctly:
% using the matlab symbolic toolbox to get the integrals:
% Chadwell eqn30 becomes:   -atanh((1 - b^2*k^2*w^2)^(1/2))/b
% and eqn 31 is OK:       -(1 - b^2*k^2*w^2)^(1/2)/(b*k)
%

% scalar version:
% k = cos(theta_check)./sound_speeds(1);  % dims: 1 x ntheta. Constant for each sample launch_angle
% w=NaN.*sound_speeds;
% w(1:nlayers,:) = sound_speeds(1:nlayers,:)./b;  % only have nlayers of b - need to specially define the final interface value:
% w(end,:) = w(nlayers,:)+layer_thicknesses(end);
%
% clear xdiff tdiff rdiff;
% for it = 1:ntheta
%     k = cos(theta_check(it))./sound_speeds(1);  % dims: 1 x ntheta. Constant for each sample launch_angle
%     for il = 1:nlayers
%         X1 = 1./(k.*b(il));
%         X2 = sqrt(1 - k.^2.*b(il).^2.*w(il+1).^2);
%         X3 = sqrt(1 - k.^2.*b(il).^2.*w(il).^2);
%         xdiff(il,it)=X1.*(-X2 + X3);
%         % using the Matlab integral result as it gives the correct numbers:
%         tdiff(il,it) = (atanh(X3)-atanh(X2))./b(il);
%     end
%     rdiff(:,it) = tdiff(:,it).*mean_layer_speeds;
% end

