function [t_oneway,r_oneway,launch_angle,x_range,layers] = find_direct_raypath(sound_speeds,layer_depths,receiver_depth,geodesic_range,source_depth,earth_radius,method)
%
% function [t_oneway,r_oneway,launch_angle,x_range] = find_direct_raypath(sound_speeds,layer_depths,receiver_depth,geodesic_range,source_depth,earth_radius,method);
%
% find_direct_raypath takes a sound_speed profile, the depths of the
% receiver and source and optional mean radius (for either spherical
% ray-tracing or earth flattening) and uses the given "method" to find to
% first/direct-path "eigen" ray (to within nominal accuracy) between the
% receiver and source.
%
if nargin < 7
    method = 'chadwell';
end

if nargin < 6
    earth_radius = 6371008; % volumetric mean
end

%% initial parameters for the search
sweep_half_range = 20;
nsweeps = 100;

% direct geometric path length between source and receiver:
R = sqrt((receiver_depth-source_depth).^2+geodesic_range.^2);
% geometric angle between source and receiver:
theta_init = asin(receiver_depth./R);

% Now define a range of thetas to search over to find the "eigen"ray:
theta_lims = theta_init - pi.*[-sweep_half_range sweep_half_range]./180;
if (theta_lims(1)>pi/2)
    theta_lims(1)=pi/2;
end
if (theta_lims(2)<=0)
    theta_lims(2)=0.1;
end

theta_check = linspace(theta_lims(1),theta_lims(2),nsweeps+1);
% define some initial parameters for the interation

convergence_check = 0.001;     % 1 mm (final step should make this even smaller)
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
    
    %     disp(['Iteration ',int2str(nloop+1),': Searching between ',num2str(theta_check(1)*180./pi,'%14.10f'),' and ',num2str(theta_check(end)*180./pi,'%14.10f'),' for best launch angle'])
    % got 4 possible approaches to choose from: hovem and chadwell should
    % be identical and they and vincent require earth flattening. Julian is
    % meant to be a spherical solution
    switch method
        case 'hovem'
            [t_oneway,r_oneway,x_range] = ray_trace_hovem(sound_speeds,layer_depths,receiver_depth,geodesic_range,source_depth,theta_check);
        case 'chadwell'
            [t_oneway,r_oneway,x_range] = ray_trace_chadwell(sound_speeds,layer_depths,receiver_depth,geodesic_range,source_depth,theta_check);
        case 'vincent'
            [t_oneway,r_oneway,x_range] = ray_trace_vincent(sound_speeds,layer_depths,receiver_depth,geodesic_range,source_depth,theta_check);
        case 'julian'
            [t_oneway,r_oneway,x_range] = ray_trace_julian(sound_speeds,layer_depths,receiver_depth,geodesic_range,source_depth,theta_check,earth_radius)
    end
    
    % or is perfectly vertical (so check for that case as well);
    %check_diffs = xdiff;
    %     check_diffs(xdiff == 0) = NaN; % first cut assume 0 means no ray:
    %     check_diffs(:,cos(theta_check)==0) = 0;    % put zeros back for that column if exactly vertical (not sure this is a practical concern!)
    
    check_ranges = x_range; % should now have NaNs for any launch angle that doesn't get to the bottom
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
    theta_init=double(theta_lo+((geodesic_range)/((xdist_lo-xdist_hi)/(theta_lo-theta_hi)))-((xdist_lo)/((xdist_lo-xdist_hi)/(theta_lo-theta_hi))));
    
    % define tighter limits on the range on launch_angles to search through
    
    nloop = nloop+1;
    sweep_half_range = double(2*dtheta);  % in degrees
    
    %     disp(['Iteration ',int2str(nloop),': Found horizontal range ',num2str(x_off,'%8.6f m'),' from target range'])
    %     disp(['Iteration ',int2str(nloop),': Guessing  ',num2str(theta_init(1)*180./pi,'%14.10f'),' for launch angle'])
    x_off_last = x_off;
        
    %% loop interior subsection
    
end

layers = [];
theta_best = theta_init;
% final run with just the best theta:
    switch method
        case 'hovem'
            [t_oneway,r_oneway,x_range] = ray_trace_hovem(sound_speeds,layer_depths,receiver_depth,geodesic_range,source_depth,theta_best);
        case 'chadwell'
            [t_oneway,r_oneway,x_range,layers] = ray_trace_chadwell(sound_speeds,layer_depths,receiver_depth,geodesic_range,source_depth,theta_best);
        case 'vincent'
            [t_oneway,r_oneway,x_range] = ray_trace_vincent(sound_speeds,layer_depths,receiver_depth,geodesic_range,source_depth,theta_best);
        case 'julian'
            [t_oneway,r_oneway,x_range] = ray_trace_julian(sound_speeds,layer_depths,receiver_depth,geodesic_range,source_depth,theta_best,earth_radius)
    end
 launch_angle = theta_best.*180./pi;
% disp(['Found horizontal range ',num2str(x_range,'%10.5f m'),' using  ',num2str(launch_angle,'%14.10f'),' for launch angle'])
