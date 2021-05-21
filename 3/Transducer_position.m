function [lat_transducer, lon_transducer, elev_transducer] = Transducer_position(lat_GPS,lon_GPS, elev_GPS,t)
    %calculation of the offset
    % translate GPS coords to UTM coords
    %compute dE and dN
   
    t_obs = [];
    if (nargin < 4)
            t_obs = [];
    else
        t_obs = t;
    end
    
    % the WaveGlider was rigged with the GPS antenna installed slightly to
    % the right of the center line, 
    
    
    
    [east_GPS,north_GPS]=d2u(lon_GPS,lat_GPS);
    nobs = length(east_GPS);
    iobs = 1:nobs;
    iwant = 0.5:nobs+.5;
   
    
   east_GPS2 = interp1(iobs,east_GPS,iwant,'linear','extrap');
   north_GPS2 = interp1(iobs,north_GPS,iwant,'linear','extrap');
   east_GPS2 = interp1(iobs,east_GPS,iwant,'pchip','extrap');
   north_GPS2 = interp1(iobs,north_GPS,iwant,'pchip','extrap');
   %up_GPS2 = interp1(iobs,elev_GPS,iwant,'cubic','extrap');
    
    
    dE = north_GPS2(2:end) - north_GPS2(1:end-1);
    dN = east_GPS2(2:end) - east_GPS2(1:end-1);
  
    % angle elevation
    % iuse_azimuth = find(~isnan(dN(1:end) + dE(1:end))); %in case there are some NaN values
    azimuth_WG = atan2(dN(1:end),dE(1:end));%*360/(2*pi) to have azymuth in degrees

    % rotate offset
    north_offset_rotate = cos(azimuth_WG(1:end))*0.485000 - sin(azimuth_WG(1:end))*0.167000;
    east_offset_rotate = sin(azimuth_WG(1:end))*0.485000 + cos(azimuth_WG(1:end))*0.167000;

    %transducer_positions
    east_transducer_offset = east_GPS + east_offset_rotate';
    north_transducer_offset = north_GPS + north_offset_rotate';
    % vertical offsets:
    % 1) from bottom of transducer head to waveglider float rail = 30 cm
    % 2) waveglider float rail to base of lift-ring = -2.5cm
    % 3) base of lift-ring to top of gps antenna spacer = 20.5 cm
    % 4) top of gps antenna spacer to base of gps antenna (height of delrin
    % standoff) = 15.24 cm
    vertical_offset = 30.0 -2.5 + 20.5 + 15.24;
    elev_transducer = elev_GPS - vertical_offset/100;   % cm -> m
    [lat_transducer,lon_transducer] = u2d(east_transducer_offset,north_transducer_offset,55);