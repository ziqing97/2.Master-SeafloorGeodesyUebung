	function [x,y,zone,zonelims] = d2u(lon,lat,zone)
%
%	Convert geographic lat-lon degrees to metres
% 
% usage: [x,y,zone,zonelims] = d2u(lon,lat,zone)
%
% d2u converts geographic lon-lat degrees in UTM x-y
% meters. If the zone is passed (string) this zone is 
% used for the conversion otherwise the formal UTM zone is
% calculated based on the median longitude (and latitude)
% are true UTM coordinates are returned.
%

if nargin < 3
    zone=[];
end

if isempty(zone)
    zone = utmzone(nanmedian(lat(:)),nanmedian(lon(:)));
end

zonelims = utmzone(zone);
mstruct = defaultm('utm');
mstruct.zone = [zone];
mstruct.geoid = [6378137.000000 0.0818191910435];   % GRS80 ellipsoid
mstruct = defaultm(utm(mstruct));
[x,y] = mfwdtran(mstruct,lat,lon);
