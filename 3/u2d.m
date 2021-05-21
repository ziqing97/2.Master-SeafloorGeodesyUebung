function    [lat,lon] = u2d(x,y,izone)
%
% [lat,lon] = u2d(x,y,izone);
%
% u2d converts the utm x-y coordinates into geographic
% lat lon based on the specified UTM zone (either an integer, or
% the string zone as returned by utmzone. Requires the
% mapping toolbox.
%
%

mstruct = defaultm('utm');
if ~isstr(izone)
    mstruct.zone = [int2str(izone)];
else
    mstruct.zone = izone;
end
%mstruct.geoid = almanac('earth','geoid','m','clarke66');
mstruct.geoid = [6378137.000000 0.0818191910435];   % GRS80 ellipsoid
mstruct = defaultm(utm(mstruct));
[lat,lon] = minvtran(mstruct,x,y);
