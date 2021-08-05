function [t_oneway,r_oneway,x_range] = ray_trace_julian(sound_speeds,layer_depths,receiver_depth,geodesic_range,source_depth,theta_check,earth_radius)
%
% [t_oneway,r_oneway,x_range] = ray_trace_julian(sound_speeds,layer_depths,receiver_depth,geodesic_range,source_depth,theta_check,earth_radius)
%
%

%% 

R = earth_radius-layer_depths;
layer_thicknesses = diff(layer_depths);
mean_layer_speeds = sound_speeds(1:end-1)+diff(sound_speeds)./2;


% now need to find the A, B that provide the correct sound_speeds at the
% interfaces:
dR2 = diff(R.^2);
dV = sound_speeds(1:end-1)-sound_speeds(2:end);
B = dV./dR2;
A = sound_speeds(1:end-1)+B.*R(1:end-1).^2;
%A = (sound_speeds(1:end-1)+B.*sound_speeds(1:end-1)

ntheta = size(theta_check,2);
nlayers = size(B,1);

k = repmat(cos(theta_check)./sound_speeds(1),nlayers,1);  % dims: nlayers x ntheta. Constant for each sample launch_angle
P=k;

TT = NaN.*ones(nlayers,ntheta);
DEL = NaN.*ones(nlayers,ntheta);
DDELP = NaN.*ones(nlayers,ntheta);

for itheta = 1:ntheta
    for ilayer = 1:nlayers
        
        [TT(ilayer,itheta),DEL(ilayer,itheta),DDELP(ilayer,itheta)]= julian(A(ilayer),B(ilayer),R(ilayer),P(itheta));
    end
end

t_oneway = sum(TT);
x_range = sum(DEL);
launch_angle =[];
rdiff = diag(mean_layer_speeds)*TT;
r_oneway = sum(rdiff);

%%

return


