function [grnd_atd,grnd_elev,atd_limits] = extractgroundpts(pcdata)
%extractgroundpts Extract ground points (atd, elev) and data limits
%   [grnd_atd,grnd_elev,atd_limits] = extractgroundpts(pcdata) extracts
%   terrain point values (alongtrack distance, grnd_atd and
%   elevation,grnd_elev) together with atd data limits, atd_limits, to
%   support interpolation of the terrain surface
%
% Lonesome Malambo 08/8/2021, Texas A&M Univeristy

atd = pcdata.alongtrack_distance; % approx alongtrack distance
elev = pcdata.elevation; % elevation
phclass = pcdata.classed_pc_flag; % photon class
tf1 = phclass == 1;
if max(tf1)==1
    grnd_atd = atd(tf1);
    grnd_elev = elev(tf1);
    atd_limits = [min(atd) max(atd)];
else
    grnd_atd = [];
    grnd_elev = [];
    atd_limits = [min(atd) max(atd)];
end

end