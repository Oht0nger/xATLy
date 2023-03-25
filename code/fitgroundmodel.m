function [grndfunc, grnd_pts] = fitgroundmodel(pcdata,varargin)
%fitGroundModel Terrain surface interpolation
%   fitgroundmodel interpolates a curve through points classied as terrain
%   in the ATL08 data
%
%   Input
%     pcdata - table containing point segment_id, delta_time, alngtrack
%            distance, lat/lon, elev etc generated from linkATL08
%     grndinterval - atd bin size or spacing used in the interpolation
%     interpolatemethod - method used for curve fitting. Takes
%     'smoothingspline','cubicspline' or 'linearinterp'
%     rubustfit - logical value indicating whether to apply robust
%     selection and fitting using RANSAC.
%
%   Output
%     grndfunc - curve fit object/function
%     grnd_pts - terrain points used in the interpolation
%
% Lonesome Malambo 08/8/2021, Texas A&M Univeristy

% set defaults for optional inputs
numvarargs = length(varargin);
optargs = {10,'smoothingspline',false};
optargs(1:numvarargs) = varargin;
 
% Place optional args in memorable variable names
[grndinterval,interpolatemethod, robustfit] = optargs{:};

% get initial ground estimate
[grnd_atd,grnd_elevation,atdlimits] = extractgroundpts(pcdata);
if or(isempty(grnd_atd),isempty(grnd_elevation))
    error('Insuffcient ground points for interpolation')
end
edg0 = min(grnd_atd):grndinterval:max(grnd_atd);
grnd_pts = zeros(length(edg0)+10,2);
kk = 1;
for b = 2:length(edg0)
    if b == length(edg0)
        fg = and(grnd_atd >= edg0(b-1), grnd_atd <= edg0(b));
    else
        fg = and(grnd_atd >= edg0(b-1), grnd_atd < edg0(b));
    end
    if sum(fg) >= 1
        mat = (edg0(b-1)+edg0(b))/2;
        bindata = grnd_elevation(fg);
        gh =  mean(bindata); %prctile(bindata,50);
        grnd_pts(kk,:) = [mat gh];
        kk = kk + 1;
    else
        continue
    end
end
grnd_pts(kk:end,:) = [];

mn = min(grnd_atd);
mx = max(grnd_atd);

% include endpoints, if necessary to cover all points
if mn > atdlimits(1)
    ga = [atdlimits(1) grnd_pts(1,2)];
    grnd_pts = [ga; grnd_pts];
end

if mx < atdlimits(2)
    ga = [max(atdlimits) grnd_pts(end,2)];
    grnd_pts = [grnd_pts;ga];
end

% fit smooth spline to ground points
if robustfit
    [grndfunc,grnd_pts] = robustgfit(grnd_pts,interpolatemethod);
else
    grndfunc = fit(grnd_pts(:,1),grnd_pts(:,2),interpolatemethod,'Normalize','on'); % fit spline
end

end

