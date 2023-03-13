function [grndFunc,grnd_pts] = robustgfit(grnd_pts,interpolatemethod)
%robustgfit Robust curve fitting using RANSAC
%
%   robustgfit(grnd_pts,interpolatemethod) fits a curve to terrain points,
%   grnd_pts, to return a  curve fit object/function,grndFunc, and terrain 
%   points,grnd_pts, used in the interpolation.
%     interpolatemethod - method used for curve fitting. Takes
%     'smoothingspline','cubicspline' or 'linearinterp'
%
% Lonesome Malambo 08/8/2021, Texas A&M Univeristy

% ransac params
sampleSize = round(0.1*size(grnd_pts,1)); % number of points to sample per trial
maxDistance = 6.25; % max allowable distance for inliers

% fit  model evaluation function
fitLineFcn = @(gpts) fit(gpts(:,1),gpts(:,2),interpolatemethod,'Normalize','on'); 
evalLineFcn = ...   % distance evaluation function
  @(model, gpts) sum((gpts(:, 2) - feval(model, gpts(:,1))).^2,2);

% fit model using ransac
[~, inlierIdx] = ransac(grnd_pts,fitLineFcn,evalLineFcn, ...
  sampleSize,maxDistance,'MaxSamplingAttempts',150);

grnd_pts = [grnd_pts(1,:);grnd_pts(inlierIdx,:);grnd_pts(end,:)];% pad with original begin/end points

% refit using only inliers
grndFunc = fit(grnd_pts(:,1),grnd_pts(:,2),interpolatemethod,'Normalize','on'); % fit ground spline

end

