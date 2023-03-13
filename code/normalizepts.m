function npcdata = normalizepts(pcdata,grndfunc, removeptsbelow)
%normalizepts Normalize points to aboveground level
%
%   normalizepts(pcdata,grndfunc, removeptsbelow) normalizes point
%   data,pcdata, to aboveground level by subtracting the elevation values
%   represented by the ground curve object/fit, grndFunc, create a 
%   normalized point data, npcdata. Points below the surface defined by 
%   grndFunc are removed from npcdata if remotePtsBelow flag is True or 
%   retained if flag is False
%
% Lonesome Malambo 08/8/2021, Texas A&M Univeristy

atd = pcdata.alongtrack_distance;
elev = pcdata.elevation;

ge = feval(grndfunc,atd); % ground estimate
nelev = elev - ge; % heights normalized to ground level
pcdata.nelevation = nelev;

if removeptsbelow
    tc = nelev > 0;
    npcdata = pcdata(tc,:);
else
    npcdata = pcdata;
end
npcdata.Properties.CustomProperties.type = 'normalized';
npcdata = addprop(npcdata,{'grndfunc'}, {'table'});
npcdata.Properties.CustomProperties.grndfunc = grndfunc;

end

