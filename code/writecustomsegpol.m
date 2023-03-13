function writecustomsegpol(hmetrics, outshp, segmentwidth, clip_aoi)
%writecustomsegpol Write custom track segment polygon shapefile
%
%   writecustomsegpol(hmetrics, outshp, segmentwidth, clip_aoi) writes a 
%   segment polygon shapefile associated with hmetrics to a path specified by
%   outshp. Output polygons are generated with a width specified with 
%   segmentwidth. Generated shapefile may be restricted to an area of interest
%   by specifying a clipping file, clip_aoi, in kml or shp format. Output
%   shapefile is referenced to WGS84(geo).
%
% Lonesome Malambo 08/8/2021, Texas A&M Univeristy

% Get bounding box for target area
if ~isempty(clip_aoi)
    if isfile(clip_aoi)
        % determine format
        if endsWith(clip_aoi,'.kml')
            kms = kml2struct(clip_aoi); % read poly struct
        elseif endsWith(clip_aoi,'.shp')
            kms = shaperead(clip_aoi);% read poly 
        else
            error('Unsupported clipper format! Please use shp or kml')
        end
    else
        error(['Clipper file does not exist or is not located ' ...
            'on the specified path'])
    end
else
    kms = [];
end

% get values
lonbeg = hmetrics.seg_lon_beg;% projected begin lon of segment
latbeg = hmetrics.seg_lat_beg;
lonend = hmetrics.seg_lon_end;% projected end lon of segment
latend = hmetrics.seg_lat_end;

% limit to points in polygon, if specified
if ~isempty(kms)
    if endsWith(clip_aoi,'.kml')
        tf = inpolygon(lonend, latend, kms.Lon, kms.Lat);
    else
        tf = inpolygon(lonend, latend, kms.X, kms.Y);
    end
    lonbeg = lonbeg(tf);% projected begin lon of segment
    latbeg = latbeg(tf);
    lonend = lonend(tf);% projected end lon of segment
    latend = latend(tf);
    hmetrics = hmetrics(tf,:);
end

if isempty(hmetrics)
    return
end

% write polygon shapes
seg_begin_coords = [lonbeg latbeg];
seg_end_coords = [lonend latend];
poly = generatepolygons(hmetrics.seg_id,seg_begin_coords, seg_end_coords, ...
                                               hmetrics, segmentwidth);
shapewrite(poly, outshp);

% write shapefile projection file
% wgs84 project wkt
prjinfo = ['GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",' ...
    '6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",' ...
    '0.0174532925199433]]'];
prj = strrep(outshp, '.shp','.prj');
fid = fopen(prj,'w');
fprintf(fid,'%s\n',prjinfo);
fclose(fid);

end

