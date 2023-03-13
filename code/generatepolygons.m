function poly = generatepolygons(sgi,seg_begin_coords, seg_end_coords, attr_table, segmentwidth)
%generatepolygons Generate shapefile object with specified attributes
%   poly = generatepolygons(sgi,seg_begin_coords, seg_end_coords, ...
%                                               attr_table, segmentWidth)
%   generates a shapefile polygon object,poly, based on:
%        - sgi, segment begin ids - serves as shape unique ids
%        - seg_begin_coords, segment begin lat/lon coordinates
%        - seg_end_coords, segment end lat/lon coordinates
%        - attr_table, shapefile attribute data
%        - segmentwidth, output segment width
%
% Lonesome Malambo 08/8/2020, Texas A&M Univeristy

% insert geometries
poly = repmat(struct('Geometry',[]),1,length(seg_begin_coords(:,1)));
[poly(1:length(seg_begin_coords(:,1))).Geometry] = deal('Polygon');

% remove duplicate uids
attrList = attr_table.Properties.VariableNames;
if isprop(attr_table.Properties.CustomProperties,'type')
    % write segment id as attribute
    flt = not(ismember(attrList,'segment_id'));
    ttyp  = 1;
else
    % write segment id as attribute
    flt = not(ismember(attrList,'segment_id_beg'));
    ttyp  = 0;
end
attrList = attrList(flt);
% fx = ismember(attr_table.beamlevel,"Stong");
% attr_table.beamlevel = fx;
attr_table = table2array(attr_table(:,flt));

for j = 1:(length(seg_begin_coords(:,1)))
    % make seg corner coords
    coordbegin = seg_begin_coords(j,:); % begin coords
    coordend = seg_end_coords(j,:); % end coords
    [seglat,seglon] = makesegmentcoords(coordbegin,coordend,segmentwidth);

    % add coordinates to polygon instance
    poly(j).Lat = seglat;
    poly(j).Lon = seglon;

    % write canopy and ground attributes where necessary
    if ttyp == 1
        % write segment id as attribute
        poly(j).seg_id = double(sgi(j));
    else
        % write segment id as attribute
        poly(j).segment_id_beg = double(sgi(j));
    end

    if ~isempty(attr_table)
        for m = 1:length(attrList)
            cda = attr_table(j,m);
            poly(j).(attrList{m}) = double(cda);
        end
    end

end

end

function [seglat,seglon] = makesegmentcoords(coordbegin,coordend,segwith)
E = wgs84Ellipsoid(); % ref eliipsoid
arclen = km2deg(segwith*0.5/1000);% half seg width

% get azimuth
ltbeg = coordbegin(2);
lnbeg = coordbegin(1);
ltend = coordend(2);
lnend = coordend(1);
az = azimuth(ltbeg,lnbeg, ltend,lnend,E);

% get four corners
az1 = az + 90; % 90 deg deflection
[lat1,lon1] = reckon(ltbeg,lnbeg,arclen,az1);
az1 = az - 90;
[lat2,lon2] = reckon(ltbeg,lnbeg,arclen,az1);
az1 = az + 90;
[lat3,lon3] = reckon(ltend,lnend,arclen,az1);
az1 = az - 90;
[lat4,lon4] = reckon(ltend,lnend,arclen,az1);

Latu = [lat1;lat2;lat3;lat4; lat1];
Lonu = [lon1;lon2;lon3;lon4; lon1];

% ensure no cris-cross polygons are created
K = convhull(Lonu,Latu);

Latu = Latu(K);
Lonu = Lonu(K);
seglat = [Latu; NaN];% complete individual polygon loop
seglon = [Lonu; NaN]; % complete individual polygon loop
end

