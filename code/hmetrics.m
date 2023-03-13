function hdata = hmetrics(npcdata, seglength, varargin)
% hmetrics Custom Canopy/Terrain height metrics calculation
%
%   hdata = hmetrics(npcdata, seglength) calculates basic canopy
%   height metrics (min, mean, max and std) and terrain metrics (min, mean,
%   max, te_interp (from terrain interpolation using grndFunc), and std)
%   from input table data, npcdata, in segments of length, seglength, along
%   ground track. Canopy height calaculations are limited to segments with
%   9 or more points with max. height >= 2.0 meter. The output, hdata,
%   also contains associated point stats (segment total, canopy and terrain
%   totals) and segment geolocation data,- mid lat/lon, begin and end
%   lat/lon.Table may also include various data flags, if included in prior
%   processes.
%
%   hdata = hmetrics(npcdata, seglength, canopy_percentiles) will
%   also include canopy height percentiles as specified. canopy_percentiles
%   can be a  number or a numeric vector with entries ranging from 0 to 100
%
%   hdata = hmetrics(npcdata, seglength, canopy_percentiles, minpts)
%   limits canopy height calculation to segments with minpts points or more
%
%   hdata = hmetrics(npcdata, seglength, canopy_percentiles,
%   minpts, min_h) limits canopy height calculation to segments with at
%   least minpts points with max canopy height equal or greater than min_h
%
%   hdata = hmetrics(npcdata, seglength, canopy_percentiles,
%   minpts, min_h, includeGround) will include points classified as ground
%   toward the minpts threshold when includeGroud is true (default) or
%   exclude them,includeGround = false
%
% Lonesome Malambo 08/8/2021, Texas A&M Univeristy

numvarargs = length(varargin);

% set defaults for optional inputs
optargs = {[], 11, 2.0, false};
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[canopy_percentiles, minpts, min_h, includegroundpts] = optargs{:};

% get variables
if isempty(npcdata)
    hdata = [];
    return
end

% output variable names
base_vars = {'seg_id', 'seg_atd', 'seg_lon','seg_lat', 'seg_lon_beg',...
    'seg_lat_beg','seg_lon_end', 'seg_lat_end','seg_npts','seg_ca_pts',...
    'seg_te_pts'}; % base variable names
ca_vars = {'h_ca_min','h_ca_mean','h_ca_max','h_ca_std'};%canopy
te_vars = {'h_te_min','h_te_mean','h_te_max','h_te_interp','h_te_std'};%terrain

% determine if flag values are included, make update field list
phx = find(ismember(npcdata.Properties.VariableNames,'classed_pc_flag'));
if (phx + 1) == size(npcdata,2)
    flag_data = [];
    flag_vars = {};
else
    tblcol = size(npcdata,2);
    flag_data = table2array(npcdata(:,(phx+2:tblcol-1)));
    flag_vars = npcdata(:,(phx+2:tblcol-1)).Properties.VariableNames;
end

if isempty(canopy_percentiles)
    pc_flag = false;
    pc_fields = {};
else
    errmgs = ['canopy_percentiles must be a scalar or non-empty vector ' ...
        'of real values > 0 and < 100'];
    if isnumeric(canopy_percentiles)
        if or(min(canopy_percentiles) <= 0, max(canopy_percentiles) >= 100)
            error(errmgs)
        end
    else
        error(errmgs)
    end
    % percentiles var names
    pc_fields = arrayfun(@(x) ['h_ca_p' num2str(x)], canopy_percentiles, 'UniformOutput', false);
    pc_flag = true;
end
varNames = horzcat(base_vars,ca_vars,pc_fields,te_vars,flag_vars);

% filter noise points
phclass = npcdata.classed_pc_flag; % photon class
t = phclass > 0;
npcdata = npcdata(t,:);
if isempty(npcdata)
    hdata = [];
    return
end

elev = npcdata.elevation;
norm_elev = npcdata.nelevation;
atd = npcdata.alongtrack_distance;
phclass = npcdata.classed_pc_flag; % photon class

% get canopy points
if includegroundpts
    tf2 = phclass >= 1; % includes terrain points in canopy height calc
else
    tf2 = phclass > 1;
end
ncpdt = atd(tf2);
ncpalt = norm_elev(tf2);

cpdt = atd(phclass == 1);
cpalt = elev(phclass == 1); % selects only ground points

if ~isempty(flag_data)
    flag_data = flag_data(tf2,:);
end

% estimate approx track direction
lat = npcdata.latitude;
lon = npcdata.longitude;
lat1 = lat(1);
lon1 = lon(1);
lat2 = lat(end);
lon2 = lon(end);
E = wgs84Ellipsoid;
az = azimuth(lat1,lon1,lat2,lon2,E);

% calculate specified height metric in each bin
mxe = ceil(range(ncpdt)/seglength)*seglength + min(ncpdt);
edg = min(ncpdt):seglength:mxe; % bin sequence
hdata = zeros(length(edg)+10,length(varNames));
k = 1;
gfunc = npcdata.Properties.CustomProperties.grndfunc;
for b = 2:length(edg)

    % get indices for segment - normalized
    if b == length(edg)
        nfg = and(ncpdt >= edg(b-1), ncpdt <= edg(b));
    else
        nfg = and(ncpdt >= edg(b-1), ncpdt < edg(b));
    end

    % get indices for segment - unnormalized (terrain)
    if b == length(edg)
        fg = and(cpdt >= edg(b-1), cpdt <= edg(b));
    else
        fg = and(cpdt >= edg(b-1), cpdt < edg(b));
    end

    % compute number of photons (terrain, canopy, top of canopy)
    seg = edg(b) - edg(b-1);
    capts = sum(and(phclass > 1,and(atd >= edg(b-1), atd <= edg(b))));
    if capts < minpts
        [latout,lonout] = reckon(lat1,lon1,seg,az,E);% end lat/lon
        lat1 = latout;
        lon1 = lonout;
        continue
    end
    totpts = sum(and(phclass >= 1,and(atd >= edg(b-1), atd <= edg(b))));
    tepts = sum(and(phclass == 1,and(atd >= edg(b-1), atd <= edg(b))));

    % calculate height for segment
    seg_alt = ncpalt(nfg);
    seg_alt = seg_alt(seg_alt > min_h);
    if ~isempty(seg_alt)
        if pc_flag
            % calculate basic canopy height stats
            basic_ca = [min(seg_alt) mean(seg_alt) max(seg_alt) std(seg_alt)];
            % calculate percentiles
            pct = prctile(seg_alt,canopy_percentiles);% tt.info = h5info(obj.filepath);
            basic_ca_h = horzcat(basic_ca, pct);
        else
            % calculate basic canopy height stats
            basic_ca_h = [min(seg_alt) mean(seg_alt) max(seg_alt) std(seg_alt)];
        end
    else
        % calculate basic canopy height stats
        if pc_flag
            basic_ca_h = NaN.*zeros(1,length(canopy_percentiles)+length(ca_vars));
        else
            basic_ca_h = NaN.*zeros(1,length(ca_vars));
        end
    end

    % calculate terrain basic stats
    satd = edg(b) - 0.5*seg;
    gnpts = sum(fg);
    if gnpts >= 3
        basic_te_h = [min(cpalt(fg)) mean(cpalt(fg)) max(cpalt(fg)) feval(gfunc,satd) std(cpalt(fg))];
    else
        basic_te_h = NaN.*zeros(1,length(te_vars));
        basic_te_h(4) = feval(gfunc,satd);
    end

    % gather all metrics
    hgt = horzcat(basic_ca_h,basic_te_h);
    [latout,lonout] = reckon(lat1,lon1,seg,az,E);% end lat/lon
    mid_atd = mean([edg(b-1),edg(b)]);
    mid_ll = [(lon1+lonout)/2 (lat1+latout)/2];% mid segment lat/lon

    % id, geoloc, flags,pt stats, heights
    if ~isempty(flag_data)
        seg_flags = mode(flag_data(nfg,:),1);
        hdata(k,:) = [k mid_atd mid_ll lon1 lat1 lonout latout totpts capts tepts hgt seg_flags];
    else
        hdata(k,:) = [k mid_atd mid_ll lon1 lat1 lonout latout totpts capts tepts hgt];
    end

    % set begin lat/lon for next segment
    lat1 = latout;
    lon1 = lonout;
    k = k+1;
end
hdata(k:end,:) = [];

% write table object and add custom attributes
hdata = array2table(hdata,'VariableNames',varNames);
hdata.Properties.Description = 'Canopy/Terrain Height Metrics';
hdata = addprop(hdata,{'type'},{'table'});
hdata.Properties.CustomProperties.type = 'hmetrics';

end


