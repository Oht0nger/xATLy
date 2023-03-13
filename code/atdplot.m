function atdplot(obj, varargin)
%atdplot Plot along-track data
%
% get plot data

% set defaults for optional inputs
numvarargs = length(varargin);
optargs = {'r',[],[],true};
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[plothis, colors, atdlim, excludenoisepts] = optargs{:};

xl = 'Along-track distance (m)';
yl = 'Elevation (m)';
isclassed = false;

switch plothis
    case 'r'
        if isa(obj,'gtrack')
            pc = obj.getrawpc();
        elseif  isa(obj,'table')
            ptypes = struct('heights',{{'r','s'}},...
                'classed_ph',{{'r','s','c'}},...
                'normalized',{{'r','s','c','n'}});
            ttype = obj.Properties.CustomProperties.type;
            if ~ismember(plothis,ptypes.(ttype))
                error('Plot type %s not supported for this input',plothis)
            end
            pc = obj;
        else
            error('Unsupported input object')
        end

        atd = pc.alongtrack_distance;
        alt = pc.elevation;
        grp = ones(length(atd),1);
        if isempty(colors)
            colors = 'b';
        end
        t = 'Raw ATL03 Point Cloud';
        subt = '';
    case 's'
        if isa(obj,'gtrack')
            pc = obj.getrawpc();
        elseif  isa(obj,'table')
            ptypes = struct('heights',{{'r','s'}},...
                'classed_ph',{{'r','s','c'}},...
                'normalized',{{'r','s','c','n'}});
            ttype = obj.Properties.CustomProperties.type;
            if ~ismember(plothis,ptypes.(ttype))
                error('Plot type %s not supported for this input',plothis)
            end
            pc = obj;
        else
            error('Unsupported input object')
        end

        atd = pc.alongtrack_distance;
        alt = pc.elevation;
        grp = pc.signal_conf;
        if isempty(colors)
            colors = 'krmbg';
        end
        t = 'ATL03 Signal Confidence';
        subt = '0 = Noise, 1 = Buffer, 2 = Low, 3 = Medium, 4 = High';
    case 'c'
        if isa(obj,'gtrack')
            pc = obj.getclassedpc();
        elseif  isa(obj,'table')
            ptypes = struct('heights',{{'r','s'}},...
                'classed_ph',{{'r','s','c'}},...
                'normalized',{{'r','s','c','n'}});
            ttype = obj.Properties.CustomProperties.type;
            if ~ismember(plothis,ptypes.(ttype))
                error('Plot type %s not supported for this input',plothis)
            end
            pc = obj;
        else
            error('Unsupported input object')
        end
        atd = pc.alongtrack_distance;
        alt = pc.elevation;
        grp = pc.classed_pc_flag;
        if isempty(colors)
            colors = 'rgbm';
        end
        t = 'Classified ATL03 Point Cloud';
        subt = '0 = Noise, 1 = Terrain, 2 = Canopy, 3 = Top of canopy';
        isclassed = true;
    case 'n'
        if isa(obj,'gtrack')
            pc = obj.getnormalizedpc([],[],10,'smoothingspline',false,false);
        elseif  isa(obj,'table')
            ptypes = struct('heights',{{'r','s'}},...
                'classed_ph',{{'r','s','c'}},...
                'normalized',{{'r','s','c','n'}});
            ttype = obj.Properties.CustomProperties.type;
            if ~ismember(plothis,ptypes.(ttype))
                error('Plot type %s not supported for this input',plothis)
            end
            pc = obj;
        else
            error('Unsupported input object')
        end

        atd = pc.alongtrack_distance;
        alt = pc.nelevation;
        grp = pc.classed_pc_flag;
        yl = 'Aboveground height (m)';
        if isempty(colors)
            colors = 'rgbm';
        end
        t = 'Aboveground Classified ATL03 Point Cloud ';
        subt = '0 = Noise, 1 = Terrain, 2 = Canopy, 3 = Top of canopy';
        isclassed = true;
    otherwise
        error('Unsupported plot type')
end

tf0 = true(size(atd));
if excludenoisepts
    if isclassed
        tf0 = grp > 0;
    end
end

% filter plot data
if ~isempty(atdlim) % filter based on atd limits
    if isnumeric(atdlim)
        if length(atdlim) == 2
            if atdlim(1) >= atdlim(2)
                error('Order Error: adtrange(2) must be > atdrange(1)')
            end
            tf = and(atd >= atdlim(1), atd<= atdlim(2));
            tf = and(tf,tf0);
        else
            error(['adtrange must be a 2-element ' ...
                'vector defining the lower and upper atd limits'])
        end
    else % filter based on aoi
        if isfile(atdlim)
            % determine format
            if endsWith(atdlim,'.kml')
                kms = kml2struct(atdlim); % read poly struct
                rlon = kms.Lon;
                rlat = kms.Lat;
            elseif endsWith(atdlim,'.shp')
                kms = shaperead(atdlim);% read poly
                rlon = kms.X;
                rlat = kms.Y;
            else
                error('Unsupported clipper format! Please use shp or kml')
            end

            lat = obj.latitude;
            lon = obj.longitude;
            tf = inpolygon(lon, lat, rlon, rlat);
            tf = and(tf,tf0);
        else
            error('Unsupported format! Please use a shp or kml file')
        end
    end
else
    tf = tf0;
end
atd = atd(tf);
alt = alt(tf);
grp = grp(tf);

% plot data
mksize = 7;
ax = gscatter(atd,alt,grp,colors,'.',mksize);
if strcmp(plothis,'r')
    legend('off')
else
    legend(ax,'Orientation','horizontal','Location','best');
end
xlabel(xl)
ylabel(yl)
title(t);
subtitle(subt)

end