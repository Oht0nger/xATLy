
classdef gtrack
    %ATL03/ATL08 Ground Track Class
    %   Defines ground track object from an ATL03/ATL08 granule

    properties
        gtrackid % gtrack id
        atl03 Atl03 % ATL03 object
        atl08 Atl08 % ATL08 object
        beamlevel
    end

    methods(Access = protected, Hidden)
        % get methods for beamlevel and a
        function beamlevel = getbeamlevel(obj)
            % Determine groud track beam strength (strong/weak)

            if ~isempty(obj.atl08)
                % try first to get from atl03
                tf = ismember(obj.atl03.groundtracks,obj.gtrackid);% check if exists
                bl = obj.atl03.beamsequence{tf};

                if strcmp(bl,'Undefined')
                    tf = ismember(obj.atl08.groundtracks,obj.gtrackid);% check in atl08
                    b2 = obj.atl08.beamsequence{tf};
                    if strcmp(b2,'Undefined')
                        beamlevel = b1;
                        return
                    end
                    beamlevel = b2;
                else
                    beamlevel = bl;
                end
            else
                tf = ismember(obj.atl03.groundtracks,obj.gtrackid);% check if exists
                beamlevel = obj.atl03.beamsequence{tf};
            end
        end

        function sig_ph = readsigph(obj,varargin)
            %read_sigph read ATL08 photon-level data
            %
            % readsigph(obj) or obj.readsigph() reads signal photon level
            % data for a ground track, obj,from a linked ATL08 granule and
            % returns a data table,sig_ph.sig_ph is mainly auxiliary for
            % linking ATL08 photon classification values to raw ATL03
            % photons
            %
            % Lonesome Malambo 03/25/2022, Texas A&M Univeristy

            % ensure ATL08 is set
            if isempty(obj.atl08)
                error('ATL08 object not set!')
            end

            % include additional attributes
            numvarargs = length(varargin);
            % set defaults for optional inputs
            optargs = {[]};
            optargs(1:numvarargs) = varargin;
            % Place optional args in memorable variable names
            inc_these_seg_attrib = optargs{:};

            tf = ismember(obj.atl08.groundtracks,obj.gtrackid);% check if exists
            if sum(tf)
                bm = ['/' obj.gtrackid];
                % read ATL08 photon-level data datasets
                segid8 = h5read(obj.atl08.filepath,[bm,'/signal_photons/ph_segment_id']); %segment ID
                pcidx8 = h5read(obj.atl08.filepath,[bm,'/signal_photons/classed_pc_indx']); %class_pc_index
                pcid = h5read(obj.atl08.filepath,[bm,'/signal_photons/classed_pc_flag']); % classed_pc_flag
                data = [segid8 pcidx8 pcid];

                varNames = {'segment_id', 'class_pc_index', 'classed_pc_flag'};
                if ~isempty(inc_these_seg_attrib)
                    attrtype = 'land_seg';
                    if ~iscellstr(inc_these_seg_attrib) %#ok<*ISCLSTR>
                        error('input must be a string cell array')
                    end
                    attdata = obj.getattrdata(attrtype, inc_these_seg_attrib);
                    attlist = attdata.Properties.VariableNames;
                    attdata = table2array(attdata);
                    % get ATL08 land-segment data
                    segment_id_beg = h5read(obj.atl08.filepath,[bm,'/land_segments/segment_id_beg']); %
                    segment_id_end = h5read(obj.atl08.filepath,[bm,'/land_segments/segment_id_end']); %

                    % join land_seg data to individual signal photons
                    land_segs = zeros(length(segid8),size(attdata,2));
                    for i = 1:length(segment_id_beg)
                        tf = and(segid8 >= segment_id_beg(i), segid8 <= segment_id_end(i));
                        land_segs(tf,:) = repmat(attdata(i,:), sum(tf),1);
                    end
                    varNames = horzcat(varNames,attlist);
                    data = horzcat(data,land_segs);
                end

                % convert to table
                data = array2table(data,'VariableNames',varNames);
                data.Properties.Description = 'ATL08 Photon-level Data';

                % add metadata
                data = addprop(data,{'type','producttype','gtrackid','groundtracks','beamLevel'}, ...
                    {'table','table','table','table','table'});
                %bmLevel = obj.atl08.beamsequence{tf};
                data.Properties.CustomProperties.type = 'signal_ph';
                data.Properties.CustomProperties.producttype = obj.atl08.producttype;
                data.Properties.CustomProperties.gtrackid = obj.gtrackid;
                data.Properties.CustomProperties.groundtracks = obj.atl08.groundtracks;
                data.Properties.CustomProperties.beamLevel = obj.beamlevel;

                sig_ph = data;
            else
                error('Specified track not found in granule!')
            end
        end

        function attr_table = getattrdata(obj, attrtype, reqattr)
            %getattrdata Read ground track attribute data in the
            % land_segments, canopy or terrain h5 data groups
            %
            % attr_table = obj.getattrdata(attrtype, reqattr) extracts
            % track attribute data for ground track object, obj. attrtype
            % is the target datagroup ('land_segments'|'canopy'|'terrain')
            % and reqattr is the cell array with list of requested
            % attributes.
            % Example 1: attr_tbl = gt.getattrdata('canopy', {'h_canopy','n_ca_photons'})
            % to extract canopy height (98th percentile height) and number
            % of canopy points in each segment for a gt ground track
            % object.
            %
            % Example 2: attr_tbl = gt.getattrdata('terrain', {'h_te_mean'})
            % to extract mean terrain elevation from terrain datagroup.
            %
            % Lonesome Malambo 08/8/2021, Texas A&M Univeristy

            % ensure ATL08 is set
            if isempty(obj.atl08)
                warning('ATL08 object not set!')
                return
            end

            bm = ['/' obj.gtrackid];
            reqattr = unique(reqattr,'stable');
            j = 0;
            switch attrtype
                case 'land_seg'
                    for i = 1:length(reqattr)
                        try
                            atdat = double(h5read(obj.atl08.filepath,[bm,'/land_segments/',reqattr{i}])); %
                            atdat(atdat>3.0e+38) = NaN;
                            if ~isvector(atdat)
                                atdat = atdat';
                            end
                        catch ME
                            warning(ME.message)
                            fprintf('***** %s not be included in output ******\n',reqattr{i})
                            continue
                        end
                        cl = size(atdat,2);
                        pcFields = makeFields(reqattr{i},cl);
                        attdata(:,j+1:j+cl) = atdat; %append
                        attlist(j+1:j+cl) = pcFields;
                        j = j + cl;
                    end
                case 'canopy'
                    for i = 1:length(reqattr)
                        try
                            atdat = double(h5read(obj.atl08.filepath,[bm,'/land_segments/canopy/',reqattr{i}])); %
                            atdat(atdat>3.0e+38) = NaN;
                            if ~isvector(atdat)
                                atdat = atdat';
                            end
                        catch ME
                            warning(ME.message)
                            fprintf('***** %s not be included in output ******\n',reqattr{i})
                            continue
                        end
                        cl = size(atdat,2);
                        pcFields = makeFields(reqattr{i},cl);
                        attdata(:,j+1:j+cl) = atdat; %append
                        attlist(j+1:j+cl) = pcFields;
                        j = j + cl;
                    end
                case 'terrain'
                    for i = 1:length(reqattr)
                        try
                            atdat = double(h5read(obj.atl08.filepath,[bm,'/land_segments/terrain/',reqattr{i}])); %
                            atdat(atdat>3.0e+38) = NaN;
                            if ~isvector(atdat)
                                atdat = atdat';
                            end
                        catch ME
                            warning(ME.message)
                            fprintf('***** %s not be included in output ******\n',reqattr{i})
                            continue
                        end
                        cl = size(atdat,2);
                        pcFields = makeFields(reqattr{i},cl);
                        attdata(:,j+1:j+cl) = atdat; %append
                        attlist(j+1:j+cl) = pcFields;
                        j = j + cl;
                    end
            end
            sgi = double(h5read(obj.atl08.filepath,[bm,...
                '/land_segments/','segment_id_beg']));
            attdata = [sgi attdata];
            attlist = [{'segment_id_beg'} attlist];
            attr_table = array2table(attdata,"VariableNames",attlist);

            function pcFields = makeFields(fld,ncol)
                if ncol == 1
                    pcFields{1} = fld;
                else
                    formatSpec = strcat(fld, "_%d");
                    sq = 1:ncol;
                    s = compose(formatSpec,sq);
                    pcFields = cell(1,length(s));
                    for k = 1:length(s)
                        pcFields{k} = char(s(k));
                    end
                end
            end

        end


        function [segment_id_beg,seg_begin_coords, seg_end_coords] = seggeostruct(obj, clip_aoi)
            %seggeostruct Extract begin/end latlon coordinates for ATL08 segments
            %
            %   [segment_id_beg,seg_begin_coords, seg_end_coords] = seggeostruct(atl03path,
            %           atl08path, gtrack, clip_aoi) extracts segmend ids, begin
            %   (seg_begin_coords) and end (seg_end_coords) latlon coordinates for a
            %   given ground track (gtrackid)based on specified ATL03 (defined by
            %   atl03path filepath) and ATL08 (atl08path) granules. Retrieved coordinates
            %   may be sub-setted to an area of interest by setting clip_aoi
            %   (filepath to shp/kml aoi file, [] for no clipping).
            %
            % Lonesome Malambo 08/8/2020, Texas A&M Univeristy

            % Read atl08 data
            bm = ['/' obj.gtrackid];
            atl03path = obj.atl03.filepath;
            atl08path = obj.atl08.filepath;
            lat8 = double(h5read(atl08path,[bm,'/land_segments/latitude'])); % latitude
            lon8 = double(h5read(atl08path,[bm,'/land_segments/longitude'])); % longitude
            segment_id_beg = h5read(atl08path,[bm,'/land_segments/segment_id_beg']);        % segment id
            dtmbeg = double(h5read(atl08path,[bm,'/land_segments/delta_time_beg']));% segment delta time begin
            dtmend = double(h5read(atl08path,[bm,'/land_segments/delta_time_end']));% segment delta time end

            % read auxiliary atl03
            lat3 = h5read(atl03path,[bm, '/heights/lat_ph']);
            lon3 = h5read(atl03path, [bm, '/heights/lon_ph']);
            dt3 = h5read(atl03path,[bm, '/heights/delta_time']);

            % get nearest segment begin and end lat/lon
            mdl = KDTreeSearcher(dt3); % nn searcher fitting
            idx = knnsearch(mdl, dtmbeg);
            latbeg = lat3(idx);
            lonbeg = lon3(idx);

            idx2 = knnsearch(mdl, dtmend);
            latend = lat3(idx2);
            lonend = lon3(idx2);

            % remove coicident points;
            tt = (latend - latbeg) ~= 0;
            ff = (lonend - lonbeg) ~= 0;
            valid = or(tt,ff);

            segment_id_beg = segment_id_beg(valid);
            lat8 = lat8(valid);
            lon8 = lon8(valid);

            latbeg = latbeg(valid); % begin lat of segment
            lonbeg = lonbeg(valid); % begin lon of segment

            latend = latend(valid);% end lat of segment
            lonend  = lonend(valid);% end lon of segment

            % restrict to clipper boundary
            if ~isempty(clip_aoi)
                if isfile(clip_aoi)
                    % determine format
                    if endsWith(clip_aoi,'.kml')
                        kms = kml2struct(clip_aoi); % read poly struct
                        tf = inpolygon(lon8, lat8, kms.Lon, kms.Lat);
                    elseif endsWith(clip_aoi,'.shp')
                        kms = shaperead(clip_aoi);% read poly
                        tf = inpolygon(lon, lat, kms.X, kms.Y);
                    else
                        error('Unsupported clipper format! Please use shp or kml')
                    end
                else
                    error('Invalid clipper file')
                end
            else
                kms = [];
            end

            if ~isempty(kms)
                % checks if lat/lon contained in clipper bounds
                if max(tf) == 1
                    segment_id_beg = segment_id_beg(tf);
                    latbeg = latbeg(tf);
                    lonbeg = lonbeg(tf);
                    latend = latend(tf);
                    lonend = lonend(tf);

                    lon8 = lon8(tf);
                    lat8 = lat8(tf);
                else
                    disp('No contained/intersecting ATL08 segments');
                end
            end

            % get smooth line through points and project begin and end points
            % onto smooth curve
            f =fit(lon8,lat8,'smoothingspline','Normalize','on');
            llbeg = [lonbeg,latbeg];
            llend = [lonend,latend];
            id = 1:length(lonbeg);

            seg_begin_coords = (cell2mat(arrayfun(@(x) projpt(f,llbeg(x,:)),id, 'UniformOutput',false)))';
            seg_end_coords = (cell2mat(arrayfun(@(x) projpt(f,llend(x,:)),id, 'UniformOutput',false)))';
            function projpt = projpt(fun, lonlat)
                % project point unto fitted line/curve

                lonlim = [lonlat(1) - 0.00045; lonlat(1) + 0.00045];
                latlim = feval(fun, lonlim);

                vector = [lonlim latlim];

                p0 = vector(1,:);
                p1 = vector(2,:);
                a = [-lonlat(1)*(p1(1)-p0(1)) - lonlat(2)*(p1(2)-p0(2)); ...
                    -p0(2)*(p1(1)-p0(1)) + p0(1)*(p1(2)-p0(2))];
                b = [p1(1) - p0(1), p1(2) - p0(2);...
                    p0(2) - p1(2), p1(1) - p0(1)];

                projpt = -(b\a);
            end
        end
    end

    methods
        % gtrack class contructor
        function obj = gtrack(trackid,atl03,varargin)
            numvarargs = length(varargin);
            % set defaults for optional inputs
            optargs = {[]};
            optargs(1:numvarargs) = varargin;
            atl08src = optargs{:};

            if numvarargs == 0
                assert(isa(atl03,'Atl03'),'Expected object of class Atl03')
                tf = ismember(trackid,atl03.groundtracks);% check if exists
                if tf
                    obj.gtrackid = trackid;
                else
                    error('Invalid ground track id! Valid ids : %s',strjoin(atl03.groundtracks,', '))
                end
            elseif numvarargs == 1
                assert(isa(atl03,'Atl03'),'Expected object of class Atl03')
                assert(isa(atl08src,'Atl08'),'Expected object of class Atl08')
                tf = and(ismember(trackid,atl03.groundtracks),...
                    ismember(trackid,atl03.groundtracks));% check if exists in both
                if tf
                    obj.gtrackid = trackid;
                else
                    error('Invalid ground track id! Valid ids : %s',strjoin(atl03.groundtracks,', '))
                end
            end
            obj.atl03 = atl03;
            obj.atl08 = atl08src;
            obj.beamlevel = obj.getbeamlevel();
        end

        % trackid get and set methods
        function gtrackid = get.gtrackid(obj)
            gtrackid = obj.gtrackid;
        end

        function obj = set.gtrackid(obj,id)
            obj.gtrackid = id;
        end
        % atl03 get and set methods
        function atl03 = get.atl03(obj)
            atl03 = obj.atl03;
        end
        function obj = set.atl03(obj,atl03obj)
            obj.atl03 = atl03obj;
        end
        % atl08 get and set methods
        function atl08 = get.atl08(obj)
            atl08 = obj.atl08;
        end
        function obj = set.atl08(obj,atl08obj)
            obj.atl08 = atl08obj;
        end
        % get methods for beamlevel and a
        function beamlevel = get.beamlevel(obj)
            beamlevel = obj.getbeamlevel();
        end

        function writergt(obj, varargin)
            % Write reference ground track as shp or kml
            %
            % writergt(obj) or obj.writergt() writes a ground track line
            % shapefile for the ground track object, obj, saving it to the
            % same folder as the source ATL03 granule. The output shapefile
            % filepath is automatically generated by appending the ground
            % track id to the ATL03 filepath.
            %
            % writergt(obj,outfile) save the ground track shapefile with
            % the specified filepath, outfile, using ATL03 as the source
            % granule
            %
            % writergt(obj,outfile,src) writes the shapefile using specified
            % source ('ATL03' or 'ATL08')
            %
            % writergt(obj,outfile,src, clip_aoi) writes the shapefile
            % clipped to an area of interest defined by clip_aoi, a
            % shapefile or kml file.
            %
            % Lonesome Malambo 03/25/2022, Texas A&M Univeristy

            % include additional attributes
            numvarargs = length(varargin);
            % set defaults for optional inputs
            outf = strrep(obj.atl03.filepath,'.h5',['_' obj.gtrackid '.shp']);
            optargs = {outf,'ATL03',[]};
            optargs(1:numvarargs) = varargin;
            % Place optional args in memorable variable names
            [outfile,src,clip_aoi] = optargs{:};

            if and(strcmp(src,'ATL08'),isempty(obj.atl08))
                src = 'ATL03';
                warning('ATL08 not set, using ATL03 as soruce...')
            end

            try
                switch src
                    case 'ATL08'
                        % read lat/lon values from h5file
                        bm = ['/' obj.gtrackid];
                        lat = double(h5read(obj.atl08.filepath,[bm,'/land_segments/latitude'])); % latitude
                        lon = double(h5read(obj.atl08.filepath,[bm,'/land_segments/longitude'])); % longitude
                    case 'ATL03'
                        % read lat/lon values from h5file
                        bm = ['/' obj.gtrackid];
                        lat = h5read(obj.atl03.filepath,[bm '/geolocation/reference_photon_lat']);
                        lon = h5read(obj.atl03.filepath,[bm '/geolocation/reference_photon_lon']);
                        ix = h5read(obj.atl03.filepath,[bm '/geolocation/reference_photon_index']);
                        bg = h5read(obj.atl03.filepath,[bm '/geolocation/ph_index_beg']);

                        % restrict to segments with data
                        ix = int64(ix) + bg - 1;
                        tf = ix > 0;
                        lat = lat(tf);
                        lon = lon(tf);

                        xx = 1:5:size(lat,1); % sample every 5th point
                        lat = lat(xx);
                        lon = lon(xx);
                    otherwise
                        error('Unsupported source file type')
                end

                % limit to points in polygon, if specified
                if ~isempty(clip_aoi)
                    if isfile(clip_aoi)
                        % determine format
                        if endsWith(clip_aoi,'.kml')
                            kms = kml2struct(clip_aoi); % read poly struct
                            tf = inpolygon(lon, lat, kms.Lon, kms.Lat);
                            lat = lat(tf);
                            lon = lon(tf);
                        elseif endsWith(clip_aoi,'.shp')
                            kms = shaperead(clip_aoi);% read poly
                            tf = inpolygon(lon, lat, kms.X, kms.Y);
                            lat = lat(tf);
                            lon = lon(tf);
                        else
                            error('Unsupported clipper format! Please use shp or kml')
                        end
                    else
                        error('Clipper file does not exist')
                    end
                end

                if and(~isempty(lat),~isempty(lon))
                    % get smooth line through points
                    f =fit(lon,lat,'smoothingspline','Normalize','on');
                    latm = feval(f,lon);

                    if endsWith(outfile,'.shp')
                        outformat = 'shp';
                    elseif endsWith(outfile,'.kml')
                        outformat = 'kml';
                    end

                    switch outformat
                        case 'shp'
                            %write shapefile output
                            waypts = [latm,lon];
                            [lttrk,lntrk] = track('rh',waypts,'degrees');
                            trk.Geometry = deal('Line');
                            trk.Lat = lttrk;
                            trk.Lon = lntrk;
                            trk.Type = obj.gtrackid;
                            shapewrite(trk, outfile);
                            % write projecttion info
                            prj = strrep(outfile, '.shp','.prj');
                            fid = fopen(prj,'w');
                            prjinfo = ['GEOGCS["GCS_WGS_1984",' ...
                                'DATUM["D_WGS_1984",SPHEROID["WGS_1984",' ...
                                '6378137.0,298.257223563]],' ...
                                'PRIMEM["Greenwich",0.0],' ...
                                'UNIT["Degree",0.0174532925199433]]'];
                            fprintf(fid,'%s\n',prjinfo);
                            fclose(fid);
                            %disp("Output successfully written to > " + outfile)
                        case 'kml'
                            %write kml
                            gtrks = {'gt1l','gt1r','gt2l','gt2r','gt3l','gt3r'};
                            i = ismember(gtrks, obj.gtrackid);
                            cols = {'red', 'blue', 'magenta','cyan',...
                                'yellow','green'};% line cols (kml)
                            kmlwriteline(outfile, latm, lon, 'Color',cols{i},'LineWidth', 3);
                            %disp("Output successfully written to > " + outfile)
                        otherwise
                            error("Wrong format specified! Output format should be 'kml' or 'shp' ")
                    end
                else
                    warning('Empty object, not output written!')
                end
            catch ME
                display(ME.message)
            end
        end

        function writesegpol(obj,varargin)
            % Write ATL08 ground track segment polygons shapefile as shp or
            % kml
            %
            % writesegpol(obj) or obj.writesegpol() writes a shapefile for
            % the ground track object, obj, based on linked ATL03 and ATL08
            % granule objects, saving it to the same folder as the source
            % ATL08 granule. The output shapefile filepath is automatically
            % generated by appending the ground track id to the ATL08
            % filepath. By default, the shapefile attribute table contains
            % the following:
            %   - n_seg_ph: total number of points in segment
            %   - n_ca_photons: number of canopy points in segment
            %   - n_toc_photons: number of top of canopy photons
            %   - n_te_photons: number of terrain photons
            %   - h_canopy: 98th percetile relative canopy height
            %   - h_canopy_abs: absolute canopy height
            %   - h_te_mean: mean terrain elevation
            %
            % writesegpol(obj, outfile) saves the shapefile with specified
            % filepath.
            %
            % writesegpol(obj, outfile, land_seg_attr) also includes
            % additional specified attributes from the land_segments group.
            % land_seg_attr is a string cell array,e.g {'night_flag'}, to
            % include number of  segment night flag in the shapefile
            % attribute table.
            %
            % writesegpol(obj, outfile, land_seg_attr,canopy_metrics) also
            % includes additional canopy height metrics from the canopy
            % data group.canopy_metrics a string cell array
            %
            % writesegpol(obj, outfile, land_seg_attr,canopy_metrics,
            % terrain_metrics) also includes additional terrain metrics
            % from the terrain data group.terrain_metrics a string cell array
            %
            % writesegpol(obj, outfile, land_seg_attr,canopy_metrics,
            % terrain_metrics, clip_aoi) also subsets the output file to an
            % area of interest as defined by clip_aoi - a shapefile/kml
            % file.
            %
            % Lonesome Malambo 03/25/2022, Texas A&M Univeristy

            % ensure ATL08 is set
            if isempty(obj.atl08)
                error('ATL08 object not set!')
            end

            % write ATL08 segment polygons as shp
            numvarargs = length(varargin);
            % create temp path
            outf = strrep(obj.atl08.filepath,'.h5',['_' obj.gtrackid '.shp']);
            % set defaults for optional inputs
            def_land_seg = {'n_seg_ph'};
            def_ca_metrics = {'n_ca_photons','n_toc_photons','h_canopy','h_canopy_abs'};
            def_te_metrics = {'n_te_photons','h_te_mean'};
            optargs = {outf, def_land_seg, def_ca_metrics, def_te_metrics, 11,[]};
            optargs(1:numvarargs) = varargin;

            % Place optional args in memorable variable names
            [outfile, land_seg_attr, canopy_metrics, terrain_metrics,...
                segment_width, clip_aoi] = optargs{:};

            % generate segment geo structure
            [segment_id_beg,seg_begin_coords,...
                seg_end_coords] = seggeostruct(obj,clip_aoi);

            % extract attribute data
            attr_lseg = obj.getattrdata('land_seg', land_seg_attr);
            attr_ca = obj.getattrdata('canopy', canopy_metrics);
            attr_ca.segment_id_beg = [];
            attr_te = obj.getattrdata('terrain', terrain_metrics);
            attr_te.segment_id_beg = [];

            attr_table = [attr_lseg attr_ca attr_te];
            itb = ismember(attr_table.segment_id_beg,segment_id_beg);
            attr_table = attr_table(itb,:);

            % generate segment shape structure
            poly = generatepolygons(segment_id_beg,seg_begin_coords, seg_end_coords, ...
                attr_table, segment_width);

            if endsWith(outfile,'.shp')
                outformat = 'shp';
            elseif endsWith(outfile,'.kml')
                outformat = 'kml';
            end

            switch outformat
                case 'shp'
                    % write shapefile projection file
                    shapewrite(poly, outfile);

                    % wgs84 project wkt
                    prjinfo = ['GEOGCS["GCS_WGS_1984",' ...
                        'DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],' ...
                        'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]'];
                    prj = strrep(outfile, '.shp','.prj');
                    fid = fopen(prj,'w');
                    fprintf(fid,'%s\n',prjinfo);
                    fclose(fid);
                case 'kml'
                    kmlwrite(outfile,poly)
                otherwise
                    error('Unsupported output format! Use shp or kml')
            end
        end

        function rawpc = getrawpc(obj)
            % Read ATL03 ground track 3D data
            %
            % getrawpc(obj) or obj.getrawpc() reads 3D photon data
            % for a ground track, obj,from a linked ATL03 granule and
            % returns a data table, gtdata
            %
            % Lonesome Malambo 03/25/2022, Texas A&M Univeristy

            tf = ismember(obj.atl03.groundtracks,obj.gtrackid);% check if exists
            if sum(tf)
                bm = ['/' obj.gtrackid];
                % geolocation attributes and sig. conf
                alt = double(h5read(obj.atl03.filepath,[bm,'/heights/h_ph'])); % reads elevation
                lat = double(h5read(obj.atl03.filepath,[bm,'/heights/lat_ph']));% reads latitude
                lon = double(h5read(obj.atl03.filepath,[bm,'/heights/lon_ph']));% reads longitude
                dt3 = double(h5read(obj.atl03.filepath,[bm,'/heights/delta_time']));% reads time
                sigc = double(h5read(obj.atl03.filepath,[bm,'/heights/signal_conf_ph']));% signal_conf land
                sigc = sigc(1,:)';

                dtff = dt3 - obj.atl03.refdt(1);
                atd = cumsum([0;diff(dtff)].*(7000)); % approx alongtrack distance

                % Form equivalent seg and ph_indx data
                segid3 = h5read(obj.atl03.filepath,[bm,'/geolocation/segment_id']); %segment ID
                phcount = h5read(obj.atl03.filepath,[bm,'/geolocation/segment_ph_cnt']); %photon count
                segrep = repelem(segid3,phcount);
                idxrep = (arrayfun(@num2seq,phcount,'UniformOutput', false));
                idxrep = [idxrep{:}]';

                % gather data
                data = [double(segrep) double(idxrep) atd dt3 lon lat alt sigc];
                varNames = {'segment_id', 'photon_index','alongtrack_distance',...
                    'deltatime','longitude','latitude', 'elevation', 'signal_conf'};
                % convert to table
                data = array2table(data,'VariableNames',varNames);
                data.Properties.Description = 'ATL03 Raw Photon Data';

                % add metadata
                %refdt = getMindt(h5file,tracks);
                data = addprop(data,{'type','producttype','gtrackid',...
                    'groundtracks','beamLevel','refdt'}, ...
                    {'table','table','table','table','table','table'});
                bmLevel = obj.atl03.beamsequence{tf};
                data.Properties.CustomProperties.type = 'heights';
                data.Properties.CustomProperties.producttype = obj.atl03.producttype;
                data.Properties.CustomProperties.gtrackid = obj.gtrackid;
                data.Properties.CustomProperties.groundtracks = obj.atl03.groundtracks;
                data.Properties.CustomProperties.beamLevel = bmLevel;
                data.Properties.CustomProperties.refdt = obj.atl03.refdt;

                rawpc = data;
            else
                error('Specified track not found in granule!')
            end

            function seq = num2seq(num)
                seq = 1:num;
            end
        end

        % link atl08 photon class method
        function classpc = getclassedpc(obj, varargin)
            % getclassedpc Trace ATL08 classified photons back to the
            % source ATL03 raw photons.
            %
            % getclassedpc(obj) or obj.getclassedpc() traces AT08 classified
            % photons back to the source ATL03 file to produce a classified
            % point cloud data table,classpc. pcdata is table containing
            % latitude, longitude, elevation values for each photon
            % together with ATL03 signal conf and photon class ids
            % (0 = Noise,1 = terrain, 2 = canopy, 3 = top of canopy). The
            % ground track object, obj, must have both the atl03 and atl08
            % properties set.
            %
            % getclassedpc(obj, clip_aoi) or obj.getclassedph(clip_aoi)
            % will only return data matching the specified area of
            % interest. clip_aoi is a filepath to a shapefile or kml file.
            %
            % Lonesome Malambo 08/8/2021, Texas A&M Univeristy

            % ensure ATL08 is set
            if isempty(obj.atl08)
                error('ATL08 object not set!')
            end

            % set defaults for optional inputs
            numvarargs = length(varargin);
            optargs = {[],[]};
            optargs(1:numvarargs) = varargin;

            % Place optional args in memorable variable names
            [inc_these_seg_attrib, clip_aoi] = optargs{:};

            % validate atl03 and atl08 source
            v1 = obj.atl03.orbitnum == obj.atl08.orbitnum;
            v2 = obj.atl03.cyclenum == obj.atl08.cyclenum;
            if ~and(v1,v2)
                error('Incompatible ATL03 and ATL08 granules')
            end

            % get gt info
            a3tracks = obj.atl03.groundtracks;
            a8tracks = obj.atl08.groundtracks;
            t1 = ~ismember(obj.gtrackid, a3tracks);
            t2 = ~ismember(obj.gtrackid, a8tracks);
            sm = (t1+t2);
            if sm == 1
                if t1
                    error('Specified track missing in ATL03 dataset')
                else
                    error('Specified track missing in ATL08 dataset')
                end
            elseif sm > 1
                error('Specified track missing in both dataset')
            end

            % read ATL08 photon-level data datasets
            sig_ph = obj.readsigph(inc_these_seg_attrib);
            a8 = table2array(sig_ph(:,1:2));

            % ATL03 datasets
            rawpc = obj.getrawpc();
            a3 = table2array(rawpc(:,1:2));

            if isempty(clip_aoi)
                % get common data
                [C,ia,ib] = intersect(a3,a8,'rows');
                if ~isempty(C)
                    % creates a table compatible with heights
                    rw = size(rawpc,1); % # rows == rows(heights)
                    cl = size(sig_ph(:,3:end),2);% # columns, exclude 1st 2
                    sig_phx = array2table(NaN*zeros(rw,cl), ...
                        'VariableNames',...
                        sig_ph(:,3:end).Properties.VariableNames);

                    % get matching ATL08 class labels
                    sig_ph = sig_ph(ib,3:end);
                    sig_phx(ia,:) = sig_ph;
                    f = isnan(sig_phx.classed_pc_flag);

                    % reassign unmatched to 0 (noise)
                    sig_phx.classed_pc_flag(f) = 0;

                    % combine photon heights and photons classification
                    classpc = [rawpc sig_phx];
                else
                    classpc = [];
                end
            else
                if ~isempty(clip_aoi)
                    if isfile(clip_aoi)
                        % determine format
                        if endsWith(clip_aoi,'.kml')
                            kms = kml2struct(clip_aoi); % read poly struct
                            tf = inpolygon(rawpc.longitude,...
                                rawpc.latitude, kms.Lon, kms.Lat);
                        elseif endsWith(clip_aoi,'.shp')
                            kms = shaperead(clip_aoi);% read poly
                            tf = inpolygon(rawpc.longitude, ...
                                rawpc.latitude, kms.X, kms.Y);
                        else
                            error('Unsupported format! Use shp or kml')
                        end
                    else
                        error('Invalid clipper file')
                    end
                end
                % get common data
                a3 = a3(tf,:);
                rawpc = rawpc(tf,:);
                [C,ia,ib] = intersect(a3,a8,'rows');
                if ~isempty(C)
                    % creates a table compatible with heights
                    rw = size(rawpc,1); % # rows == rows(heights)
                    cl = size(sig_ph(:,3:end),2);% # columns, exclude 1st 2
                    sig_phx = array2table(NaN*zeros(rw,cl), ...
                        'VariableNames',sig_ph(:,3:end).Properties.VariableNames);

                    % get matching ATL08 class labels
                    sig_ph = sig_ph(ib,3:end);
                    sig_phx(ia,:) = sig_ph;
                    f = isnan(sig_phx.classed_pc_flag);

                    % reassign unmatched to 0 (noise)
                    sig_phx.classed_pc_flag(f) = 0;

                    % combine photon heights and photons classification
                    classpc = [rawpc sig_phx];
                else
                    classpc = [];
                end
            end
            if ~isempty(classpc)
                classpc.Properties.CustomProperties.type = 'classed_ph';
            end
        end

        function npcdata = getnormalizedpc(obj,varargin)
            % getnormalizedpc extract classified point cloud normalized to
            % aboveground level. The routine combines steps of 1) matching
            % ATL08 photon classes to raw ATL03 points, 2) estimating the
            % ground surface and 3) normalizing the point cloud
            %
            % getnormalizedpc(obj) or obj.getnormalizedpc() extracts
            % normalized classified point cloud, npcdata as a table, from
            % the ground track object, obj, by using default parameters.
            % Default parameters include:
            %     inc_these_seg_attrib = [] i.e no additional segment
            %     attributs included
            %     clip_aoi = [] , i.e. no clipping to area of interest
            %     grndInterval = 10,  10 m binning distance used in ground
            %     surface fitting
            %     interpolatemethod = 'smoothingspline'
            %     robustfit = false, i.e no robust ground surface fitting
            %     removeptsbelow = true, i.e. removes points
            %
            % getnormalizedpc(obj, inc_these_seg_attrib) or
            % obj.getnormalizedpc(inc_these_seg_attrib) will only return
            % the normalized data including attributes specified in
            % inc_these_seg_attrib e.g.
            % inc_these_seg_attrib = {'night_flag'};, will include the
            % night flag in the output. Optionally,a clip_aoi (kml/shp) can
            % be set to return data clipped to area of interest (aoi).
            %
            % obj.getnormalizedpc(...,grndInterval) enables control of
            % ground interplation process by setting the binning interval,
            % grndInterval
            %
            % obj.getnormalizedpc(...,grndInterval, interpolatemethod)
            % further allows specifying interpolation method. See
            % fitgroundmodel for supported methods
            %
            % obj.getnormalizedpc(..., removeptsbelow) removes any points
            % with elevations below the fitted ground surface from the
            % returned point cloud when removeptsbelow is true, includes
            % them if removeptsbelow is false.
            %
            % Lonesome Malambo 08/8/2022, Texas A&M Univeristy

            % set defaults for optional inputs
            numvarargs = length(varargin);
            optargs = {[],[],10,'smoothingspline',false,true};
            optargs(1:numvarargs) = varargin;

            % Place optional args in memorable variable names
            [inc_these_seg_attrib, clip_aoi, grndInterval,...
                interpolatemethod,robustfit, removeptsbelow] = optargs{:};

            % get classified pc
            pcdata = obj.getclassedpc(inc_these_seg_attrib,clip_aoi);

            % fit ground surface
            [grndFunc, ~] = fitgroundmodel(pcdata,grndInterval,...
                interpolatemethod,robustfit);
            % normalize to aboveground level
            npcdata = normalizepts(pcdata,grndFunc, removeptsbelow);
        end

        function attrdata = getATL08data(obj,varargin)
            % getATL08data Read ATL08 segnemt level attribute data in the
            % land_segments, canopy or terrain h5 data groups. Set the
            % following inputs to read to control what the method returns
            %    - landseg_attr, cell array, to retrieve attributes from
            %    the land_segments group. Use obj.showattr('land_seg') to
            %    list available attributes in this group.
            %    - ca_attr, cell array, to retrieve attributes under the
            %    canopy group. Use obj.showattr('canopy') to list available 
            %    attributes in this group.
            %    - te_attr, cell array, to retrieve attributes under the
            %    terrain group.Use obj.showattr('terrain') to list available 
            %     attributes in this group.
            % 
            %    attrdata = obj.getATL08data(). This a default call and will 
            %    returns the 'longitude','latitude','n_seg_ph'
            %    from the land_segments group,'h_canopy' and 'n_ca_photons' from
            %    the canopy group and 'h_te_best_fit' and 'n_te_photons' from
            %    the terrain group. 
            %    It is the same as obj.getATL08data(landseg_attr, ca_attr, te_attr),
            %    where landseg_attr = {'longitude','latitude','n_seg_ph'},
            %          ca_attr = {'h_canopy','n_ca_photons'} and
            %          te_attr = {'n_te_photons','h_te_best_fit'}
            %
            %    attrdata = obj.getATL08data({'longitude','latitude'})
            %    returns the longitude and latitude values from the land
            %    segments group, setting the ca_attr and te_attr to their
            %    default values
            %    
            %   Note: Each column in attributes that are not vectors e.g. longitude_20m
            %   for the land_segments group or h_canopy_20m from the canopy
            %   group will be represented separately e.g. longitude_20m_1,
            %   longitude_20m_2,... longitude_20_5, for the five columns from the
            %   longitude_20m attribute.
            %
            % Lonesome Malambo 08/8/2021, Texas A&M Univeristy

            % set defaults for optional inputs
            numvarargs = length(varargin);
            optargs = {{'longitude','latitude','n_seg_ph'},{'h_canopy','n_ca_photons'},...
                {'n_te_photons','h_te_best_fit'}};
            optargs(1:numvarargs) = varargin;
            % Place optional args in memorable variable names
            [landseg_attr, ca_attr, te_attr] = optargs{:};
            
            attrdata = [];
            isset = 0;
            if ~isempty(landseg_attr)
                attrdata = obj.getattrdata('land_seg',landseg_attr);
                isset = 1;
            end

            if isset > 0
                if ~isempty(ca_attr)
                    tdt = obj.getattrdata('canopy',ca_attr);
                    attrdata = join(attrdata, tdt);
                    isset = 1;
                end
            else
                if ~isempty(ca_attr)
                    attrdata = obj.getattrdata('canopy',ca_attr);
                    isset = 1;
                end
            end

            if isset > 0
                if ~isempty(te_attr)
                    tdt = obj.getattrdata('terrain',te_attr);
                    attrdata = join(attrdata, tdt);
                end
            else
                if ~isempty(te_attr)
                    attrdata = obj.getattrdata('terrain',te_attr);
                end
            end
        end

        function showattr(obj, attrtype)
            % List attributes under land_segments, canopy and terrain
            % groups
            %
            % showattr(obj, attrtype) lists ATL08 attributes under a
            % specified attribute group. attrype may either take 'land_seg'
            % ,for attributes under the land_segments group, 'canopy', for
            % attributes under the canopy group and 'terrain for attributes
            % under the terrain group.
            %
            % Lonesome Malambo 04/01/2022, Texas A&M Univeristy

            if nargin == 1
                error("Attribute group not specified!\n%s " , ...
                    "Use 'land_seg','canopy' or 'terrain'")
            end

            % ensure ATL08 is set
            if isempty(obj.atl08)
                warning('ATL08 object not set! Existing...')
                return
            end

            h1 = obj.atl08.info;
            bm = ['/' obj.gtrackid];
            gtdx = find(ismember(extractfield(h1.Groups,'Name'),bm));
            switch attrtype
                case 'land_seg'
                    attr = extractfield(h1.Groups(gtdx).Groups(1).Datasets,'Name');
                case 'canopy'
                    attr = extractfield(h1.Groups(gtdx).Groups(1).Groups(1).Datasets,'Name');
                case 'terrain'
                    attr = extractfield(h1.Groups(gtdx).Groups(1).Groups(2).Datasets,'Name');
                otherwise
                    error("Unsupported attribute group - use 'land_seg', 'canopy', 'terrain")
            end
            dispattr(attr, attrtype)

            function dispattr(attr, attrty)
                v = sprintf('------------------------------- %s group attributes -------------------------------',attrty);
                disp(v)

                chunksize = 4;
                ind = 1:chunksize:length(attr);
                if ind(end) < length(attr)
                    ind(end+1) = length(attr);
                end
                k = 2;
                for i = 1:length(attr)
                    if i ==1
                        s = "[" + int2str(ind(k-1)) + "-" + ...
                            int2str(ind(k))+ "]" + " | " + ...
                            convertCharsToStrings(attr{i});
                    else
                        if i < ind(k)
                            s = s + " | " + convertCharsToStrings(attr{i});
                        else
                            s = s + " | " + convertCharsToStrings(attr{i});
                            disp(s)
                            if k < length(ind)
                                k = k+1;
                                s = "[" + int2str(ind(k-1)+1) + ...
                                    "-" + int2str(ind(k))+ "]";
                            end
                        end
                    end
                end
                disp('*** Check the ATL08 ATBD for attribute descriptions ***')
                disp(repelem('-',length(v)))
            end
        end
    end

end

