classdef Atl08 < IC2 & matlab.mixin.CustomDisplay
    %ICESat-2 ATL08 Granule Class

    properties
        groundtracks
    end

    properties(SetAccess = private, Hidden)
        gtsegs = [];
        metas = struct([]);
    end

    methods(Hidden)
        function value = properties(obj)
            value = builtin("properties", obj);
            ln = length(value);
            ids = [3:ln 1:2];
            value = value(ids); %custom sort -hard coded
            if nargout == 0
                disp(value); %disp( builtin("properties", obj));
            end
        end
    end

    methods(Access = protected)
        function group = getPropertyGroups(obj)
            props = properties(obj);
            group = matlab.mixin.util.PropertyGroup(props);
        end

        function gtdat = gettrackdata(obj)
            % Get per track data/attributes
            %
            gtdat = zeros(length(obj.groundtracks),14);
            for i = 1:length(obj.groundtracks)
                bm = ['/' obj.groundtracks{i}];
                try
                    dt3b = double(h5read(obj.filepath,[bm,'/land_segments/delta_time_beg']));% reads time
                    dt3e = double(h5read(obj.filepath,[bm,'/land_segments/delta_time_end']));% reads time
                    lat = double(h5read(obj.filepath,[bm,'/land_segments/latitude'])); % latitude
                    lon = double(h5read(obj.filepath,[bm,'/land_segments/longitude'])); % longitude

                    nca = sum(double(h5read(obj.filepath,[bm,'/land_segments/canopy/n_ca_photons']))); % n ca
                    ntoc = sum(double(h5read(obj.filepath,[bm,'/land_segments/canopy/n_toc_photons']))); % n toc
                    nte = sum(double(h5read(obj.filepath,[bm,'/land_segments/terrain/n_te_photons']))); % n te

                    hcam = double(h5read(obj.filepath,[bm,'/land_segments/canopy/h_min_canopy'])); % h ca
                    hcam(hcam > 1e38) = [];
                    hcam = min(hcam);
                    hcax = double(h5read(obj.filepath,[bm,'/land_segments/canopy/h_max_canopy'])); % h ca
                    hcax(hcax > 1e38) = [];
                    hcax = max(hcax);

                    htem = double(h5read(obj.filepath,[bm,'/land_segments/terrain/h_te_min'])); % h te
                    htem(htem > 1e38) = [];
                    htem = min(htem);
                    htex = double(h5read(obj.filepath,[bm,'/land_segments/terrain/h_te_max'])); % h ca
                    htex(htex > 1e38) = [];
                    htex = max(htex);

                    gtdat(i,:) = [dt3b(1) dt3e(end) lon(1) lat(1) lon(end) lat(end) length(dt3b) nca ntoc nte hcam hcax htem htex];
                catch ME
                    disp(ME.message)
                    continue
                end
            end

        end
    end

    methods
        function obj = Atl08(fpath)
            %ATL03 Construct an instance of this class
            %   Detailed explanation goes here
            if nargin == 0
                args{1} = '';
            else
                if isfile(fpath)
                    args{1} = fpath;

                else
                    error('Invalid file path!')
                end
            end
            obj = obj@IC2(args{:});
            if ~strcmp(obj.producttype,'ATL08')
                error('Invalid ATL08 file')
            end
        end

        function groundtracks = get.groundtracks(obj)
            groundtracks = obj.gettracklist();
        end

        function gtsegs = get.gtsegs(obj)
            ds = obj.gettrackdata();
            gtsegs = ds(:,7)';
        end
        
        function metas = get.metas(obj)
            % Extracts some metadata attributes to provide input for
            % showinfo()

            if isempty(obj.metas) % get the attr
                obstr = obj.info.Groups(1).Groups(3).Attributes;
                metas(1).id = obstr(ismember(extractfield(obstr,'Name'),'fileName')).Value;% id; acquisition
                metas(1).verid = obstr(ismember(extractfield(obstr,'Name'),'VersionID')).Value;
                metas(1).orbitnum = obj.orbitnum; % orbitnum
                metas(1).cyclenum = obj.cyclenum; % cyclenum
                metas(1).sc_orient = obj.sc_orient;% sc_orient
                metas(1).gtracks = obj.groundtracks; % tracks
                metas(1).beamseq = obj.beamsequence;% beam sequence
                metas(1).gtnsegs = obj.gtsegs;% npts tracks
            end
        end

        function tracks = gettracklist(obj)
            % Extract available tracks from ATL08 granule

            % get file info and extract name field
            h1 = obj.info;
            tracks = repmat({''},1,6);
            k = 1;
            for i = 1:length(h1.Groups)
                if startsWith(h1.Groups(i).Name,'/gt')
                    grp = extractfield(h1.Groups(i).Groups,'Name');
                    ht0 = max(contains(grp,'land_segments'));
                    ht1 = max(contains(grp,'signal_photons'));
                    if (ht0 + ht1) == 2 % enure has land_seg and sigph
                        tracks(k) = strrep(extractfield(h1.Groups(i),'Name'),'/','');
                        k = k+1;
                    end
                end
            end
            tracks(k:end) = [];
        end

        function TF = isintersecting(obj, aoifile)
            % isintersecting Determine if ATL08 granule covers area
            % of interest specified by aoifile
            % TF = obj.isintersecting(aoifile) determines whether
            % ATL08 obj, intersects an area of interest defined by
            % aoifile (a kml/shp format in WGS84 (geo) reference frame).
            % Returns true if intersecting, false otherwise.

            % Lonesome Malambo 08/8/2021, Texas A&M Univeristy

            % limit to points in polygon, if specified
            if ~isempty(aoifile)
                if isfile(aoifile)
                    % determine format
                    if endsWith(aoifile,'.kml')
                        kms = kml2struct(aoifile); % read poly struct
                        x = kms.Lon;
                        y = kms.Lat;
                    elseif endsWith(aoifile,'.shp')
                        kms = shaperead(aoifile);% read poly
                        x = kms.X;
                        y = kms.Y;
                    else
                        error('Unsupported clipper format! Please use shp or kml')
                    end
                else
                    error('File does not exist')
                end
            end

            latpath = '/land_segments/latitude';
            lonpath = '/land_segments/longitude';

            polyg = polyshape(x,y);
            tracklist = obj.gettracklist();
            sel = 0; %zeros(length(tracklist),1);
            for i = 1:length(tracklist)
                bm = ['/' tracklist{i}];
                lat = double(h5read(obj.filepath,[bm,latpath])); % latitude
                lon = double(h5read(obj.filepath,[bm,lonpath])); % longitude

                gtline = [lon(1) lat(1); lon(end) lat(end)];
                [in,~] = intersect(polyg,gtline);
                if ~isempty(in)
                    sel = 1;
                    break
                end
            end

            if sel == 1
                TF = true;
            else
                TF = false;
            end
        end

        function showinfo(obj)
            dtm = gettrackdata(obj);
            disp(' <--------------------------- Dataset Info ------------------------->');
            disp('                                                         ')
            fprintf('   Granule ID: %s\n',obj.metas(1).id);
            fprintf('   Version ID: %s\n',obj.metas(1).verid);
            fprintf('   Orbit/Cycle number: %d/%d\n',obj.metas(1).orbitnum,obj.metas(1).cyclenum);
            if obj.metas(1).sc_orient == 0
                fprintf('   Spacecraft orientation: %d [Backward] \n',obj.metas(1).sc_orient);
            else
                fprintf('   Spacecraft orientation: %d [Forward]\n',obj.metas(1).sc_orient);
            end
            disp('                                                         ')
            fprintf('   Track IDs: %s\n',makestr(obj.metas(1).gtracks));
            fprintf('   Beam sequence: %s\n',makestr(obj.metas(1).beamseq));
            fprintf('   No. segments: %s\n',makestr(obj.metas(1).gtnsegs));
            fprintf('   No. canopy pts: %s\n',makestr(dtm(:,8)));
            fprintf('   No. top.canopy pts: %s\n',makestr(dtm(:,9)));
            fprintf('   No. terrain pts: %s\n',makestr(dtm(:,10)));
            disp('                                                         ')
            fprintf('   Min/max Lat/Lon: [%0.4f, %0.4f] [%0.4f, %0.4f]\n',min(dtm(:,4)),min(dtm(:,3)),max(dtm(:,4)),max(dtm(:,3)));
            fprintf('   Approx. distance: %0.2f km\n',(max(dtm(:,1)) - min(dtm(:,1)))*7);
            fprintf('   Min/max canopy height (m): [%0.1f, %0.1f]\n',min(dtm(:,11)),max(dtm(:,12)));
            fprintf('   Min/max terrain height (m): [%0.1f, %0.1f]\n',min(dtm(:,13)),max(dtm(:,14)));
            disp('                                                         ')
            disp(' <--------------------------- End of Info ------------------------->');

            function st = makestr(cs)
                if isnumeric(cs)
                    cnv = true;
                else
                    cnv = false;
                end
                st = '| ';
                for k = 1:length(cs)
                    if cnv
                        itm = int2str(cs(k));
                    else
                        itm = cs{k};
                    end
                    st = append(st, itm, ' | ');
                end
            end

        end
    end
end