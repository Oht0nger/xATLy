classdef Atl03 < IC2 & matlab.mixin.CustomDisplay
    %ICESat-2 ATL03 Granule Class

    properties
        groundtracks        
    end

    properties(SetAccess = private, Hidden)
        refdt  = []; % min delta time
        gtpoints = [];
        metas = struct([]);
    end

    methods(Hidden)
        % for managing class properties display order
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

        function [refdt,dtm] = getmindt(obj)
            % Determine reference minimum delta time lon/lat values
            %
            % refdt = getMindt(h5file,tracks) returns the minimum delta time value,
            % refdt, from ATL03 granule, h5file, across all available tracks

            dtm = zeros(length(obj.groundtracks),7);
            for i = 1:length(obj.groundtracks)
                bm = ['/' obj.groundtracks{i}];
                try
                    dt3 = double(h5read(obj.filepath,[bm,'/heights/delta_time']));% reads time
                    lat = double(h5read(obj.filepath,[bm,'/heights/lat_ph']));% reads latitude
                    lon = double(h5read(obj.filepath,[bm,'/heights/lon_ph']));% reads longitude
                    dtm(i,:) = [dt3(1) dt3(end) lon(1) lat(1) lon(end) lat(end) length(dt3)];
                catch ME
                    disp(ME.message)
                    continue
                end
            end

            % check if valid
            if min(dtm) == 0
                error('Invalid ATL03 file')
            end
            [~,idx] = min(dtm(:,1),[],1);
            refdt = dtm(idx,:); % minimum delta time value across tracks
        end
    end

    methods
        % constructor
        function obj = Atl03(fpath)
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
            if ~strcmp(obj.producttype,'ATL03')
                error('Invalid ATL03 file')
            end
        end

        function groundtracks = get.groundtracks(obj)
            if isempty(obj.groundtracks)
                groundtracks = obj.gettracklist();
            else
                groundtracks = obj.groundtracks;
            end
        end

        function refdt = get.refdt(obj)
            [rf,~] = obj.getmindt();
            refdt = rf;
        end

        function gtpoints = get.gtpoints(obj)
            [~,ds] = obj.getmindt();
            gtpoints = ds(:,end)';
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
                metas(1).gtnpts = obj.gtpoints;% npts tracks
            end
        end
        function tracks = gettracklist(obj)
            % Extract available tracks from ATL03ngranule

            % get file info and extract name field
            h1 = obj.info;
            tracks = repmat({''},1,6);
            k = 1;
            for i = 1:length(h1.Groups) %length(tracks)
                if startsWith(h1.Groups(i).Name,'/gt')
                    grp = extractfield(h1.Groups(i).Groups,'Name');
                    sm = sum(contains(grp,'heights'));
                    if sm
                        tracks(k) = strrep(extractfield(h1.Groups(i),'Name'),'/','');
                        k = k+1;
                    end
                end
            end
            tracks(k:end) = [];
        end

        function TF = isintersecting(obj, aoifile)
            % isintersecting Determine if ATL03/ATL08 granule covers area
            % of interest specified by aoifile
            % TF = obj.isintersecting(aoifile) determines whether
            % ATL03obj, intersects an area of interest defined by
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

            latpath = '/heights/lat_ph';
            lonpath = '/heights/lon_ph';

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
            [~,dtm] = getmindt(obj);
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
            fprintf('   Number points: %s\n',makestr(obj.metas(1).gtnpts));
            disp('                                                         ')
            fprintf('   Min/Max Lat/Lon: [%0.4f %0.4f] [%0.4f %0.4f]\n',min(dtm(:,4)),min(dtm(:,3)),max(dtm(:,4)),max(dtm(:,3)));
            fprintf('   Approx. Range(km): %0.2f\n',(max(dtm(:,1)) - min(dtm(:,1)))*7);
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