classdef IC2
    %ICESat-2 Granule Super Class

    properties
        beamsequence % beam strength sequence
        filepath % granule filepath
        producttype % product type ATL03 or ATL08
        orbitnum % orbit_number
        cyclenum % cycle_number
        sc_orient % sc_orient
    end

    properties(Abstract)
        groundtracks % list of available ground tracks
    end

    properties(SetAccess = private,Hidden) %, GetAccess = private
        info % h5info struct
    end

    methods(Abstract)
        gettracklist(obj)
        isintersecting(obj, aoifile)
        showinfo(obj)
    end

    methods
        % constructor
        function obj = IC2(filepath)
            %ATL03 Construct an instance of this class
            if nargin == 1
                if isfile(filepath)
                    obj.filepath = filepath;
                    obj.info = h5info(obj.filepath);
                end
            end
        end

        % get/set methods
        function filepath = get.filepath(obj)
            filepath = obj.filepath;
        end
        function obj = set.filepath(obj,filepath)
            if isfile(filepath)
                obj.filepath = filepath;
                % trick set to update info also, improves performance
                tt = obj;
                tt.info = h5info(obj.filepath);
                obj = tt;
            end
        end

        function orbitnum = get.orbitnum(obj)
            orbitnum = h5read(obj.filepath,'/orbit_info/orbit_number');
        end

        function cyclenum = get.cyclenum(obj)
            cyclenum = h5read(obj.filepath,'/orbit_info/cycle_number');
        end

        function sc_orient = get.sc_orient(obj)
            try
                sc_orient = h5read(obj.filepath, '/orbit_info/sc_orient');
            catch
                sc_orient = 'missing';
            end
        end

        function producttype = get.producttype(obj)
            producttype = obj.getproducttype();
        end

        function beamsequence = get.beamsequence(obj)
            beamsequence = getbeamseq(obj);
        end

        function obj = set.beamsequence(obj, bmseq)
            tt = obj;
            if length(bmseq) ~= length(tt.groundtracks)
                error('Beam sequence must be same length as ground tracks')
            else
                obj.beamsequence = bmseq;
            end
        end
  
        % helper functions
        function prodtype = getproducttype(obj)
            % Determine product type based on short_name id
            %
            prodtype = h5readatt(obj.filepath,'/','short_name');
            prodtype = prodtype(1:5);
        end

        function bmseq = getbeamseq(obj)
            % Get granule beam strength sequence
            gtrks = {'gt1l','gt1r','gt2l','gt2r','gt3l','gt3r'};% full list
            tf = ismember(gtrks,obj.groundtracks); % restrict to available
            sco = obj.sc_orient; %sc_orient
            if sco == 1
                bmseq = {'Weak','Strong', 'Weak','Strong', 'Weak', 'Strong'};
            elseif sco == 0
                bmseq = {'Strong', 'Weak','Strong', 'Weak','Strong', 'Weak'};
            else
                bmseq = {'Undefined','Undefined', 'Undefined','Undefined',...
                    'Undefined', 'Undefined'};
            end
            bmseq = bmseq(tf);
        end

    end
end