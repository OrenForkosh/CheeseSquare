classdef Segs
    % Segs working with segments of an array. The individual segments are
    % reffered to as events
    properties (SetAccess = private)
        Events = struct('data', [], 'beg', 0, 'end', 0);
        Length = 0;
    end
    properties (Dependent)
        Map
        Data
    end
    
    properties (GetAccess = private)
        data_ = [];
        map_ = [];
        events_ = [];
        eventsset = false;
    end
    
    methods
        function obj = Segs(varargin)
            % Constructor create a segs object
            %   obj = Segs(map) creates a segs object from the logical map
            %   where adjacent 'true' entries are assigned to the same
            %   segment.
            %   For examples, Segs([1 0 0 1 1 1 0 0 1 1]) creates a 3
            %   segment (or events) object, where the 1st is of length 1, 
            %   the 2nd of length 3 and the 3rd is of two.
            %
            %   obj = Segs(map, data) Also assigs data to the segments.
            %   'data' should be a vector of the same length as map
            %
            %   obj = Segs(segs) initialize the object from another segs
            %   object
            
            if nargin < 2
                if isa(varargin{1}, 'Segs')
                    obj = varargin{1};
                else
                    obj.map_ = varargin{1};
                    obj.eventsset = true;
                    obj.events_ = obj.FindEvents(obj.map_);
                end
            else
                %obj.events_ = Segs.FindEvents(varargin{1}, varargin{2});
                obj.map_ = varargin{1};
                obj.data_ = varargin{2};
                obj.eventsset = true;
                obj.events_ = obj.FindEvents(obj.map_, obj.data_);
            end
            obj.Length = length(obj.map_);
        end
        
        function value = NumberOfEvents(segs)
            % returns the number of events (or segments)
            value = length(segs.Events);
        end
        
        function str = ToString(segs)
            % write segmentation information as a string
            str = '';
            for i=1:length(segs.Events)
                if i > 1
                    str = sprintf('%s, ', str);
                end
                if segs.Events(i).end - segs.Events(i).beg > 0
                    str = sprintf('%s%s-%s', str, num2str(segs.Events(i).beg), num2str(segs.Events(i).end));
                else
                    str = sprintf('%s%s', str, num2str(segs.Events(i).beg));
                end
            end
        end
        
        function Plot(segs)
            % Plot segmentation to a figure
            plot(segs.Map, '-');
            ylim([0, 1.2]);
            xlim([1, segs.Length]);
        end
        
        function obj = DeleteEvents(obj, idx)
            % delete specific event by index
            if islogical(idx)
                idx = find(idx);
            end
            if isempty(idx)
                return;
            end
            for i=idx
                obj.map_(obj.Events(i).beg:obj.Events(i).end) = false;
            end
            obj.eventsset = false;
        end
        
        function obj = Foreach(obj, func)
            % Run a function on each event
            %   obj = obj.Foreach(func) runs the function handle 'func' on
            %   each event. Places the result in the event's data field
            for i=1:length(obj.Events)
                obj.data_(obj.Events(i).beg:obj.Events(i).end) = func(obj.Events(i).beg, obj.Events(i).end, obj.Events(i).data);
            end
            obj.events_ = obj.FindEvents(obj.map_, obj.data_);
        end
        
        function res = Run(obj, func)
            % Run a function on each event 
            %   res = obj.Foreach(func) runs the function handle 'func' on
            %   each event. Places the result in the 'res' variable.
            res = [];
            for i=1:length(obj.Events)
                r = func(obj.Events(i).beg, obj.Events(i).end, obj.Events(i).data);
                if i == 1
                    res = zeros(length(obj.Events), length(r));
                end
                res(i, :) = r;
            end
        end
        
        function value = get.Map(obj)
            % return the map of events (segments)
            if isempty(obj.map_)
                obj.map_ = false(1, obj.Length);
                for i=1:length(obj.Events)
                    e = obj.Events(i);
                    obj.map_(e.beg:e.end) = true;
                end
            end
            value = obj.map_;
        end

        function obj = set.Map(obj, map)
            % set the map of events (segments)            
            obj.map_ = map;
            if ~isempty(obj.data_) && length(obj.data_) ~= length(map)
                error('sizes of ''data'' and ''map'' do not match');
            end
            obj.Length = length(map);
            obj.eventsset = true;
            obj.events_ = obj.FindEvents(obj.map_, obj.data_);
        end
        
        
        function value = get.Data(obj)
            % returns the data vector
            value = obj.data_;
        end
        
        function obj = set.Data(obj, data)
            % set the data vector
            if length(obj.map_) ~= length(data)
                error('sizes of ''data'' and ''map'' do not match');
            end
            obj.data_ = data;
            obj.eventsset = true;
            obj.events_ = obj.FindEvents(obj.map_, obj.data_);
        end
        
        
        function value = get.Events(obj)
            % Seturn a struct of events. Each entry contains:
            %   beg - starting index
            %   end - ending index
            %   data - data associated with event
            %   Length - length of segment
            
            if ~obj.eventsset
                if ~isempty(obj.data_)
                    obj.events_ = obj.FindEvents(obj.map_, obj.data_);
                else
                    obj.events_ = obj.FindEvents(obj.map_);
                end
                obj.Length = length(obj.map_);
                obj.eventsset = true;
            end
            value = obj.events_;
        end
        
        
        function obj = Close(obj, maxgap)
            % Remove small gaps between segments
            %   obj = obj.Close(maxgaps) all gaps which are smaller or
            %   equal maxgaps will be closed (the segments surrounding it
            %   would be connected)
            obj.eventsset = false;
            b = [obj.Events.beg];
            e = [obj.Events.end];
            merge = (b(2:end) - e(1:end-1) - 1) <= maxgap;
            map = obj.map_;
            for i=find(merge)
                map(e(i)+1:b(i+1)-1) = true;
            end
            obj.Map = map;
        end
        
        function segs = Open(segs, minsize)
            % Remove small segments
            %   obj = obj.Open(minsize) all segs which are shorter than
            %   minsize would be removed
            obj.eventsset = false;
            b = [segs.Events.beg];
            e = [segs.Events.end];
            l = e - b + 1;
            b = b(l>minsize);
            e = e(l>minsize);
            segs.Map = Q.segstomap(segs.Length, b, e);
        end
        
        function segs = Expand(segs, minsize)
            % Connect small segments to their closest neighbor
            %   obj = obj.Expand(minsize) all segs which are shorter than
            %   minsize would be connected to their closes neighbor
            obj.eventsset = false;
            b = [segs.Events.beg];
            e = [segs.Events.end];
            l = e - b + 1;
            m = segs.Map;
            for i=find(l<minsize)
                if i>1
                    d1 = b(i) - e(i-1);
                else
                    d1 = inf;
                end
                if i<length(l)
                    d2 = b(i+1) - e(i);
                else
                    d2 = inf;
                end
                if d1 < d2 && isfinite(d1)
                    m(e(i-1):b(i)) = true;
                elseif d2 < d1 && isfinite(d2)
                    m(e(i):b(i+1)) = true;
                end
            end
            segs.Map = m;
        end
        
        function segs = Pad(segs, padsize)
            % Increases the size of all segments
            %   obj = obj.Pad(padsize) adds 'padsize' to the begining and
            %   end of each segment
            padsize = abs(padsize);
            map = segs.Map;
            for i=1:length(segs.Events)
                e = segs.Events(i);
                map(max(e.beg-padsize, 1):e.beg-1) = true;
                map(e.end+1:min(e.end+padsize, segs.Length)) = true;
            end
            segs.Map = map;
        end
        
        function Patch(segs, c, factor, varargin)
            % Plot the segs as a patch
            if nargin < 3
                factor = 1;
            end
            b = [segs.Events.beg];
            e = [segs.Events.end];
            ax = axis;
            y = ax(3:4);
            y(1) = y(1) + factor * (y(2)-y(1)) / 100;
            y(2) = y(2) - factor * (y(2)-y(1)) / 100;
            for i=1:length(e)
                patch([b(i) e(i) e(i) b(i)], [y(1) y(1) y(2) y(2)], 'w', 'FaceColor', c, 'EdgeColor', c, 'FaceAlpha', .2, varargin{:});
            end
        end
        
        function obj = Resize(obj, newlen, method)
            % Resizes the segments vector
            %   obj = obj.Resize(newlen) changed length of all the segments
            %   vector to newlen (scaling all segments and gaps
            %   accordingly)
            obj.eventsset = false;
            if nargin < 3
                method = 'any';
            end
            if ~isempty(obj.data_)
                obj.data_ = interp1(linspace(0, 1, length(obj.data_)), double(obj.data_), linspace(0, 1, newlen), 'linear');
            end
            if newlen > obj.Length
                map = interp1(linspace(0, 1, obj.Length), double(obj.map_), linspace(0, 1, newlen), 'nearest');
                obj.Map = map;
            else
                %                 res = ceil(obj.Length / newlen);
                %                 map = Q.padright(obj.Map, res * newlen - obj.Length);
                %                 map = reshape(map, res, length(map) / res);
                
                switch method
                    case 'any'
                        map = interp1(linspace(0, 1, obj.Length), double(obj.map_), linspace(0, 1, newlen), 'linear') > 0;
                        %map = any(map);
                    otherwise
                        error(['unknown method ''' method '''']);
                end
                obj.Map = map;
            end
        end
    end
    
    methods (Static = true)
        
        function seg = FromBegEnd(len, begf, endf)
            % obj = FromBegEnd(len, begf, endf) Create a segs object of a
            % length 'len' with segments starting at 'begf' and ending at
            % 'endf'. 'begf' and 'endf' can be vectors
            map = false(1, len);
            index = zeros(1, len);
            for i=1:length(begf)
                map(begf(i):endf(i)) = true;
                index(begf(i):endf(i)) = i;
            end
            seg = Segs(map, index);
        end
        
        function events = FindEvents(map, data)
            % events = FindEvents(map, data) creates a list of events from
            % a segment map. Each entry in the list contains:
            %   beg - starting index
            %   end - ending index
            %   data - data associated with event
            c = diff([0 map(:)' 0]);
            begf = find(c > 0);
            endf = find(c < 0) - 1;
            
            events = repmat(struct('data', [], 'beg', 0, 'end', 0), length(begf), 1);
            for i=1:length(begf)
                if nargin >= 2 && ~isempty(data)
                    if isvector(data)
                        events(i).data = data(begf(i):endf(i));
                    else
                        events(i).data = data(begf(i):endf(i), :);
                    end
                else
                    events(i).data = [];
                end
                events(i).beg = begf(i);
                events(i).end = endf(i);
            end
        end
        
        function varargout = Find(map, data)
            % Creates a vector of events. Deprecated!!! use 'FindEvents'
            % instead
            c = diff([0 map(:)' 0]);
            begf = find(c > 0);
            endf = find(c < 0) - 1;
            
            if nargin >= 2
                events = repmat(struct('data', [], 'beg', 0, 'end', 0), length(begf), 1);
                for i=1:length(begf)
                    if isvector(data)
                        events(i).data = data(begf(i):endf(i));
                    else
                        events(i).data = data(begf(i):endf(i), :);
                    end
                    events(i).beg = begf(i);
                    events(i).end = endf(i);
                end
                varargout{1} = events;
            else
                varargout{1} = begf;
                varargout{2} = endf;
            end
        end
        
        
        function map = ToMap(len, b, e)
            % obj = ToMap(len, begf, endf) Create a segs map of a
            % length 'len' with segments starting at 'begf' and ending at
            % 'endf'. 'begf' and 'endf' can be vectors
            map = false(1, len);
            for i=1:length(b)
                map(b(i):e(i)) = true;
            end
        end
        
        function map = CloseMap(map, maxgap)
            % Remove small gaps between segments
            %   obj = CloseMap(map, maxgaps) all gaps in map which are smaller or
            %   equal maxgaps will be closed (the segments surrounding it
            %   would be connected)
            obj = Segs(map);
            b = [obj.Events.beg];
            e = [obj.Events.end];
            merge = (b(2:end) - e(1:end-1) - 1) <= maxgap;
            map = obj.map_;
            for i=find(merge)
                map(e(i)+1:b(i+1)-1) = true;
            end
        end
        
        function map = OpenMap(map, minsize)
            % Remove small segments
            %   obj = OpenMap(map, minsize) all segs in map which are shorter than
            %   minsize would be removed
            obj = Segs(map);
            b = [obj.Events.beg];
            e = [obj.Events.end];
            l = e - b + 1;
            b = b(l>minsize);
            e = e(l>minsize);
            map = Q.segstomap(obj.Length, b, e);
        end
        
    end
end