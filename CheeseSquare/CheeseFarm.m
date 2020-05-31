classdef CheeseFarm
    % CheeseFarm combine several cheese square file into one object
    properties
        Sources = {}; % all CheeseSquare objects
        DateTime = []; % timestamp of each CheeseSquare object
        Meta = struct(); % aux data
    end
    properties (Dependent)
        Hierarchy; % hierarchy data
    end
    
    properties (Hidden)
        Hierarchy_ = [];
    end
    
    methods
        function farm = CheeseFarm(cheese, varargin)
            % CHEESEFORM the constructor loads the cheese square files into
            % memory
            %   farm = CheeseFarm(cheese, days) loads the experiment in 
            %   'cheese' and all additional 'days' of the experiment
            %   ('cheese' can either be an CheeseSquare object or video
            %   filename)
            
            % load all objects
            farm.Sources = CheeseSquare.LoadGroup(cheese, varargin{:});
            % set timestamp for each object
            for i=1:length(farm.Sources)
                try
                    farm.DateTime(i) = farm.Sources{i}.DateTime;
                catch
                    farm.DateTime(i) = nan;
                end
            end
            % copies meta data from objects
            farm.Meta = Q.cpfield(farm.Sources{1}.Meta, farm.Meta, {'GroupType', 'GroupId', 'CameraId', 'CreationDate', 'Remarks'});
        end
        
        function value = get.Hierarchy(farm)
            % return hierarchies of all objects in the farm
            if isempty(farm.Hierarchy_)
                farm.Hierarchy_.EloRating = Hierarchy.EloRating(farm);
                ce = 0;
                for i=1:length(farm.Sources)
                    farm.Hierarchy_.ChaseEscape(:, :, i) = farm.Sources{i}.Hierarchy.AggressiveChase.ChaseEscape;
                    ce = ce + farm.Sources{i}.Hierarchy.AggressiveChase.ChaseEscape;
                    farm.Hierarchy.CumulDavidScore(:, i) = Hierarchy.DavidScore(ce);
                    farm.Hierarchy.DailyDavidScore(:, i) = Hierarchy.DavidScore(farm.Sources{1}.Hierarchy.AggressiveChase.ChaseEscape);
                end
            end
            value = farm.Hierarchy_;
        end
        
        function s = Sheltered(farm, dati)
            % is sheltered 
            s = farm.GetTrackingData(dati, 'sheltered');
        end

        function h = Hidden(farm, dati)
            % is hidden
            h = farm.GetTrackingData(dati, 'hidden');
        end
        
        function [x, y] = Position(farm, dati)
            % Position(dati) get coordinates of all mice in time
            % 'dati' of type DateTime
            x = farm.GetTrackingData(dati, 'x');
            y = farm.GetTrackingData(dati, 'y');
        end

        function z = Zone(farm, dati)
            % Zone(dati) get zones of all mice in time
            % 'dati' of type DateTime
            z = farm.GetTrackingData(dati, 'zones');
        end
        
        function Show(farm, varargin)
            % UNDER CONSTRUCTION
            opt = Q.defaultargs(false, varargin ...
                , 'DateTime', [] ...
                );
            
            [~, src, ~] = farm.GetTrackingData(opt.DateTime, 'x');
            farm.Sources{src}.Show(varargin{:});
        end
        
        function [value, src, idx] = GetTrackingData(farm, dati, field)
            % [value, src, idx] = GetTrackingData(farm, dati, field) aux
            % function to get tracking data at specific 'dati' date (of
            % object DateTime). The returned 'value' is according to the
            % requested 'field'. The variables 'src' and 'idx' denote the
            % index of the CheeseSquare file and the corresponding frame.
            d = bsxfun(@minus, dati(:), repmat(farm.DateTime, length(dati), 1)); 
            d(d < 0) = inf;
            src = Q.argmin(d, 2);
            idx = round(DateTime.ToFrames(dati(:)' - farm.DateTime(src), farm.Sources{1}.Video.FrameRate) + 1)'; %#ok<PROP>
            %
            usrc = unique(src);
            if length(usrc) > 1
                valid = false(size(src));
                value = zeros(length(src), length(farm.Sources{usrc(1)}.Tracking.(field)(:, idx(1))));
                for u = usrc(:)'
                    map = src == u;
                    valid(map   ) = farm.Sources{u}.Tracking.valid  (   idx(map));
                    value(map, :) = farm.Sources{u}.Tracking.(field)(:, idx(map))';
                end
            else
                valid = farm.Sources{usrc}.Tracking.valid(idx)';
                value = double(farm.Sources{usrc}.Tracking.(field)(:, idx)');
            end
            value(~valid, :) = nan;
        end
    end
end
