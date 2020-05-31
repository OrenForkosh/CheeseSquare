classdef CheeseSquare < Serializable
    % CHEESESQUARE Main object type in system. Provides access to tracking
    % data, video file, and preprocessing data.
    properties
        Prefix = [];            % file prefix (such as 'SC.exp0001.day01.cam02')
        Tracking = [];          % Tracking data
        Video = [];             % Video file ('MyVideoReader' object)
        Meta = struct();        % Video meta data
        Hierarchy = struct();   % Domincance hierarchy
        Background = struct();  % Video background data (for background subtraction)
        nSubjects = 0;          % Number of mice
        Colors = struct();      % Mouse colors 
        Analysis = struct();    % Additional analysis (not used anymore)
        ROI = [];               % Regions of interenst 
        UID = char(java.rmi.server.UID().toString); % Unique ID
    end
    
    properties (Dependent)
        DateTime;               % Date of video recording
        Profile;                % Behavioral profile (see 'CheeseNewProfiler')
    end
    
    properties (Hidden)
        % List of class fields that we always copy
        TrackingFields = {'nSubjects', 'Hierarchy', 'ROI', 'Colors', 'Profile'}; 
        Filename = MyFilename('');  % Source file for object
        Profile_ = [];
    end
    
    methods (Static = true)
        function cmap = ZoneColors()
            % ZoneColors Colormap for different ROIs
            z = CheeseSquare.ZoneColormap;
            idx = CheeseSquare.ZoneIdx(z(1:2:end));
            cmap = cat(1, z{idx*2});
        end
        
        function [z, c] = ZoneIdx(name)
            % [index, category] = ZoneIdx(name) Gives the ROI 'name' an index and 
            % category. The categories are: 
            %   1-open, 2-food\water, 3-shelter, 4-high place, 5-block
            if nargin == 0
                z = 13; % max number of zones
                c = 5; % max number of categories
                return;
            end
            if iscell(name)
                z = zeros(1, length(name));
                c = zeros(1, length(name));
                for i=1:length(name)
                    [z(i), c(i)] = CheeseSquare.ZoneIdx(name{i});
                end
                return;
            end
            switch lower(name)
                case lower('Open');         z =  1; c = 1;
                case lower('Feeder1');      z =  2; c = 2;
                case lower('Feeder2');      z =  3; c = 2;
                case lower('Water');        z =  4; c = 2;
                case lower('Water2');       z =  5; c = 2;
                case lower('SmallNest');    z =  6; c = 4;
                case lower('(SmallNest)');  z =  7; c = 3;
                case lower('Labyrinth');    z =  8; c = 4;
                case lower('BigNest');      z =  9; c = 2;
                case lower('(BigNest)');    z = 10; c = 3;
                case lower('Block');        z = 11; c = 5;
                case lower('[Ramp1]');      z = 12; c = 4;
                case lower('[Ramp2]');      z = 13; c = 4;
                otherwise
                    Console.Warning('unrecognized zone: %s', name);
                    z = nan;
                    c = nan;
            end
        end
        
        function z = ZoneColormap()
            % ZoneColormap Colormap for different ROIs
            z = {...
                'Open',         Colors.ParseHex('D94D57'),...
                'Feeder1',      Colors.ParseHex('3EB049'),...
                'Feeder2',      Colors.ParseHex('3A6731'),...
                'Water',        Colors.ParseHex('33A6B3'),...
                'Water2',       Colors.ParseHex('19888D'),...
                'SmallNest',    Colors.ParseHex('CA852A'),...
                '(SmallNest)',  Colors.ParseHex('D6D6D6'),...
                'Labyrinth',    Colors.ParseHex('A5274D'),...
                'BigNest',      Colors.ParseHex('AA6615'),...
                '(BigNest)',    Colors.ParseHex('FFFFFF'),...
                'Block',        Colors.ParseHex('D150B8'),...
                '[Ramp1]',      Colors.ParseHex('9E7EB9'),...
                '[Ramp2]',      Colors.ParseHex('715599') ...
                };
        end
        
        function z = MiceColormap
            % Colormap for the mice
            z = {...
                'Purple', [202   152   227]/255, ...
                'Red', [196    77    88]/255, ...
                'Orange', [253   189    51]/255, ...
                'Yellow', [253   189    51]/255, ...
                'Green', [179   224   110]/255, ...
                'Blue', [0   177   229]/255 ...
                'White', [236 229 206]/255 ...
                };
        end
        
        function c = MiceColors(colors)
            % cmap = MiceColors(colors) Generates a colormap for the mice 
            % by a specific scheme. For example set 'colors' to 'PRBY' for
            % mice which are: purple, red, blue and yellow
            z = CheeseSquare.MiceColormap;
            c = zeros(length(colors), 3);
            for i=1:length(colors)
                idx = find(cellfun(@(x) ~isempty(x), regexp(z(1:2:end), colors(i))), 1);
                c(i, :) = z{2 * (idx-1) + 2};
            end
        end
        
        function s = GeneratePrefix(varargin)
            % Generate a file prefix (such as 'SC.exp0001.day01.cam02') from
            % meta data
            if isstruct(varargin{1})
                GroupType = varargin{1}.GroupType;
                GroupId = varargin{1}.GroupId;
                DayId = varargin{1}.DayId;
                CameraId = varargin{1}.CameraId;
            else
                 if nargin >= 1; GroupType =    varargin{1}; end
                 if nargin >= 2; GroupId =      varargin{2}; end
                 if nargin >= 3; DayId =        varargin{3}; end
                 if nargin >= 4; CameraId =     varargin{4}; end
            end
            if isempty(DayId)
                s = sprintf('%s.exp%04d.day%%02d.cam%02d', GroupType, GroupId, CameraId);
                return;
            end
            s = sprintf('%s.exp%04d.day%02d.cam%02d', GroupType, GroupId, DayId, CameraId);
        end
        
        function label = Ethogram(cheese, varargin)
            % label = Ethogram(cheese, options) Create both location and
            % behavioral ethograms and placed them in 'label'. For a list 
            % of valid 'options' see 'SocialEthogram' and 'LocationEthogram' 
            label = struct();
            [~, label.social] = CheeseSquare.SocialEthogram(cheese, varargin{:});
            [~, label.location] = CheeseSquare.LocationEthogram(cheese, varargin{:}, 'UseCategories', true);
        end
        
        function SocialEthogram(cheese, varargin)
            % SocialEthogram Display behavioral ethograms
            %
            %   SocialEthogram(options) Plots the behavioral ethogram for 
            %   all mice. Options can be
            %       frames - which frames to use
            %       size - size of created image
            %       brightness - increase brightness
            %% parse arguments
            p = inputParser;
            addParameter(p, 'frames', []);
            addParameter(p, 'size', [200 1000]);
            addParameter(p, 'brightness', 0.2);
            p.parse(varargin{:});
            opt = p.Results;
            %%
            if isempty(opt.frames)
                opt.frames = 1:size(cheese.Tracking.zones, 2);
            end
            %% interactions
            duration = length(opt.frames) / cheese.Video.FrameRate;
            Zones.Nest = find(cellfun(@(x) ~isempty(x), strfind(cheese.ROI.ZoneNames, '(BigNest)')));
            width  = floor(duration * cheese.Video.FrameRate / opt.size(2));
            height = opt.size(1) / cheese.nSubjects;
            
            background = zeros([cheese.nSubjects,  opt.size(2)]);
            for subj=1:cheese.nSubjects
                nest  = blockproc(ismember(cheese.Tracking.zones(subj, :), Zones.Nest), [1, width], @(x) mean(x.data));
                valid = blockproc(cheese.Tracking.valid, [1, width], @(x) mean(x.data)) > .5;
                nest = [nest(1:min(opt.size(2), length(nest))) nan(1, opt.size(2) - length(nest))];
                nest(~valid(1:min(length(nest), length(valid)))) = nan;
                background(subj, :) = nest;
            end
            background = repelem(background, height, 1);
            background = repmat((1-opt.brightness) + background*opt.brightness, [1 1 3]);
            background(isnan(background)) = 0;
            imagesc(background)
            hold on
            %
            cmap = [...
                Colors.Black; 
                Colors.RetroGreen; 
                Colors.RetroYellow; 
                Colors.RetroBlue; 
                Colors.RetroRed; 
                Colors.RetroOrange; 
                ];
            for subj=1:cheese.nSubjects
                map = any([cheese.Hierarchy.Contacts.List.subjects] == subj);
                list = cheese.Hierarchy.Contacts.List(map);
                %%
                for i=1:length(list)
                    curr = list(i);
                    subidx = curr.subjects == subj;
                    range = curr.beg(subidx):curr.end(subidx);
                    idx = floor((range - 1) / width) + 1;
                    for j=[4 2 3 5 6]
                        beg = min(idx(curr.data{subidx} == j));
                        if ~isempty(beg)
                            Patches.Rect(beg, height * (subj - 1), max(max(idx(curr.data{subidx} == j)) - beg, 1), height, cmap(j, :))
                        end
                    end
                    % being approached
%                     Map.appraoch(idx(curr.data{subidx} == 2)) = true;
%                     Map.beingappraoched(idx(curr.data{subidx} == 3)) = true;
%                     Map.contact(idx(curr.data{subidx} == 4)) = true;
%                     Map.chase(idx(curr.data{subidx} == 5 & agg(i))) = true;
%                     Map.escape(idx(curr.data{subidx} == 6 & agg(i))) = true;
                end
                %%
            end
            axis off
        end        

        function [im, label] = SocialEthogramOld(cheese, varargin)
            % SocialEthogramOld Old version of social ethogram
            %
            %   [im, label] = SocialEthogramOld(options) returns the
            %   created ethogram as an image ('im') and label matrix.
            %   Options can be:
            %       frames - which frames to use            
            %% parse arguments
            p = inputParser;
            addParameter(p, 'frames', []);
            p.parse(varargin{:});
            opt = p.Results;
            %%
            if isempty(opt.frames)
                opt.frames = 1:size(cheese.Tracking.zones, 2);
            end
            showmouse = false;
            %% interactions
            Zones.Nest = find(cellfun(@(x) ~isempty(x), strfind(cheese.ROI.ZoneNames, '(BigNest)')));
            label = zeros(cheese.nSubjects, size(cheese.Tracking.x, 2));
            for subj=1:cheese.nSubjects
                map = any([cheese.Hierarchy.Contacts.List.subjects] == subj);
                list = cheese.Hierarchy.Contacts.List(map);
                %
                agg = cheese.Hierarchy.Contacts.Behaviors.AggressiveChase.Map;
                agg = [agg(map(1:length(agg))) false(1,length(map) - length(agg))];
                %
                Map.nest = ismember(cheese.Tracking.zones(subj, :), Zones.Nest);                
                Map.appraoch = false(1, size(label, 2));
                Map.beingappraoched = false(1, size(label, 2));
                Map.chase  = false(1, size(label, 2));
                Map.escape = false(1, size(label, 2));
                Map.contact = false(1, size(label, 2));
                for i=1:length(list)
                    curr = list(i);
                    subidx = curr.subjects == subj;
                    range = curr.beg(subidx):curr.end(subidx);
                    Map.appraoch(range(curr.data{subidx} == 2)) = true;
                    Map.beingappraoched(range(curr.data{subidx} == 3)) = true;
                    Map.contact(range(curr.data{subidx} == 4)) = true;
                    Map.chase(range(curr.data{subidx} == 5 & agg(i))) = true;
                    Map.escape(range(curr.data{subidx} == 6 & agg(i))) = true;
                end
                %%
                label(subj, Map.nest) = 1;
                label(subj, Map.appraoch) = 2;
                label(subj, Map.beingappraoched) = 3;
                label(subj, Map.contact) = 4;
                label(subj, Map.chase) = 5;
                label(subj, Map.escape) = 6;
            end
            label = label(:, opt.frames);
            %%
            cmap = [...
                ones(1,3)*.95;
                Colors.White; 
                Colors.RetroGreen; 
                Colors.RetroYellow; 
                Colors.RetroBlue; 
                Colors.RetroRed; 
                Colors.RetroOrange; 
                ];
            if showmouse
                subplot(1,10,2:10);
            end
            im = ind2rgb(label+1, cmap);
            if nargout > 0
                return;
            end
            imagesc(im);
            axis off
            Fig.Fix
            colormap(cmap);
            nc = size(cmap, 1);
            %set(h, 'TickLabels', {'Unknown' obj.ROI.ZoneNames{:}});
            if showmouse
                subplot(1,10,1);
                m = obj.Colors.Marks;
                for i=1:obj.nSubjects
                    plot(0, i, 'o', 'MarkerFaceColor', mean(m.color(m.id == i, :)/255), 'MarkerEdgeColor', 'none', 'MarkerSize', 20);
                    Fig.Hon;
                end
                Fig.Hoff;
                Fig.YAxis(0.5, obj.nSubjects+.5);
                set(gca, 'YDir', 'reverse')
                axis off;
            end
        end        
                
        function [im, h]  = LocationHistogram(obj, varargin)
            % LocationHistogram Compute location histograms
            %
            %   [im, h]  = LocationHistogram(obj, options) Create a
            %   histogram of the mice locations and return it as an image
            %   ('im') or an histogram ('h'). See 'LocationEthogram' for a
            %   list of available options.
            %%
            [~, zlabel] = CheeseSquare.LocationEthogram(obj, varargin{:});
            label = zlabel - 1;
            idx = CheeseSquare.ZoneIdx(obj.ROI.ZoneNames);
            label(label ~= 0) = idx(label(label ~= 0));
            h = histc(label', 1:CheeseSquare.ZoneIdx);
            h = bsxfun(@rdivide, h, sum(h));
            %%
            n = 100;
            ch = round(n * cumsum(h));
            im = ones(n, obj.nSubjects);
            for s=1:obj.nSubjects
                for i=size(ch, 1):-1:1
                    im(1:ch(i, s), s) = i;
                end
            end
        end
            
        function [im, label] = LocationEthogram(obj, varargin)
            % SocialEthogram Compute location ethograms
            %
            %   [im, label] = LocationEthogram(options) Computes the location 
            %   ethograms for all mice. The ethograms are returned as an
            %   image ('im') or a matrix of all locations ('label')
            %   Available options are
            %       HighPlaces - list of elevatd areas listed as a cell
            %       FoodWater - list of feeding areas listed as a cell            
            %       UseCategories - set to true if should use categories
            %       rather than locations themselves (several ROIs might
            %       belong to the same category. See 'ZoneIdx')
            %% parse arguments
            p = inputParser;
            %addParameter(p, 'frames', []);
            p.addOptional('HighPlaces', {'BigNest', 'Block', '[Ramp1]', '[Ramp2]'});
            p.addOptional('FoodWater', {'Feeder1', 'Feeder2', 'Water', 'Water2'});
            
            addParameter(p, 'UseCategories', false);
            p.parse(varargin{:});
            opt = p.Results;
            
            %%
            showmouse = false;
            zmap = CheeseSquare.ZoneColormap;
            map = containers.Map(zmap(1:2:end), zmap(2:2:end));
            zones = obj.Tracking.zones;
            label = ones(size(zones ));
            cmap = zeros(obj.ROI.nZones + 1, 3);
            idx = 2;
            for z=1:obj.ROI.nZones
                %%
                if ~map.isKey(obj.ROI.ZoneNames{z})
                    warning('no color defined for ''%s''', obj.ROI.ZoneNames{z});
                else
                    label(zones == z) = idx;
                    cmap(idx, :) = map(obj.ROI.ZoneNames{z});
                    idx = idx + 1;
                end
            end
            %%
            if showmouse
                subplot(1,10,2:10);
            end
            if opt.UseCategories
                cmap = lines(3);
                l = zeros(size(label));
                %%
                l(~obj.Tracking.sheltered) = 1;
                %%
                HighPlacesNames = opt.HighPlaces;
                HighPlaces = [];
                for i=1:length(HighPlacesNames)
                    idx = find(strcmp(obj.ROI.ZoneNames, HighPlacesNames{i}));
                    if ~isempty(idx)
                        HighPlaces = [HighPlaces, idx];
                    end
                end
                l(ismember(obj.Tracking.zones, HighPlaces)) = 3;
                %%
                FoodWaterNames = opt.FoodWater;
                FoodWater = [];
                for i=1:length(FoodWaterNames)
                    idx = find(strcmp(obj.ROI.ZoneNames, FoodWaterNames{i}));
                    if ~isempty(idx)
                        FoodWater = [FoodWater, idx];
                    end
                end
                l(ismember(obj.Tracking.zones, FoodWater)) = 2;
                %%
                l(:, ~obj.Tracking.valid) = nan;
                %%
                label = l;
                
            end
            
            im = ind2rgb(label, cmap);
            if nargout > 0
                return;
            end
            imagesc(im);
            axis off
            Fig.Fix
            colormap(cmap);
            h = colorbar;
            nc = size(cmap, 1);
            set(h, 'Ticks', linspace(0+.5/nc, 1-.5/nc, nc));
            set(h, 'TickLabels', {'Unknown' obj.ROI.ZoneNames{:}});
            if showmouse
                subplot(1,10,1);
                m = obj.Colors.Marks;
                for i=1:obj.nSubjects
                    plot(0, i, 'o', 'MarkerFaceColor', mean(m.color(m.id == i, :)/255), 'MarkerEdgeColor', 'none', 'MarkerSize', 20);
                    Fig.Hon;
                end
                Fig.Hoff;
                Fig.YAxis(0.5, obj.nSubjects+.5);
                set(gca, 'YDir', 'reverse')
                axis off;
            end
        end

        function out = ProfileManager(name, value)
            % ProfileManager Keeps profiles in memory (to avoid computing
            % them repetitively)
            %
            % Note! Don't use this function directly
            persistent manager;
            if isempty(manager)
                manager = containers.Map();
            end
            if nargin > 1;
                manager(name) = value;
            elseif ~manager.isKey(name)
                manager(name) = [];
            end
            out = manager(name); 
        end
        
    end
    
    methods
        function s = Sheltered(cheese, varargin)
            % Sheltered(begf, endf) Returns a matrix of size nxF (where n
            % is number of mice and F number of frames) which indicates in
            % each frame and for each mouse whether he was inside a
            % shelter
            s = cheese.GetTrackingData('sheltered', varargin{:});
        end

        function h = Hidden(cheese, varargin)
            % Hidden(begf, endf) Returns a matrix of size nxF (where n
            % is number of mice and F number of frames) which indicates in
            % each frame and for each mouse whether he was hidden from the
            % camera.
            h = cheese.GetTrackingData('hidden', varargin{:});
        end
        
        function [x, y] = Position(cheese, varargin)
            % [x,y] = Position(begf, endf) Returns the matrices of size nxF (where n
            % is number of mice and F number of frames) that hold the x and
            % y coordinates of the mice.
            x = cheese.GetTrackingData('x', varargin{:});
            y = cheese.GetTrackingData('y', varargin{:});
        end

        function z = Zone(cheese, varargin)
            % Zone(begf, endf) Returns a matrix of size nxF (where n
            % is number of mice and F number of frames) which indicates for
            % each frame the region that each mouse is in.
            z = cheese.GetTrackingData('zones', varargin{:});
        end
        
        function out = SetZone(cheese, subj, zone, frames)
            % SetZone Manually fixing tracking
            %
            %   cheese = SetZone(cheese, subj, zone, frames) Set the location
            %   of mouse 'subj' as 'zone' in the frames specified in
            %   'frames'
            %
            %%
            centers = zeros(cheese.ROI.nRegions, 2);
            idx = 2;
            zid = zeros(1, cheese.ROI.nZones);
            for i=1:cheese.ROI.nRegions
                [x,y] = find(cheese.ROI.Regions{i});
                centers(i, :) = [mean(x), mean(y)];
                zid(idx) = i; idx = idx + 1;
                if cheese.ROI.IsSheltered(i)
                    zid(idx) = i; idx = idx + 1;
                end
            end
            zid = zid + 1;
            centers = [size(cheese.ROI.Regions{1})/2; centers];
            %%
            out = cheese;
            out.Tracking.zones(subj, frames) = zone;
            out.Tracking.x(subj, frames) = centers(zid(zone), 2);
            out.Tracking.y(subj, frames) = centers(zid(zone), 1);
            out.Tracking.sheltered(subj, frames) = cheese.ROI.IsSheltered(zid(zone));
            %out.Tracking.sheltered(subj, frames) = cheese.ROI.IsSheltered(zone;
        end
        
        function value = GetTrackingData(cheese, field, begdt, enddt)
            % GetTrackingData Get tracking data
            %
            % value = GetTrackingData(cheese, field, begdt, enddt) return
            % tracking data of type 'field' during the frames between
            % 'begdt' and 'enddt'
            %
            if nargin >= 4
                idx = round(DateTime.ToFrames(begdt - cheese.DateTime, cheese.Video.FrameRate) + 1) : round(DateTime.ToFrames(enddt - cheese.DateTime, cheese.Video.FrameRate) + 1);
            else
                idx = round(DateTime.ToFrames(begdt - cheese.DateTime, cheese.Video.FrameRate) + 1);
            end
            valid = idx >= 1 & idx <= length(cheese.Tracking.valid);
            value = nan(length(valid), size(cheese.Tracking.(field), 1));
            value(valid, :) = cheese.Tracking.(field)(:, idx(valid))';
        end
        
        function obj = CheeseSquare(varargin)
            % Constructor
            %
            %   obj = CheeseSquare(cheese) initialize from existing 
            %   CheeseSquare object 'cheese'
            %
            %   obj = CheeseSquare(videofilename) initialize from a video
            %   file
            %
            %   obj = CheeseSquare(cheesefilename) initialize from a
            %   CheeseSquare file
            %
            %   obj = CheeseSquare(oldmat) initialize from an old version
            %   CheeseSquare file
            %
            obj = obj.Load(varargin{:});
        end

        function delete(obj)
            % destructor
            CheeseSquare.ProfileManager(obj.UID, []);
        end
        
        function obj = Load(obj, varargin)
            % Load data from file
            %
            %   obj = Load(obj, cheese) initialize from existing 
            %   CheeseSquare object 'cheese'
            %
            %   obj = Load(obj, videofilename) initialize from a video
            %   file
            %
            %   obj = Load(obj, cheesefilename) initialize from a
            %   CheeseSquare file
            %
            %   obj = Load(obj, oldmat) initialize from an old version
            %   CheeseSquare file
            %
            if nargin >= 1
                arg = varargin{1};
                if isa(arg, 'CheeseSquare')
                    obj = arg;
                elseif isa(arg, 'char')
                    filename = MyFilename(varargin{1});
                    if strcmp(filename.Ext, '.mat')
                        classfield = fieldnames(obj);
                        warning('off', 'MATLAB:load:variableNotFound');
                        %%
                        m = load(varargin{1});
                        if isfield(m, 'social') % compressed format
                            obj.Tracking = m.social;
                            src = Q.cpfield(m.social, struct(), classfield, true);
                            warning('on', 'MATLAB:load:variableNotFound');
                            loaded = fieldnames(src);
                            for i=1:length(loaded)
                                obj.(loaded{i}) = src.(loaded{i});
                            end
                            obj.Filename = filename;
                        elseif isfield(m, 'Tracking') % CheeseSquare format
                            src = struct();
                            src = Q.cpfield(m, src, obj.TrackingFields, true);
                            obj.Tracking = src;
                            src = Q.cpfield(m, src, classfield, true);
                            warning('on', 'MATLAB:load:variableNotFound');
                            loaded = fieldnames(src);
                            for i=1:length(loaded)
                                obj.(loaded{i}) = src.(loaded{i});
                            end
                            obj.Filename = filename;
                        else % old format
                            %%
                            obj.Prefix = m.FilePrefix;
                            obj.Tracking = Q.cpfield(m, struct(), {'x', 'y', 'zones', 'hidden', 'sheltered', 'valid', 'time'}, true);
                            obj.Video = m.VideoFile;
                            obj.Filename = filename;
                            obj = Q.cpfield(m, obj, {'Hierarchy'}, true);
                            obj.Hierarchy = Q.cpfield(m, obj.Hierarchy, {'Contacts'}, true);
                            obj = Q.cpfield(m, obj, {'Background'}, true);
                            obj = Q.cpfield(m, obj, {'nSubjects'}, true);
                            obj = Q.cpfield(m, obj, {'Colors'}, true);
                            if isfield(m, 'PixelsPerCM')
                               obj.Meta.Scale.PixPerCM = m.PixelsPerCM;
                            end
                            %src = Q.cpfield(m, src, {'Analysis'}, true);
                            obj.Analysis = struct();
                            obj = Q.cpfield(m, obj, {'ROI'}, true);
                            obj.Profile_ = [];
                            %%
                            obj = CheeseSquare(obj);
                            %% if video file is missing, try to complete missing information
                            if obj.Video.FrameRate == 0;        obj.Video.FrameRate = m.FrameRate;    end
                            if obj.Video.NumberOfFrames == 0;   obj.Video.NumberOfFrames = m.nFrames; end
                            if obj.Video.Height == 0;           obj.Video.Height = m.VideoHeight;     end
                            if obj.Video.Width == 0;            obj.Video.Width = m.VideoWidth;       end
                            if obj.Video.Duration == 0;         obj.Video.Width = m.nFrames / m.FrameRate; end
                        end
                        
                        %%
%                         m = matfile(varargin{1});
%                         if isprop(m, 'social')
%                             obj.Tracking = m.social;
%                             src = Q.cpfield(m.social, struct(), classfield, true);
%                         else
%                             src = load(varargin{1}, obj.TrackingFields{:});
%                             obj.Tracking = src;
%                             src = load(varargin{1}, classfield{:});                            
%                         end
                    else % .avi
                        objmat = [filename.Full, '.obj.mat'];
                        if exist(objmat, 'file')
                            obj = obj.Load(objmat, filename.Full, varargin{2:end});
                        else
                            obj.Video = varargin{1};
                        end
                        obj.Filename = MyFilename(objmat);
                    end
                elseif isa(arg, 'CheeseSquare')
                    obj = varargin{1};
                elseif isa(arg, 'struct')
                    obj.Tracking = arg;
                    for i=1:length(obj.TrackingFields)
                        if isfield(arg, obj.TrackingFields{i})
                            %                             %obj.addprop(obj.TrackingFields{i});
                            obj.(obj.TrackingFields{i}) = arg.(obj.TrackingFields{i});
                        end
                    end
                    obj.Filename = obj.Tracking.Source;
                    obj.Video = obj.Tracking.VideoFile;
                end
                obj.Filename = MyFilename(obj.Filename);
                obj.Prefix = obj.GetPrefix;
                if length(varargin) > 1
                    obj.Video = varargin{2};
                end
                if ischar(obj.Video)
                    alt = MyFilename(obj.Video).SetPath(obj.Filename.Path).Full;
                    if exist(obj.Video, 'file')
                        obj.Video = MyVideoReader(obj.Video);
                    elseif exist(alt, 'file')
                        obj.Video = MyVideoReader(alt);
                    else
                        obj.Video = MyVideoReader(alt);
                    end
                else
                    alt = [MyFilename(obj.Prefix).SetPath(obj.Filename.Path).Full '.avi'];
                    if isempty(obj.Video)
                        obj.Video = MyVideoReader(alt);
                    elseif isa(obj.Video, 'MyVideoReader')
                        obj.Video = MyVideoReader(obj.Video);
                    end
                end
                obj.Meta = Q.cpfield(obj.LoadMetaData, obj.Meta);
                if (isempty(obj.nSubjects) || obj.nSubjects == 0) && Q.isfield(obj, 'Colors.Histrogram.Count')
                    obj.nSubjects = length(obj.Colors.Histrogram.Count);
                end
            end
            obj.Profile_ = [];
            obj = obj.Fix();
        end
        
        function obj = Fix(obj)
            % Check file for errors and try to fix them
            if obj.Video.FrameRate == 0
                % try to fix files where the frame rate is set to 0
                fps = 0;
                try
                    fps = median(1./diff(obj.Tracking.time));
                catch
                end
                if fps == 0
                    warning('cannot determine framerate; assuming 25 fps');
                    fps = 25;
                end
                obj.Video.FrameRate = fps;
            end

        end
        
        function meta = LoadMetaData(obj)
            % Load additional data from xml file that comes with avi file
            % (no longer used)
            meta = struct();
            try
                meta = CheeseSquare.ParseName(obj.Video.Filename);
            catch
            end
            %%
            try
                videofile = MyFilename(obj.Video.Filename);
                videofile = videofile.SetExt([videofile.Ext, '.xml']);
                meta.CreationDate = addtodate(MyFilename(obj.Video.Filename).CreationDate, -round(obj.Video.NumberOfFrames / obj.Video.FrameRate*1000), 'millisecond');
                meta.Remarks = '';
                if exist(videofile.Full, 'file')
                    metaxml = xmlread(videofile.Full);
                    cdate = char(metaxml.getElementsByTagName('Start').item(0).getTextContent);
                    meta.CreationDate = datenum(cdate, 'yyyy-mm-ddTHH:MM:SS.FFF');
                    meta.Remarks =  char(metaxml.getElementsByTagName('Remarks').item(0).getTextContent);
                end
            catch
            end
        end
        
        function track = ToClassical(obj)
            % Convert CheeseSqure object to classical structure. Needed for
            % older code segments that don't support CheeseSquare.
            if ischar(obj.Video)
                track = CheeseObject.Create(obj.Video);
            else
                track = CheeseObject.Create(obj.Video.Filename);
            end
            track.nSubjects = obj.GetField('nSubjects');
            track.Colors = obj.GetField('Colors');
            track.Colors.Centers = CheeseSquare.MiceColors('PROG');
            track.BkgImage = obj.GetField('Background.im');
            track.ROI = obj.GetField('ROI');
            track.Hierarchy = obj.GetField('Hierarchy');
            
            %% try to fix possible issues (using redundencies)
            if track.nSubjects == 0
                try
                    track.nSubjects = size(obj.Colors.Centers, 1);
                catch
                    try
                        track.nSubjects = size(obj.Tracking.x, 1);
                    catch
                    end
                end
            end
            
            if track.nFrames == 0
                if obj.Video.NumberOfFrames == 0
                    try
                       track.nFrames = size(obj.Tracking.x, 2);
                    catch
                    end
                else
                    track.nFrames = obj.Video.NumberOfFrames;
                end
            end
            
            %%
            track.time = linspace(0, obj.Video.Duration, track.nFrames);
            
            try
                track = Q.cpfield(obj.Meta, track, 'Scale');
                track.PixelsPerCM = mean([obj.Meta.Scale.ArenaCoord(3) / obj.Meta.Scale.ArenaWidth, obj.Meta.Scale.ArenaCoord(4) / obj.Meta.Scale.ArenaHeight]);
            catch
                track.PixelsPerCM = [];
            end
            try
                
            catch
            end
            
            if isstruct(obj.Tracking)
                f = fieldnames(obj.Tracking);
                for i=1:length(f)
                    if ~strcmpi(f{i}, 'Properties')
                        track.(f{i}) = obj.Tracking.(f{i});
                    end
                end
            end
        end
        
        function obj = SaveOldFormat(obj)
            % Save CheeseSquare object in the older format (for back
            % compatibility)
            if isstruct(obj.Tracking)
                f = fieldnames(obj.Tracking);
                track = struct();
                for i=1:length(f)
                    if ~strcmpi(f{i}, 'Properties')
                        track.(f{i}) = obj.Tracking.(f{i});
                    end
                end
            end
            track.Background = obj.Background;
            track.Colors = obj.Colors;
            track.ROI = obj.ROI;
            track.OutputPath = [obj.Filename.Path filesep];
            track.FilePrefix = regexprep(obj.Filename.Name, '\.obj', '', 'ignorecase');
            track.VideoFile = obj.GetField('Video.Filename');
            TrackSave(track);
        end
        
        function obj = ExportProfile(obj)
            % Write the profile to an external file. The filename is
            % similar to the video file with a '.profile.mat' extension
            %%
            profile = struct();
            %profile.x            = obj.GetField('Tracking.x');
            %profile.y            = obj.GetField('Tracking.y');
            %profile.Hierarchy    = obj.GetField('Hierarchy');
            %profile.Profile      = obj.GetField('Profile', 0);
            profile.Background   = im2uint8(obj.GetField('Background.im'));
            profile.Regions      = im2uint8(obj.GetField('Background.im'));
            profile.IsTracked    = obj.IsTracked();
            profile.IsReady      = obj.IsReady();
            profile.Remarks      = obj.GetField('Meta.Remarks', '');
            profile.CreationDate = datestr(obj.GetField('Meta.CreationDate'), 'yyyy-mm-ddTHH:MM:SS.FFF');
            profile.Duration     = obj.GetField('Video.Duration');
            profile.nSubjects    = obj.GetField('nSubjects', 0);
            profile.Colors       = [];
            save([obj.Video.Filename '.profile.mat'], '-struct', 'profile', '-v6');
            Console.Message('writing profile to ''%s''', regexprep([obj.Video.Filename '.profile.mat'], '\', '/'));
        end
        
        function obj = Save(obj)
            % Save CheeseSquare to file
            out = obj.SerializeOut(); %#ok<NASGU>
            save(obj.Filename.Full, '-struct', 'out', '-v7.3');
            try
                obj.ExportProfile();
            catch
            end
        end
        
        function y = IsTracked(obj)
            % True if file has been tracked
            y = true;
            y = y && EnsureField(obj, 'Tracking.x');
            y = y && EnsureField(obj, 'Tracking.y');
            y = y && EnsureField(obj, 'Tracking.zones');
        end
        
        function y = IsReady(obj)
            % True if file had been marked for tracking (and ready to be
            % tracked)
            y = true;
            y = y && EnsureField(obj, 'Background.im');
            y = y && EnsureField(obj, 'Colors.Histrogram.Count');
            y = y && EnsureField(obj, 'ROI.Regions');
        end
       
        
        function Show(cheese, varargin)
            % Show Generate a plot of the group behavior
            %
            %   Show(cheese, options) - the options can be:
            %       ShowVideo - should use and show frames from video file
            %       ShowROI - should overlay ROIs
            %       DateTime - show data for specific times
            %
            opt = Q.defaultargs(false, varargin ...
                , 'ShowVideo', true ...
                , 'ShowROI', true ...
                , 'DateTime', [] ...
                );
            
            if opt.ShowVideo
                imagesc(cheese.Video.CurrentFrame);
            else
            end
            %
            if opt.ShowROI
                if opt.ShowVideo
                    cmap = lines;
                else
                    cmap = ones(50, 3)*.7;
                end
                for r=find(cheese.ROI.IsAvail)
                    h = patch(cheese.ROI.Points{r}(:, 1), cheese.ROI.Points{r}(:, 2), cmap(r, :), 'facealpha', .3, 'edgecolor', 'none');
                end
            end
            %
            %axis off;
            set(gca, 'XTick', [], 'YTick', [], 'XColor', [1 1 1]*.2, 'YColor', [1 1 1]*.2);
            axis equal;
            box on;
            w = 15;
            if Q.isfield(cheese, 'Tracking.x')
                Fig.Hon
                if isempty(opt.DateTime)
                    frame = cheese.Video.FrameNumber;
                else
                    frame = round(DateTime.ToSecs(opt.DateTime - cheese.Meta.CreationDate) * cheese.Video.FrameRate);
                end
                if frame > 0 && frame <= cheese.Video.NumberOfFrames
                    for s=1:cheese.nSubjects
                        x = double(cheese.Tracking.x(s, frame));
                        y = double(cheese.Tracking.y(s, frame));
                        color = uint8(mean(cheese.Colors.Marks.color(cheese.Colors.Marks.id == s, :)));
                        patch([x-w x+w x+w x-w], [y-w y-w y+w y+w], 'w', 'edgecolor', 'w', 'facecolor', color);
                        text(x+w, y+w, num2str(s), 'color', color, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'fontsize', w);
                    end
                end
                Fig.Hoff
            end
            xlim([0, size(cheese.Background.im, 2)]);
            ylim([0, size(cheese.Background.im, 1)]);
        end
        
        function Tracjectories(cheese, varargin)
            % Plot mouse trajectories
            if ~isempty(varargin) > 0 && isnumeric(varargin{1})
                range = varargin{1};
                varargin = varargin(2:end);
            else
                range = [1, size(cheese.Tracking.x, 2)];
            end
            %% parse arguments
            p = inputParser;
            addParameter(p, 'Colors', 'PRBY', @isstr);
            p.parse(varargin{:});
            opt = p.Results;
            %%
            im = cheese.Background.im;  
            cmap = CheeseSquare.MiceColors(opt.Colors); 
            imagesc(im);
            hold on;
            for i=1:cheese.nSubjects
                scatter(cheese.Tracking.x(i, range(1):range(2)), cheese.Tracking.y(i, range(1):range(2)), 'MarkerFaceColor', cmap(i, :), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); 
            end
            hold off;
        end
            
    end
    
    %% Getter and setter functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function value = get.Profile(obj)
            value = CheeseSquare.ProfileManager(obj.UID);
        end
        
        function obj = set.Profile(obj, value)
            if isempty(value)
                [~, profile] = CheeseNewProfiler(obj);
                value = CheeseSquare.ProfileManager(obj.UID, profile);
            end
            CheeseSquare.ProfileManager(obj.UID, value);
        end

        function value = get.Profile_(obj)
            value = obj.Profile_;
        end
        
        function obj = set.Profile_(obj, value)
            obj.Profile_ = value;
        end
        
        function value = get.DateTime(obj)
            value = obj.Meta.CreationDate + DateTime.FromSecs(obj.Video.Time);
        end
        
        function obj = set.DateTime(obj, value)
            if ~ischar(obj.Video)
                obj.Video.Time = (value - obj.Meta.CreationDate) * 24 * 60 * 60;
            end
        end
    end
    
    %% Static aux functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static = true)
        function Rename(filename, prefix)
            % Used to rename an object file. Don't use the system's rename
            % function because some information is stored inside the
            % files!!
            %
            %   Rename(filename, prefix) Rename object in 'filename' to
            %   have a new prefix 'prefix'
           try
                source = CheeseSquare.ParseName(filename);
            catch me
                Console.Warning(me, 'failed to parse source name');
                return;
           end
            %%
            try
                target = CheeseSquare.ParseName(prefix);
            catch me
                Console.Warning(me, 'failed to parse target name');
                return;
            end
            %%
            [folder, ~, ~] = fileparts(filename);
            sourcePrefix = CheeseSquare.GeneratePrefix(source);
            targetPrefix = CheeseSquare.GeneratePrefix(target);
            videoFile = fullfile(folder, [sourcePrefix, '.avi']);
            profileFile = fullfile(folder, [sourcePrefix, '.avi.profile.mat']);
            cheeseFile = fullfile(folder, [sourcePrefix, '.avi.obj.mat']);
            todel = {};
            
            %% video file
            if exist(videoFile, 'file')
                targetVideo = fullfile(folder, [targetPrefix, '.avi']);
                if exist(targetVideo, 'file')
                    error('target file %s already exists; exiting', targetVideo);
                end
                Console.Write('writing video file to -> %s...', regexprep(targetVideo, '\\', '\\\\'));                
                try
                    copyfile(videoFile, targetVideo);
                catch me
                    Console.Failed
                    Console.Warning(me);
                    return;
                end
                todel{end+1} = videoFile;
                Console.Done;
            else
                Console.Message('could not find video file [skip]');                
            end
            %% Cheese square
            if exist(cheeseFile, 'file')
                obj = load(cheeseFile);
                obj.Prefix = targetPrefix;
                obj.Video = videoFile;
                targetCheese = fullfile(folder, [targetPrefix, '.avi.obj.mat']);
                if exist(targetCheese, 'file')
                    error('target file %s already exists', targetCheese);
                end
                Console.Write('writing CheeseSquare object file to -> %s...', regexprep(targetCheese, '\\', '\\\\'));                
                try
                    save(targetCheese, '-struct', 'obj');
                catch me
                    Console.Failed
                    Console.Warning(me);
                    return;
                end
                todel{end+1} = cheeseFile;
                Console.Done;
            else
                Console.Message('could not find CheeseSquare object file [skip]');
            end
            %% Profile file
            if exist(profileFile, 'file')
                targetProfile = fullfile(folder, [targetPrefix, '.avi.profile.mat']);
                if exist(targetProfile, 'file')
                    error('target profile file %s already exists; exiting', targetProfile);
                end
                Console.Write('writing profile file to -> %s...', regexprep(targetProfile, '\\', '\\\\'));                
                try
                    copyfile(profileFile, targetProfile);
                catch me
                    Console.Failed
                    Console.Warning(me);
                    return;
                end
                todel{end+1} = profileFile;
                Console.Done;
            else
                Console.Message('could not find profile file [skip]');                
            end
            %%
            todelstr = regexprep(todel, '\\', '\\\\');
            str = sprintf(['Are you sure you want to delete these files?' '\n    ' strjoin(todelstr, '\n    ')]);
            choice = questdlg(str, 'Delete files', 'Yes', 'No', 'No');
            switch choice
                case 'Yes'
                    for i=1:length(todel)
                        Console.Write('deleting [x] %s...', todelstr{i});
                        delete(todel{i});
                        Console.Done;
                    end
                case 'No'
                    Console.WriteLine('not deleting source files');
            end
        end
        
        
        function meta = ParseName(name)
            % Parse prefix name and break it to 'GroupType', 'GroupID',
            % 'Day', 'Camera'. 
            %
            % For example, 'SC.exp0001.day02.cam04' is broken
            % to: 
            %   GroupType = SC
            %   GroupID = 1
            %   Day = 2
            %   Camera = 4
            expression = '(?<grouptype>\w+)\.exp(?<groupid>\d+)\.day(?<dayid>\d+)\.cam(?<camid>\d+)';
            res = regexp(name, expression, 'names');
            meta.GroupType = res.grouptype;
            meta.GroupId = str2double(res.groupid);
            meta.DayId = str2double(res.dayid);
            meta.CameraId = str2double(res.camid);
        end
        %%
        
        function tt = CreateProfileTable(filelist, varargin)
            % Create a single profile table from several files
            % Syntax: tt = CreateProfileTable(filelist, options) where
            % options can be:
            %   path - path to object files
            %   ext - file extension
            %   grouptype - set group type in the table
            %
            %% parse arguments
            p = inputParser;
            addParameter(p, 'path', '', @isstr);
            addParameter(p, 'ext', '', @isstr);
            addParameter(p, 'grouptype', '', @isstr);
            addParameter(p, 'usegroups', true, @islogical);
            p.parse(varargin{:});
            opt = p.Results;
            %%
            ischeesefarm = false;
            if opt.usegroups
                groups = containers.Map;
                nfiles = 0;
                for i=1:length(filelist)
                    if isstruct(filelist)
                        fname = filelist(i).name;
                        if isempty(opt.path)
                            opt.path = filelist(i).folder;
                        end
                    elseif iscell(filelist)
                        if ischar(filelist{i})
                            fname = filelist{i};
                        elseif isa(filelist{i}, 'CheeseSquare')
                            ischeesefarm = true;
                            fname = filelist{i}.Prefix;
                        end
                    else
                        error
                    end
                    try
                        p = CheeseSquare.ParseName(fname);
                        nfiles = nfiles + 1;
                    catch
                        warning(sprintf('failed parsing file named ''%s''', fname));
                    end
                    key = sprintf('%s.exp%04d.cam%04d', p.GroupType, p.GroupId, p.CameraId);
                    if groups.isKey(key)
                        curr = groups(key);
                    else
                        curr = {};
                    end
                    curr{p.DayId} = fname;
                    groups(key) = curr;
                end
                %%
                opt.path = [regexprep(opt.path, '\\', '/'), '/'];
                %%
                T = {};
                groupnames = groups.keys;
                failedfiles = {};
                for i=1:length(groupnames)
                    curr = groups(groupnames{i});
                    Console.Message('processing group ''%s''', groupnames{i});
                    cheese = cell(1, length(curr));
                    tgroup = cell(1, length(curr));
                    days = zeros(1, length(curr));
                    failed = false(1, length(curr));
                    parfor idx=1:length(curr) % parfor
                        filename = curr{idx};
                        if isempty(filename)
                            continue;
                        end
                        %%
                        try
                            p = [];
                            if ~ischeesefarm
                                Console.Message(1, 'loading ''%s''', [opt.path, filename opt.ext]);
                                obj = CheeseSquare([opt.path, filename opt.ext]);
                                obj.Profile = [];
                                p = obj.Profile;
                            else
                                Console.Message(1, 'processing ''%s''', filename);
                                obj = [];
                                for f=1:length(filelist)
                                    if strcmp(filelist{f}.Prefix, filename)
                                        obj = filelist{f};
                                        p = obj.Profile;
                                        break;
                                    end
                                end
                            end
%                             cheese{idx}.Prefix = filename;
%                             cheese{idx} = CheeseSquare([opt.path, filename opt.ext]); %#ok<PFBNS>
%                             cheese{idx}.Profile = [];
%                             p = cheese{idx}.Profile;
                            t = p.Tables.Individual;
                            %
                            s = cell(1, size(t, 1));
                            g = cell(1, size(t, 1));
                            for j=1:size(t, 1)
                                if ~isempty(opt.grouptype)
                                    grouptype = opt.grouptype;
                                else
                                    grouptype = obj.Meta.GroupType;
                                end
                                s{j} = sprintf('%s.exp%04d.day%02d.cam%02d.sub%02d', grouptype, obj.Meta.GroupId, obj.Meta.DayId, obj.Meta.CameraId, t.MouseID(j));
                                g{j} = sprintf('%s.exp%04d.cam%02d', grouptype, obj.Meta.GroupId, obj.Meta.CameraId);
                                t.GroupType{j} = grouptype;
                            end
                            t.Properties.RowNames = s;
                            tgroup{idx} = t;
                            cheese{idx} = obj;
                            days(idx) = cheese{idx}.Meta.DayId;
                        catch me
                            Console.Warning(me, 'Error in ''%s''', filename);
                            failedfiles{idx} = filename;
                            failed(idx) = true;
                        end
                    end
                    %%
                    if sum(failed) > 0
                        Console.Warning('failed on %d files', sum(failed));
                        %failedfiles = cat(2, failedfiles, cellfun(@(x) x.Prefix, cheese(failed), 'UniformOutput', false));
                        %cheese(failed) = [];
                    end
                    
                    %% arrange data by days
                    [days, o] = sort(days);
                    tgroup = tgroup(o);
                    cheese = cheese(o);
                    %%
                    tgroup(cellfun(@(x) isempty(x) || ~isprop(x, 'Tracking') || isempty(x.Tracking), cheese)) = [];
                    cheese(cellfun(@(x) isempty(x) || ~isprop(x, 'Tracking') || isempty(x.Tracking), cheese)) = [];
                    %%
                    try
                        %%
                        Elo = Hierarchy.EloRating(cheese);
                        prevs = [];
                        for day=1:length(cheese)
                            %%
                            s = Elo.Segments(day).Scores;
                            %T{startidx+day-1} = [T{startidx+day-1}, array2table(s', 'VariableNames', {'TotalElo1', 'TotalElo2', 'TotalElo3', 'TotalElo4', 'TotalElo5', 'TotalElo6'})];
                            if isempty(s)
                                if isempty(prevs)
                                    s = zeros(6, cheese{day}.nSubjects);
                                else
                                    s = prevs;
                                end
                            end
                            s = s(1:6, :);
                            if size(s, 1) < 6
                                s(end+1:6, :) = nan;
                            end
                            tgroup{day} = [tgroup{day}, table(s', 'VariableNames', {'CumEloRating'})];
                            prevs = s;
                        end
                        for day=1:length(cheese)
                            tgroup{day} = [tgroup{day}, array2table(s(end, :)', 'VariableNames', {'TotalEloRating'})];
                        end
                    catch me
                        Console.Warning(me);
                    end
                    %%
                    try
                        Sij = 0;
                        for day=1:length(cheese)
                            curr = cheese{day};
                            Sij = Sij + curr.Hierarchy.AggressiveChase.ChaseEscape;
                            NormDS = Hierarchy.DavidScore(Sij);
                            tgroup{day} = [tgroup{day}, array2table(NormDS(:), 'VariableNames', {'CumNormDavidScore'})];
                        end
                        for day=1:length(cheese)
                            tgroup{day} = [tgroup{day}, array2table(NormDS(:), 'VariableNames', {'TotalNormDavidScore'})];
                        end
                    catch me
                        Console.Warning(me);
                    end
                    %%
                    T(end+1:end+length(tgroup)) = tgroup;
                end
            else
                %% load CheeseSquare files
                T = cell(1, length(filelist));
                for i=1:length(filelist) % parfor
                    if isstruct(filelist) % created using 'dir'?
                        filename = filelist(i).name;
                    else % should be cell
                        filename = filelist{i}
                    end
                    cheese = CheeseSquare([opt.path, filename]); %#ok<PFBNS>
                    %%
                    try
                        cheese.Profile = [];
                        p = cheese.Profile;
                        t = p.Tables.Individual;
                        %
                        s = cell(1, size(t, 1));
                        g = cell(1, size(t, 1));
                        for j=1:size(t, 1)
                            if ~isempty(opt.grouptype)
                                grouptype = opt.grouptype;
                            else
                                grouptype = cheese.Meta.GroupType;
                            end
                            s{j} = sprintf('%s.exp%04d.day%02d.cam%02d.sub%02d', grouptype, cheese.Meta.GroupId, cheese.Meta.DayId, cheese.Meta.CameraId, t.MouseID(j));
                            g{j} = sprintf('%s.exp%04d.cam%02d', grouptype, cheese.Meta.GroupId, cheese.Meta.CameraId);
                            t.GroupType{j} = grouptype;
                        end
                        t.Properties.RowNames = s;
                        T{i} = t;
                    catch me
                        Console.Warning(me, 'Error in ''%s''', filename);
                    end
                    
                end
            end
            %% create merged table
            tt = [];
            for i=1:length(T)
                if isempty(tt)
                    tt = T{i};
                else
                    curr = T{i};
                    if isempty(curr)
                        continue;
                    end
                    count = 1;
                    while any(ismember(curr.Properties.RowNames, tt.Properties.RowNames))
                        for j=1:size(T{i}, 1)
                            curr.Properties.RowNames{j} = [T{i}.Properties.RowNames{j} '_' num2str(count)];
                        end
                        count = count + 1;
                    end
                    tt = [tt; curr]; %#ok<AGROW>
                end
            end
            %% fix indices for 'GroupNumber' and 'MouseNumber'
            %[~,~,grouptype] = unique(tt.GroupType);
            [~, ~, groupnum] = unique(regexprep(regexprep(tt.Properties.RowNames, 'day[0-9]*', ''), 'sub[0-9]*', ''));
            tt.SourceGroupNumber = tt.GroupNumber;
            tt.SourceMouseNumber = tt.MouseNumber;
            tt.GroupNumber = groupnum;
            tt.MouseNumber = (groupnum - 1)*max(tt.MouseID)+tt.MouseID;
        end
        
        function nn = TableFingerprint(source, target, varargin)
            % Create a unique identifier for a table. Used to check if
            % tables are identical (same group same day)
            p = inputParser;
            addParameter(p, 'fields', {'ForagingCorrelation', 'FractionOfTimeOutside', 'VisitsOutsideRate', 'Entropy'});
            addParameter(p, 'thresh', 1e-6);
            p.parse(varargin{:});
            opt = p.Results;
            %%
            [~, tgtidx] = ismember(opt.fields, target.Properties.VariableNames);
            tgtm = table2array(target(:, tgtidx));
            [~, srcidx] = ismember(opt.fields, source.Properties.VariableNames);
            srcm = table2array(source(:, srcidx));
            [nn, D] = knnsearch(srcm, tgtm, 'k', 1);
            matched = D(:,1) < opt.thresh;
            nn(~matched) = nan;
            
        end
        
        function [ds, rank] = Table2DavidScore(tt)
            % Uses data from profile table to compute the total david score
            % rather than the CheeseSquare file
            [~, ~, gn] = unique(tt.GroupNumber);
            ds = zeros(size(tt, 1), 1);
            rank = zeros(size(tt, 1), 1);
            for g=1:max(gn)
                currgroup = tt(gn == g, :);
                [~, ~, mousenum] = unique(currgroup.MouseNumber);
                dstable = Q.accumrows(mousenum, currgroup.PairwiseChaseCount, @sum);
                currds = Hierarchy.DavidScore(dstable);
                currrank = Q.rank(-currds);
                ds(gn == g) = currds(mousenum);
                rank(gn == g) = currrank(mousenum);
            end
        end
       
        %%
        function t = CombineProfiles(p1, p2)
            % Merge two profile tables into one
            %%
            if istable(p1)
                t1 = p1;
            elseif isstruct(p1) || isobject(p1)
                t1 = p1.Tables.Individual;
            end
            if istable(p2)
                t2 = p2;
            elseif isstruct(p1) || isobject(p1)
                t2 = p2.Tables.Individual;
            end
            
            if isempty(p1)
                t = t2;
                return;
            elseif isempty(p2)
                t = t1;
                return;
            end

            %%
            
            if size(t1, 2) ~= size(t2, 2)
                warning('Profile tables have different number of columns');
                common1 = ismember(t1.Properties.VariableNames, t2.Properties.VariableNames);
                common2 = ismember(t2.Properties.VariableNames, t1.Properties.VariableNames);
                t1 = t1(:, common1);
                t2 = t2(:, common2);
            else
                if ~strcmp([t1.Properties.VariableNames{:}], [t2.Properties.VariableNames{:}])
                    ME = MException('CheeseSquare:ProfileTablesMismatch', 'Profile tables mismatch');
                    throw(ME);
                end
            end

            if ~isempty(t1)
                t1 = CheeseSquare.CombineProfiles(table(), t1);
            end
            %% update numbers
            tgtfields = {'GroupNumber', 'MouseNumber'};
            for j=1:length(tgtfields)
                field = tgtfields{j};
                
                [~, ~, u1] = unique(t1.(field));
                [~, ~, u2] = unique(t2.(field));
                
                t1.(field) = u1;
                t2.(field) = max(u1) + u2;
            end
            
            %% ensure no duplicate row names
            retry = true;
            count = 1;
            while retry
                retry = false;
                for i=1:length(t2.Properties.RowNames)
                    if ismember(t2.Properties.RowNames{i}, t1.Properties.RowNames)
                        %%
                        retry = true;
                        try
                            t2.Properties.RowNames{i} = [t2.Properties.RowNames{i} '_' num2str(count)];
                        catch
                        end
                    end
                end
                count = count + 1;
            end
            
            %%
            t = [t1; t2];
        end
        
        %%
        function objs = RunOnGroup(obj, cmd, issave)
            % RunOnGroup(obj, cmd, issave) Runs a specific command on all
            % days of a specific group
            if nargin < 3
                issave = true;
            end
            objs = CheeseSquare.LoadGroup(obj);
            for i=1:length(objs)
                Console.Message('running ''%s'' on ''%s''', func2str(cmd), objs{i}.Prefix);
                objs{i} = cmd(objs{i});
                if issave
                    objs{i}.Save;
                end
            end
        end
        
        function objs = LoadGroup(obj, days)
            % Load all days of a specific group
            % Syntax: objs = LoadGroup(obj, days)
            if ischar(obj)
                fname = MyFilename(obj);
                meta = CheeseSquare.ParseName(fname.Full);
            elseif ischar(obj.Video)
                fname = MyFilename(obj.Video);
                meta = obj.Meta;
            else
                fname = MyFilename(obj.Video.Filename);
                meta = obj.Meta;
            end
            template = CheeseSquare.GeneratePrefix(meta.GroupType, meta.GroupId, [], meta.CameraId);
            if nargin < 2
                days = 1:99;
            end
            objs = cell(1, length(days));
            parfor i=1:length(days)
                currfile = fname.SetName(sprintf(template, days(i))).Full;
                if exist(currfile, 'file') || exist([currfile '.obj.mat'], 'file')
                    Console.Message('loading ''%s''', sprintf(template, days(i)));
                    objs{i} = CheeseSquare(currfile);
                end
            end
            objs = objs(cellfun(@(x) ~isempty(x), objs));
        end
    end
    
    methods (Hidden)
        function prefix = GetPrefix(obj)
            % Get object's prefix (such as 'SC.exp0001.day02.cam02')
            filename = obj.Filename.Full;
            [begi, endi] = regexp(filename, '[^\\\/]*\.exp\d+\.day\d+\.cam\d+');
            prefix = filename(begi:endi);
        end
        
        
        function y = GetField(obj, field, default)
            % Get a specific field from object file
            if nargin < 3
                default = [];
            end
            if obj.EnsureField(field)
                try
                    y = eval(['obj.' field]);
                catch
                    y = default;
                end
            else
                y = default;
            end
        end
        
        function y = EnsureField(obj, field) %#ok<INUSL>
            % ensure a field exists
            try
                y = eval(['~isempty(obj.' field ')']);
            catch
                y = false;
            end
        end
        
    end
end