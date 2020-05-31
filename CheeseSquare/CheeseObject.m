classdef CheeseObject
    % CHEESEOBJECT Similar to CheeseSquare file used for back compatibility
    % with older files and with older code
    
    methods (Static = true)
        
        function obj = Save(obj)
            % Save object to file
            atomic = true;
            obj = TrackOptimize(obj);
            
            filename = [obj.OutputPath obj.FilePrefix '.obj.mat'];
            fprintf('# (>) saving tracking object to ''%s''\n', filename);
            if atomic
                save([filename '.temp'], '-struct', 'obj', '-v7');
                movefile([filename '.temp'], filename);
            else
                save(filename, '-struct', 'obj', '-v7');
            end
        end
        
        function obj = Create(videoFile)
            % Create(videoFile) creates a struct that contains the video
            % filename, metadata about the video, as well as default
            % parameters to various stages of the tracking algorithm
            % 
            if nargin == 0
                videoFile = '';
            end
            obj.VideoFile = videoFile;
            obj.Name = '';
            obj.nSubjects = 4;          % the number of subjects in the arena
            
            %obj.Unix.VideoPath = '../Movies/';
            %obj.Windows.VideoPath = 'x:\tests\SocialMice\Movies\';
            %%
            obj.Colors = struct();
            %%
            obj.StartTime = -1; % video start time [sec]
            obj.EndTime = -1;   % video end time [sec]
            %%
            obj.OutputPath = 'Res/';
            obj.Output = false;
            obj.OutputToFile = true;
            obj.OutputInOldFormat = true;
            obj.FilePrefix = regexprep(regexprep(obj.VideoFile, '\.[^\.]*$', ''), '.*[\\/]', '');
            %%
            obj.ValidROI = [];
            %% meta
            obj.NoiseThresh = 10;        %
            obj.UseAdaptiveThresh = true;        %
            %   - background:
            obj.nBkgFrames = 40;        % no. of frames to use when computing background image
            obj.ErodeRadios = 2;
            %   - colors:
            obj.ColorsFile = '';
            obj.MinNumOfObjects = 4;
            obj.nColorFrames = inf;      % no. of frames to use to determine the subject colors
            obj.nColorBins = 20;         % no. of bins in color histogram
            obj.nKMeansRepeats = 20;     % no. of time to run k-means algorithm for optimal clustring of colors
            %% segementation:
            obj.MinNumOfPixels = 40;
            obj.MinNumOfPixelsAfterErode = 10;
            obj.MinSolidity = .8;          % minimal solidity of segments
            %
            obj.VideoScale = .25;            % working resolution relative to original movie
            obj.MaxNumOfCents = 10;
            obj.MaxNumOfObjects = 20;
            obj.MinSegmentGap = 3;             % important if there are small hidden areas (like wallthrough walls)
            %% path tracking
            obj.ColorMatchThresh = 0.1; % minimum certainty probability required to assign label to segment
            obj.MaxHiddenDuration = 25; % [frames] maximal time the tracking can lose a subject
            obj.NumOfProbBins = 20; %
            obj.MinJumpProb = 0.01; %
            obj.AllowHiddenTransitions = false; %
            obj.MaxIsolatedDistance = 12; % remove tracks that were isolated from other by more than this
            %%
            obj = CheeseObject.VideoFileProperties(obj);
            if nargin > 0
                obj = CheeseObject.Load(obj);
            end
            try
                obj.VideoFileSize = GetFileSize(obj.VideoFile);
            catch
                obj.VideoFileSize = [];
            end
            
        end
        
        function obj = Load(obj, subs)
            % TrackLoad(obj, subs)
            %   obj is either a struct or a filename
            % TrackLoad({"SC", 1, 1}, subs)
            
            o = CheeseObject.LoadInBackground(obj);
            if ~isempty(o)
                obj = o;
                o = [];
            end
            
            if iscell(obj)
                %SocialExperimentData;
                params = {obj{:}};
                if ischar(params{1})
                    list = {GroupsData.Experiments{GroupsData.(params{1}).idx}};
                    params = {params{2:end}};
                else
                    list = GroupsData.Experiments;
                end
                curr = list{params{1}};
                if length(params) > 1
                    days = params{2};
                else
                    days = 1:GroupsData.nDays;
                end
                obj = struct();
                for day = days
                    if day > 0
                        prefix = sprintf(curr, day);
                    elseif day == 0
                        prefix = regexprep(curr, '%02d', sprintf('%02d', 1:GroupsData.nDays));
                    elseif day == -1
                        prefix = regexprep(curr, '%02d', sprintf('%02d', 2:GroupsData.nDays));
                    end
                    file = [GroupsData.OutputPath prefix '.obj.mat'];
                    if exist('subs', 'var')
                        obj = TrackCatObject(obj, TrackLoad(file, subs));
                    else
                        obj = TrackCatObject(obj, TrackLoad(file));
                    end
                end
                return;
            end
            
            commonSubs = {'time', 'x', 'y', 'zones', 'ROI', 'valid', 'hidden', 'sheltered', 'nFrames', 'Colors', 'RecordingStart', 'RecordingEnd', 'FrameRate', 'dt'};
            
            if ischar(obj)
                filename = obj;
                if exist('subs', 'var')
                    if ischar(subs)
                        str = subs;
                        subs = cell(1);
                        subs{1} = str;
                    end
                    fprintf(['# (<) loading tracking object from ''' filename ''': ']);
                    for i=1:length(subs)
                        if i==1
                            fprintf('(');
                        else
                            fprintf(', ');
                        end
                        fprintf(subs{i});
                        if i==length(subs)
                            fprintf(')');
                        end
                    end
                    fprintf('\n');
                    if any(strcmpi(subs, 'common'))
                        subs = {subs{~strcmpi(subs, 'common')}, commonSubs{:}};
                    end
                    defaultSubs = fields(CheeseObject.Create());
                    allSubs = {subs{:} defaultSubs{:}};
                    warning off MATLAB:load:variableNotFound;
                    obj = load(filename, allSubs{:});
                    warning on MATLAB:load:variableNotFound;
                    if isempty(fields(obj)) % for back compatibility
                        load(filename);
                        newobj = struct();
                        for i=1:length(allSubs)
                            if isfield(obj, allSubs{i})
                                newobj.(allSubs{i}) = obj.(allSubs{i});
                            end
                        end
                        obj = newobj;
                    end
                else
                    %fprintf(['# (<) loading tracking object from ''' filename '''\n']);
                    obj = load(filename);
                    % for back compatibility:
                    if isfield(obj, 'obj')
                        obj = obj.obj;
                    end
                end
                obj.Source = filename;
            end
            
            if isfield(obj, 'valid');
                obj.nValid = length(obj.valid);
            end
            
            if isfield(obj, 'ROI') && isfield(obj.ROI, 'ZoneNames')
                obj.ROI.nZones = length(obj.ROI.ZoneNames);
            end
            
            %if isunix
            if isfield(obj, 'Unix') && ~isempty(obj.Unix.VideoPath)
                obj.VideoFile = [obj.Unix.VideoPath regexprep(obj.VideoFile, '.*[\\/]', '')];
            end
            % else
            %     if isfield(obj, 'Windows') && ~isempty(obj.Windows.VideoPath)
            %         obj.VideoFile = [obj.Windows.VideoPath regexprep(obj.VideoFile, '.*[\\/]', '')];
            %     end
            % end
            
            if ~isfield(obj, 'Source')
                obj.Source = [obj.OutputPath obj.FilePrefix '.obj.mat'];
            end
            
            try
                obj = TrackTimeData(obj);
            catch me
                %    warning(['Unable to determine recording times (' me.message ')']);
            end
            
            convert = {'regions', 'zones', 'sheltered', 'speed', 'hidden'};
            for i=1:length(convert)
                if isfield(obj, convert{i})
                    obj.(convert{i}) = double(obj.(convert{i}));
                end
            end
            
            %SocialExperimentData;
            try
                %%
                prop = regionprops(obj.ROI.Regions{strcmp(obj.ROI.RegionNames, GroupsData.Proportions.Region)}, 'MajorAxisLength', 'MinorAxisLength');
                obj.PixelsPerCM = mean([prop.MinorAxisLength / min(GroupsData.Proportions.Size), prop.MajorAxisLength / max(GroupsData.Proportions.Size)]);
            catch
            end
        end
        
        function obj = LoadInBackground(bkgobj)
            obj = [];
            if isfield(bkgobj, 'Type') && strcmp(bkgobj.Type, 'bkgobj')
                %%
                fprintf('# (<) - loading from temp file...\n');
                count = 1;
                first = true;
                running = true;
                while running
                    strud = ReadStrudels([bkgobj.SourceObj.OutputPath bkgobj.SourceObj.FilePrefix '.output']);
                    if count > 1
                        fprintf('\b\b\b\b\b\b');
                    end
                    if isfield(strud, 'status')
                        switch strud.status
                            case 'done'
                                if count > 1
                                    fprintf('\n');
                                end
                                obj = TrackLoad([bkgobj.SourceObj.OutputPath bkgobj.SourceObj.FilePrefix '.obj.mat']);
                                obj.FilePrefix = bkgobj.FilePrefix;
                                obj.OutputPath = bkgobj.OutputPath;
                                running = false;
                                fprintf('# (<) done!\n');
                            case 'failed'
                                if count > 1
                                    fprintf('\n');
                                end
                                throw(MException('Tracking:InvalidBackgroundObject', ['background object generation failed for ' bkgobj.FilePrefix ' (' bkgobj.SourceObj.FilePrefix '), due to ' strud.message]));
                                running = false;
                            case 'running'
                                running = true;
                            otherwise
                                throw(MException('Tracking:InvalidBackgroundObject', ['background object generation failed for ' bkgobj.FilePrefix ' (' bkgobj.SourceObj.FilePrefix '), due to ' strud.message]));
                        end
                    end
                    if running
                        if first
                            fprintf(['# waiting for background object: ' bkgobj.FilePrefix ' (' bkgobj.SourceObj.FilePrefix ')   ']);
                            logid = fopen(bkgobj.Process.logfile);
                        end
                        tline = fgetl(logid);
                        while ischar(tline)
                            fprintf('\n# (<)   . %s', tline);
                            tline = fgetl(logid);
                        end
                        switch mod(count-1,4)
                            case 0
                                fprintf('[   ] ');
                            case 1
                                fprintf('[#  ] ');
                            case 2
                                fprintf('[## ] ');
                            case 3
                                fprintf('[###] ');
                        end
                        count = count + 1;
                        pause(.5);
                    end
                    first = false;
                end
            end
        end
        
        function obj = VideoFileProperties(obj)
            %%
            obj.VideoHeight = 0;
            obj.VideoWidth = 0;
            obj.nFrames = 0;
            obj.dt = 0;
            obj.FrameRate = 0;
            %
            if ~isempty(obj.VideoFile) && exist(obj.VideoFile, 'file')
                try
                    xyloObj = MyVideoReader(obj.VideoFile);
                    obj.VideoHeight = xyloObj.Height;
                    obj.VideoWidth = xyloObj.Width;
                    obj.nFrames = xyloObj.NumberOfFrames;
                    
                    try % check for a file eventually used by 'MovieSequence' in order to know the number of frames in the user-selected range
                        slashes = strfind(obj.VideoFile, '\');
                        actDir = obj.VideoFile(1: slashes(end)-1);
                        if exist([actDir '\oldFullMatFiles'], 'dir') && exist([actDir '\oldFullMatFiles\' ], 'file')
                            fid = fopen([actDir '\oldFullMatFiles\' obj.VideoFile(slashes(end)+1:end-length('.avi')) '.txt'], 'r');
                            obj.nFrames = fscanf(fid, '%d');
                            fclose(fid);
                        end
                    catch
                    end
                    if obj.nFrames == 0
                        temp = mmread(obj.VideoFile, 1);
                        obj.nFrames = abs(temp.nrFramesTotal);
                    end
                    obj.dt = 1/xyloObj.FrameRate;
                    obj.FrameRate = xyloObj.FrameRate;
                catch me
                    fprintf('# WARNING! unable to open video file "%s": %s', obj.VideoFile, me.message);
                end
            end
        end
        
        function arg = AuxParseStructArguments(defaults, arg)
            
            f = fieldnames(defaults);
            for i=1:length(f)
                if ~isfield(arg, f{i})
                    arg.(f{i}) = defaults.(f{i});
                end
            end
        end
        
    end
end