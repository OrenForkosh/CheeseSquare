classdef MyVideoReader < Serializable
    % MyVideoReader handle video files
    properties
        FrameRate = 0;              % Frame rate
        Height = 0;                 % Frame height
        Width = 0;                  % Frame width
        RealNumberOfFrames = 0;     % Actual number of frames in file
        NumberOfFrames = 0;         % Estimated number of frames in file (from duration x framerate)
        Duration = 0;               % Duration of the video (in secs)
        StartDateTime = [];         % Recording data
        Meta = struct();            % Meta data
    end
    
    properties (Dependent)
        FrameNumber = 0;            % Current frame number
        Time = 0;                   % Current frame time
        Filename;                   % Video filename
    end
    properties (Access=private)
        % Video buffer which stores several frames for quicker access
        Buffer_ = [];        
        %
        BufferFrameRange_ = [];     % Video buffer frame range
        BufferSize_ = 5;            % Size of video buffer
        Filename_ = '';             % Video filename
        Time_ = 0;                  % Current frame number
        Source_ = [];               
        Current = [];
    end
    
    properties (Constant, Hidden)
        UseOld = false;
    end
    
    
    methods
        function obj = MyVideoReader(varargin)
            % Init the MyVideoReader object
            %   MyVideoReader(filename) Create object associated with the
            %   video file 'filename'
            %   
            %   MyVideoReader(another) create object from another MyVideoReader
            %   object
            %
            %   MyVideoReader(struct) create object from a structure
            %   containing in the 'Filename' of 'VideoFile' fields the name
            %   of the video as well as additional data (like frame-rate,
            %   height, etc.)
            %   
            if nargin < 1
                error 'Not enough input arguments.'
            end
            switch class(varargin{1})
                case 'char'
                    obj.Filename = varargin{1};
                case 'struct'
                    if isfield(varargin{1}, 'VideoFile')
                        obj = Q.cpfield(varargin{1}, obj, {'FrameRate', 'Height', 'Width', 'NumberOfFrames', 'RealNumberOfFrames', 'Duration', 'StartDateTime', 'Meta'}, true);
                        obj.Filename = varargin{1}.VideoFile;
                    else
                        obj = Q.cpfield(varargin{1}, obj, {'FrameRate', 'Height', 'Width', 'NumberOfFrames', 'RealNumberOfFrames', 'Duration', 'StartDateTime', 'Meta'}, true);
                        obj.Filename = varargin{1}.Filename;
                    end
                case 'MyVideoReader'
                    obj = varargin{1};
                otherwise
                    error(['MyVideoReader cannot except input argument of type ' class(varargin{1})]);
            end
            if nargin >= 2
                obj.FrameNumber = varargin{2};
            end
        end
        
        function obj = set.FrameNumber(obj, value)
            % Set the current frame number and load the frame
            obj.Time = value / obj.FrameRate;
            obj = UpdateBuffer(obj);
        end

        function value = get.FrameNumber(obj)
            % Get the current frame number
            value = max(round(obj.Time_ * obj.FrameRate + 1), 1);
        end

        function obj = set.Time(obj, value)
            % Set the current frame time and load the frame
            obj.Time_ = value;
            try
                obj.Source_.CurrentTime = obj.Time_;
                obj.Current = obj.Source_.readFrame;
            catch
                obj.Current = [];
            end
        end
        
        function value = get.Time(obj)
            % Get the current frame time
            value = obj.Time_;
        end
        
        function value = get.Filename(obj)
            % Get the current filename
            value = obj.Filename_;
        end
        
        function obj = set.Filename(obj, value)
            % Set the current filename and load the new file
            obj.Filename_ = value;
            try
                if ~MyVideoReader.UseOld
                    meta = VideoReader(obj.Filename_);
                    obj.Height = meta.height;
                    obj.Width = meta.width;
                    obj.FrameRate = meta.FrameRate;
                    obj.Duration = meta.Duration;
                    obj.NumberOfFrames = floor(meta.Duration * meta.FrameRate);
                    obj.RealNumberOfFrames = floor(obj.NumberOfFrames);
                    obj.Source_ = meta;
                else
                    meta = mmread(obj.Filename_, 1);
                    obj.Height = meta(end).height;
                    obj.Width = meta(end).width;
                    obj.FrameRate = ceil(meta(end).rate);
                    obj.Duration = meta(end).totalDuration;
                    obj.NumberOfFrames = floor(obj.FrameRate * meta(end).totalDuration);
                    obj.RealNumberOfFrames = floor(abs(meta(end).nrFramesTotal));
                end
            catch
            end
            
            % trying to read meta-data from assorted xml file
            try
                videofile = [obj.Filename_, '.xml'];
                meta.CreationDate = addtodate(MyFilename(obj.Video.Filename).CreationDate, -round(obj.Video.NumberOfFrames / obj.Video.FrameRate*1000), 'millisecond');
                meta.Remarks = '';
                if exist(videofile.Full, 'file')
                    metaxml = xmlread(videofile.Full);
                    cdate = char(metaxml.getElementsByTagName('Start').item(0).getTextContent);
                    obj.StartDateTime = datenum(cdate, 'yyyy-mm-ddTHH:MM:SS.FFF');
                    meta.Remarks =  char(metaxml.getElementsByTagName('Remarks').item(0).getTextContent);
                end
                Meta = meta;
            catch
            end
        end

        function frame = CurrentFrame(obj)
            % Get current video frame
            if ~MyVideoReader.UseOld
                frame = obj.Current;
            else
                currFrame = obj.FrameNumber;
                obj = UpdateBuffer(obj);
                num = currFrame - obj.BufferFrameRange_(1) + 1;
                if length(obj.Buffer_) >= num
                    frame = obj.Buffer_(num).cdata;
                else
                    frame = zeros(obj.Height, obj.Width, 3);
                end
            end
        end
        
        function [obj, frame] = NextFrame(obj)
            % Load the next video frame from file
            if ~MyVideoReader.UseOld
                frame = obj.Source_.readFrame;
                obj.Time_ = obj.Source_.CurrentTime;
            else
                obj.FrameNumber = obj.FrameNumber + 1;
                frame = obj.CurrentFrame;
            end
        end
        
        function Show(obj)
            % Show current frame in a figure
            imagesc(obj.CurrentFrame);
            axis off;
        end
            
        function Play(obj)
            % Play the video file in a figure
            h = tic;
            rt = 0;
            rtfactor = 0.95;
            for i=1:obj.NumberOfFrames
                obj.FrameNumber = i;
                imagesc(obj.CurrentFrame);
                time = toc(h);
                if rt == 0
                    rt = obj.Time/time;
                else
                    rt = rtfactor * rt + (1 - rtfactor) * obj.Time/time;
                end
                title(sprintf('%d (x%.1f RT)', i, rt));
                axis off;
                drawnow;
            end
        end
    end
    
    methods (Access = protected)
        function obj = UpdateBuffer(obj)
            % Update video buffer if needed
            if ~MyVideoReader.UseOld
                return;
            end
            currFrame = obj.FrameNumber;
            if ~isempty(obj.Buffer_) && ...
                    obj.BufferFrameRange_(1) <= currFrame && ...
                    obj.BufferFrameRange_(2) >= currFrame
            else
                meta = mmread(obj.Filename, [], [obj.Time_ obj.Time_+(obj.BufferSize_ + .5)/obj.FrameRate]);
                if isempty(meta(end).frames)
                    pause(.2);
                    meta = mmread(obj.Filename, [], [obj.Time_ obj.Time_+(obj.BufferSize_ + .5)/obj.FrameRate]);
                end
                if ~isempty(meta(end).frames)
                    obj.Buffer_ = meta(end).frames;
                    obj.BufferFrameRange_ = [currFrame, currFrame + length(obj.Buffer_) - 1];
                else
                    obj.Buffer_ = [];
                    obj.BufferFrameRange_ = [-1 -1];
                end
            end
        end
    end
end