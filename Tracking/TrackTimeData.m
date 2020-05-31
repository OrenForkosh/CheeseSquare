function obj = TrackTimeData(obj)
% Add recording times meta-data to CheeseObject
obj.RecordingEnd = datevec(GetFileCreationDate(obj.VideoFile));
obj.RecordingDuration = datevec(num2str(obj.nFrames / obj.FrameRate), 'SS');
obj.RecordingDuration(1:3) = 0;
obj.RecordingStart = datevec(datenum(obj.RecordingEnd) - datenum(obj.RecordingDuration));
