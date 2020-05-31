function obj = TrackReadVideoFileProperties(obj)
%%
obj.VideoHeight = 0;
obj.VideoWidth = 0;
obj.nFrames = 0;
obj.dt = 0;
obj.FrameRate = 0;
%
if ~isempty(obj.VideoFile) && exist(obj.VideoFile, 'file')
    try
        xyloObj = myMMReader(obj.VideoFile);
        obj.VideoHeight = xyloObj.Height;
        obj.VideoWidth = xyloObj.Width;
        obj.nFrames = xyloObj.NumberOfFrames;
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
