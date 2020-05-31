function obj = myMMReader(filename, frame, bkg)
% Read video file

persistent videoObj;

if isempty(videoObj) || ~strcmp(videoObj.Filename, filename)
    videoObj = MyVideoReader(filename);
end

if exist('frame', 'var')
    videoObj.FrameNumber = frame;
    obj = videoObj.CurrentFrame;
else
    obj = videoObj;
end

