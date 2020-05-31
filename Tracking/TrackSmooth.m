function obj = TrackSmooth(obj, span)
% Smooth mouse trajectories with a sliding windows of length 'span'
if nargin == 1
    span = 5;
end
win = ones(1, span) / span;
for s=1:obj.nSubjects
    x = conv(obj.x(s, obj.valid), win);
    y = conv(obj.y(s, obj.valid), win);
    obj.x(s, obj.valid) = x(1:end-span+1);
    obj.y(s, obj.valid) = y(1:end-span+1);
end
