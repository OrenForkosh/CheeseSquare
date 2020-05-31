function obj = TrackCatObject(obj, other)
% Concatenate two CheeseSquare objects
if isempty(obj) || isempty(fieldnames(obj))
    obj = other;
    return;
end

%%
data = TrackParse(obj);
if ~isfield(obj, 'Days')
    obj.Days = data.Day;
end
if ~isfield(other, 'Days')
    odata = TrackParse(other);
    other.Days = odata.Day;
end
%% append
fields = {'zones', 'valid', 'x', 'y', 'hidden', 'sheltered', 'speed', 'Days'};
for i=1:length(fields)
    if isfield(other, fields{i})
        if ~isfield(obj, fields{i})
            obj.(fields{i}) = [];
        end
        obj.(fields{i}) = [obj.(fields{i}), other.(fields{i})];
    end
end

%% append texts
fields = {'VideoFile', 'Source'};
for i=1:length(fields)
    if isfield(other, fields{i})
        if ~isfield(obj, fields{i})
            obj.(fields{i}) = {};
        elseif ischar(obj.(fields{i}))
            obj.(fields{i}) = {obj.(fields{i})};
        end
        obj.(fields{i}) = {obj.(fields{i}){:}, other.(fields{i})};
    end
end

%%
data.Day = obj.Days;
obj.FilePrefix = sprintf('%s.exp%04d.day%s.cam%02d', data.Type, data.GroupID, sprintf('%02d', data.Day), data.CameraID);
%% add
if isfield(other, 'nFrames')
    if ~isfield(obj, 'nFrames')
        obj.nFrames = 0;
    end
    obj.nFrames = obj.nFrames + other.nFrames;
end

%%
offset = 0;
if isfield(obj, 'RecordingEnd') && isfield(other, 'RecordingStart') && ~isempty(other.RecordingEnd) && ~isempty(other.RecordingStart)
    offset = etime(other.RecordingStart, obj.RecordingEnd);
end
if ~isfield(obj, 'time') || isempty(obj.time)
    obj.time = [];
    last = 0;
else
    last = obj.time(end) + offset;
end

if isfield(other, 'time')
    obj.time = [obj.time, other.time + last];
end

%%
preserve = {'ROI', 'Colors'};
f = fieldnames(obj);
for i=1:length(f)
    if isstruct(obj.(f{i}))
        if ~any(strcmp(f{i}, preserve))
            obj = rmfield(obj, f{i});
        end
    end
end

