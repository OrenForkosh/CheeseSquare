function [objs, idx] = CheeseGroupLoad(obj)
%% Load all days of an experiment. Should be replaced by CheeseFarm
%
%       Created by OREN FORKOSH
%
if isa(obj, 'CheeseSquare')
    error 'unsupported yet...';
else
    p1 = CheeseSquare.ParseName(obj.VideoFile);
    idx = 0;
    [filepath, filename, ext] = fileparts(obj.VideoFile);
    files = dir(fullfile(filepath, [regexprep(filename, 'day[0-9]*', 'day*'), ext]));
    objs = struct();
    for i=1:length(files)
        p2 = CheeseSquare.ParseName(files(i).name);
        if strcmp(p1.GroupType, p2.GroupType) && p1.DayId == p2.DayId && p1.GroupId == p2.GroupId && p1.CameraId == p2.CameraId
            idx = i;
            curr = obj;
        else
            curr = CheeseSquare(fullfile(filepath, files(i).name)).ToClassical;
        end
        id = regexprep(curr.FilePrefix, '\.', '_');
        objs.(id) = curr;
    end
end