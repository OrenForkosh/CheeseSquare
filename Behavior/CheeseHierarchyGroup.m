function source = CheeseHierarchyGroup(obj)
% Compute dominance hierarchy for all days of an experiment
%   source = CheeseHierarchyGroup(obj) where obj is a CheeseSquare object.
%   Returns all the CheeseSquare objects of a group (all days) with the
%   hiearchical data in the .hierarchy.group field

% get the video filename
if ischar(obj.Video)
    fname = MyFilename(obj.Video);
else
    fname = MyFilename(obj.Video.Filename);
end
% parse video filename
template = CheeseSquare.GeneratePrefix(obj.Meta.GroupType, obj.Meta.GroupId, [], obj.Meta.CameraId);
objs = {};
idx = 1;
sourceidx = 1;
% look for all available days of an experiment (from day 1 till day 99)
for i=1:99
    currfile = fname.SetName(sprintf(template, i)).Full;
    % if file exist, load it...
    if exist(currfile, 'file') || exist([currfile '.obj.mat'], 'file')
        Console.Message(1, 'loading ''%s''', sprintf(template, i));
        objs{idx} = CheeseSquare(currfile);
        succ = false;
        try
            % do behavioral analysis on each loaded file
            %if ~Q.isfield(objs{idx}, 'Hierarchy') || ~Q.isfield(objs{idx}.Hierarchy, 'ChaseEscape') || ~Q.isfield(objs{idx}.Hierarchy, 'AggressiveChase')  || ~Q.isfield(objs{idx}.Hierarchy, 'Contacts')
                track = objs{idx}.ToClassical;
                track = CheeseHierarchy(track);
                objs{idx}.Hierarchy = track.Hierarchy;
                objs{idx}.Hierarchy.Contacts = track.Contacts;
            %end
            succ = true;
        catch me
            succ = false;
            Console.Warning(me, 'unable to classify hierarchy in ''%s''', objs{idx}.Prefix);
            objs(idx) = [];
        end
        if succ
            if strcmp(sprintf(template, i), obj.Prefix)
                sourceidx = idx;
            end
            idx = idx + 1;
        end
    end
end
%% Create the .Hierarchy.Group field in the CheeseSquare object
%  The field contains two subfields: 'ChaseEscape' and 'AggressiveChase'
fields = {'ChaseEscape', 'AggressiveChase'};
for fidx=1:length(fields)
    type = fields{fidx};
    % initialize
    for i=1:length(objs)
        objs{i}.Hierarchy.Group.(type) = struct();
    end
    % assign hierarchy data from all days
    for j=1:length(objs)
        for i=1:length(objs)
            objs{j}.Hierarchy.Group.(type).ChaseEscape(:, :, i) = objs{i}.Hierarchy.(type).ChaseEscape;
            objs{j}.Hierarchy.Group.(type).Duration(i, 1) = sum(objs{i}.Tracking.valid) / objs{i}.Video.FrameRate;
            objs{j}.Hierarchy.Group.(type).DayID(i, 1) = objs{i}.Meta.DayId;
        end
        objs{j}.Hierarchy.Group.(type).DavidScore = DavidScore(sum(objs{1}.Hierarchy.Group.(type).ChaseEscape, 3))';
        %%
        [~, o] = sort(objs{1}.Hierarchy.Group.(type).DavidScore);
        seq = 1:objs{1}.nSubjects; seq(o) = seq;
        objs{j}.Hierarchy.Group.(type).Order = seq(:)';
    end
end
%% save all objects to files
for i=1:length(objs)
    Console.Message(1, 'saving ''%s''', objs{i}.Prefix);
    objs{i}.Save;
end
%% return all the loaded objects
source = obj;
try
    source = objs{sourceidx};
catch
end

