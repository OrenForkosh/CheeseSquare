function obj = TrackLoad(obj, subs)
% TrackLoad(obj, subs) Loads older versions of CheeseSquare objects
%   obj is either a struct or a filename
% TrackLoad({"SC", 1, 1}, subs)

o = TrackLoadBkgObj(obj);
if ~isempty(o)
    obj = o;
    o = [];
end

if iscell(obj)
    SocialExperimentData;
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

SocialExperimentData;
try
    %%
    prop = regionprops(obj.ROI.Regions{strcmp(obj.ROI.RegionNames, GroupsData.Proportions.Region)}, 'MajorAxisLength', 'MinorAxisLength');
    obj.PixelsPerCM = mean([prop.MinorAxisLength / min(GroupsData.Proportions.Size), prop.MajorAxisLength / max(GroupsData.Proportions.Size)]);
catch
end
