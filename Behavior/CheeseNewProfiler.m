function [obj, Profile, succ] = CheeseNewProfiler(obj, varargin)
% CHEESENEWPROFILER Compute behavioral profile for each mouse (and group)
%
%   [obj, Profile] = CheeseNewProfiler(obj) Computed the
%   behavioral profile of individual mice as well as group properties from 
%   CheeseSquare object 'obj'. The function returns the object, Profile
%   struct (details below).
%
%   [obj, Profile] = CheeseNewProfiler(obj, varargin) Computed the
%   behavioral profile with the optional parameter in varargin. The
%   possible parameters are:
%       Locations - Names of valid ROIs in the arena
%       HighPlaces - Elevated areas (such as a block)
%       FoodWater - Feeding sources
%
%       TrainMaxEnt - Should train MaxEnt model
%
%       DistanceNoiseThreshold - Minimal movement distance (in cm).
%           The algorithm ignores movements that were below this
%           thereshold. Needed for robustness.
%       DistanceNoiseDuration - The temporal resolution (in sec) by which the
%           position is sampled when determining behaviors such as speed or
%           distance. Needed for robustness.
%       SpeedWinThreshold - Like 'DistanceNoiseThreshold' for angular and
%           tangential velocities.
%       SpeedWinDuration - Like 'DistanceNoiseDuration' for angular and
%           tangential velocities.
%       MinTimeInZone - Minimum duration a mouse can spend in a region
%       (otherwise the region is ignored)
%
%   The profile fields:
%       Most useful:
%           Tables.Individual - Behaviors per animal organized in a table
%       In addition:
%           Tables.Group - Group properties organized in a table
%           Individual - List of behaviors per individual
%           Group - List of group behaviors
%           Pairwise - Square matrices of pairwise behaviors (no. of chases,
%               etc.)
%           Segments - List of behaviors per day segment (typically 6 segments
%               per day)
%           Events - List of all events and their time (such as, chase-escape)
%       Meta data:
%           HighPlaces - List of elevated ROIs
%           FoodWater - List of feeding ROIs
%%
p = inputParser;
% names of valid ROIs in the arena:
p.addOptional('Locations', {'Open', 'Feeder1', 'Feeder2', 'Water', 'Water2', 'SmallNest', 'Labyrinth', '[Ramp1]', '[Ramp2]'});
% elevated areas (such as a block):
p.addOptional('HighPlaces', {'BigNest', 'Block', '[Ramp1]', '[Ramp2]'}); 
% feeding sources:
p.addOptional('FoodWater', {'Feeder1', 'Feeder2', 'Water', 'Water2'});

% should train MaxEnt model:
p.addOptional('TrainMaxEnt', true);

% Minimal movement distance (in cm). The algorithm ignores movements that 
% were below thiS thereshold. Needed for robustness:
p.addOptional('DistanceNoiseThreshold', 2); % [cm]
% The temporal resolution (in sec) by which the position is sampled when 
% determining behaviors such as speed or distance. Needed for robustness.
p.addOptional('DistanceNoiseDuration', 1); % [sec]
% Like 'DistanceNoiseDuration' for angular and tangential velocities:
p.addOptional('SpeedWinDuration', 1/4); % [sec]
% Like 'DistanceNoiseThreshold' for angular and tangential velocities:
p.addOptional('SpeedWinThreshold', 5); % [cm]

% Minimum duration a mouse can spend in a region:
p.addOptional('MinTimeInZone', 1/4); % [sec]
p.parse(varargin{:});
opt = p.Results;

succ = true;

%% Read scaling from pixels to cm
% Either directly from preprocessing data or by infering it from the size
% of the arena
try
    if Q.isfield(obj, 'Meta.Scale.PixPerCM')
        PixelsPerCM = obj.Meta.Scale.PixPerCM;
    else
        PixelsPerCM = mean([obj.Meta.Scale.ArenaCoord(3) / obj.Meta.Scale.ArenaWidth, obj.Meta.Scale.ArenaCoord(4) / obj.Meta.Scale.ArenaHeight]);
    end
catch
end

%% Load CheeseSquare object
obj = CheeseSquare(obj);
%obj = TrackLoad(obj);

%% Initialize Profile structure
Profile = struct;
Profile.Individual = struct;
Profile.Group = struct;
Profile.Pairwise = struct;
Profile.Segments = struct;
Profile.Events = struct;

%% Determine video frame rate 
% Either directly from video or from the timestamps. If cannot determine
% fps then assume the frame rate is 25fps.
fps = obj.Video.FrameRate;
if fps == 0
    try
        fps = median(1./diff(obj.Tracking.time));
    catch
    end
    if fps == 0
        fps = 25;
    end
    obj.Video.FrameRate = fps;
end

%% Initialize data properties
try
    TotalNumberOfFrames = sum(obj.Tracking.valid); % number of frames in video
    VideoDuration = sum(obj.Tracking.valid) / obj.Video.FrameRate / 60 /60; % duration of video in hours
    TimeOutside = ... % amount of time each mouse spent outside the nest (in hours)
        sum(~obj.Tracking.sheltered(:, obj.Tracking.valid), 2) / obj.Video.FrameRate / 60 /60; % 
    Missing = TimeOutside == 0;
catch
    succ = false;
    warning('couldn''t compute profile: file not tracked');
    return;
end

try
    ContactList = obj.Hierarchy.Contacts.List; % List of all contacts
catch
    warning('no contacts list');
end

%% Initialize pairwise behavioral properties
try
    % Determine number of contacts between each pair:
    nc = hist3(sort([ContactList.subjects], 1)', 'Edges', {1:obj.nSubjects, 1:obj.nSubjects});
    nc = nc + nc'; 
    % Find mice that had no contacts
    NoContacts = all(nc == 0);
catch
    NoContacts = false(1, obj.nSubjects);
end
NoContacts = Q.tocol(NoContacts);

%% Find ROI indices of special regions
HighPlacesNames = opt.HighPlaces;
HighPlaces = [];
for i=1:length(HighPlacesNames)
    idx = find(strcmp(obj.ROI.ZoneNames, HighPlacesNames{i}));
    if ~isempty(idx)
        HighPlaces = [HighPlaces, idx];
    end
end

FoodWaterNames = opt.FoodWater;
FoodWater = [];
for i=1:length(FoodWaterNames)
    idx = find(strcmp(obj.ROI.ZoneNames, FoodWaterNames{i}));
    if ~isempty(idx)
        FoodWater = [FoodWater, idx];
    end
end

LocationNames = opt.Locations;
Locations = [];
for i=1:length(LocationNames)
    idx = find(strcmp(obj.ROI.ZoneNames, LocationNames{i}));
    if ~isempty(idx)
        Locations = [Locations, idx];
    end
end

%% Group Behaviors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for GroupSection = 1
    %% Video duration
    Profile.Group(1).Name = sprintf('Video duration [hours]');
    try
        Profile.Group(end).Val = sum(obj.Tracking.valid) / fps / 60/ 60;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% Number of hierarchy levels
    Profile.Group(end+1).Name = sprintf('Number of hierarchy levels');
    try
        Profile.Group(end).Val = length(unique(obj.Hierarchy.AggressiveChase.rank));
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% Hierarchy shape
    Profile.Group(end+1).Name = sprintf('Hierarchy shape (alpha-beta / beta-gamma)');
    try
        DS = sort(Hierarchy.DavidScore(obj.Hierarchy.AggressiveChase.ChaseEscape));
        Profile.Group(end).Val = (DS(end)-DS(end-1)) / (DS(end-1)-DS(end-2));
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% Hierarchy inbalance
    Profile.Group(end+1).Name = sprintf('Hierarchy inbalance mean(ds)/(alpha-zeta)');
    try
        DS = sort(Hierarchy.DavidScore(obj.Hierarchy.AggressiveChase.ChaseEscape));
        Profile.Group(end).Val = mean(DS) / (DS(end) - DS(1));
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% Hierarchy difference
    Profile.Group(end+1).Name = sprintf('Hierarchy difference (alpha-zeta)');
    try
        DS = sort(Hierarchy.DavidScore(obj.Hierarchy.AggressiveChase.ChaseEscape));
        Profile.Group(end).Val = (DS(end) - DS(1));
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% Number of contacts per hour
    Profile.Group(end+1).Name = sprintf('Number of contacts per hour');
    try
        nc = hist3(sort([ContactList.subjects], 1)', 'Edges', {1:obj.nSubjects, 1:obj.nSubjects});
        nc = nc + nc';
        
        Profile.Group(end).Val = mean(sum(nc, 2) / VideoDuration);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% Number of contacts per time outside
    Profile.Group(end+1).Name = sprintf('Number of contacts per time outside [hour-1]');
    try
        nc = hist3(sort([ContactList.subjects], 1)', 'Edges', {1:obj.nSubjects, 1:obj.nSubjects});
        nc = nc + nc';
        
        Profile.Group(end).Val = mean(sum(nc, 2) ./ TimeOutside);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% Number of chases per hour
    Profile.Group(end+1).Name = sprintf('Number of chases per hour');
    try
        p = sum(obj.Hierarchy.AggressiveChase.ChaseEscape, 2);
        Profile.Group(end).Val = mean(p / VideoDuration);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% Number of chases per time outside
    Profile.Group(end+1).Name = sprintf('Number of chases per time outside [hour-1]');
    try
        p = sum(obj.Hierarchy.AggressiveChase.ChaseEscape, 2);
        Profile.Group(end).Val = mean(p ./ TimeOutside);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% Chases per joint time
    Profile.Group(end+1).Name = sprintf('Chases per joint time [hour-1]');
    try
        %%
        ce = obj.Hierarchy.AggressiveChase.ChaseEscape;
        timeout = Hierarchy.TimeOutside(obj) / 3600;
        Profile.Group(end).Val = nanmean(ce ./ timeout, 2);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% Escapes per joint time
    Profile.Group(end+1).Name = sprintf('Escapes per joint time [hour-1]');
    try
        %%
        ce = obj.Hierarchy.AggressiveChase.ChaseEscape;
        timeout = Hierarchy.TimeOutside(obj) / 3600;
        Profile.Group(end).Val = nanmean(ce ./ timeout, 1)';
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% Follow per joint time
    Profile.Group(end+1).Name = sprintf('Follow per joint time [hour-1]');
    try
        %%
        ce = obj.Hierarchy.ChaseEscape.ChaseEscape;
        timeout = Hierarchy.TimeOutside(obj) / 3600;
        Profile.Group(end).Val = nanmean(ce ./ timeout, 2);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% NAChase per joint time
    Profile.Group(end+1).Name = sprintf('NA chase per joint time [hour-1]');
    try
        %%
        ce = obj.Hierarchy.ChaseEscape.ChaseEscape - obj.Hierarchy.AggressiveChase.ChaseEscape;
        timeout = Hierarchy.TimeOutside(obj) / 3600;
        Profile.Group(end).Val = nanmean(ce ./ timeout, 2);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% Being followed per joint time
    Profile.Group(end+1).Name = sprintf('Being followed per joint time [hour-1]');
    try
        %%
        ce = obj.Hierarchy.ChaseEscape.ChaseEscape;
        timeout = Hierarchy.TimeOutside(obj) / 3600;
        Profile.Group(end).Val = nanmean(ce ./ timeout, 1)';
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% NA-Escape per joint time
    Profile.Group(end+1).Name = sprintf('NA escape per joint time [hour-1]');
    try
        %%
        ce = obj.Hierarchy.ChaseEscape.ChaseEscape - obj.Hierarchy.AggressiveChase.ChaseEscape;
        timeout = Hierarchy.TimeOutside(obj) / 3600;
        Profile.Group(end).Val = nanmean(ce ./ timeout, 1)';
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    
    %% Approaches per joint time
    Profile.Group(end+1).Name = sprintf('Approaches per joint time [hour-1]');
    try
        %%
        ap = Hierarchy.Approaches(obj);
        timeout = Hierarchy.TimeOutside(obj) / 3600;
        Profile.Group(end).Val = nanmean(ap ./ timeout, 2);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% Being appraoched per joint time
    Profile.Group(end+1).Name = sprintf('Being approached per joint time [hour-1]');
    try
        %%
        ap = Hierarchy.Approaches(obj);
        timeout = Hierarchy.TimeOutside(obj) / 3600;
        Profile.Group(end).Val = nanmean(ap ./ timeout, 1)';
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    
    %% Number of approaches per hour
    Profile.Group(end+1).Name = sprintf('Number of approaches per hour');
    try
        c = zeros(obj.nSubjects);
        for l=1:length(ContactList)
            curr = ContactList(l);
            s1 = curr.subjects(1);
            s2 = curr.subjects(2);
            c(s1, s2) = c(s1, s2) + curr.states(1, 2);
            c(s2, s1) = c(s2, s1) + curr.states(2, 2);
        end
        p = sum(c, 2);
        Profile.Group(end).Val = mean(p / VideoDuration);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% Number of approaches per time outside
    Profile.Group(end+1).Name = sprintf('Number of approaches per time outside [hour-1]');
    try
        c = zeros(obj.nSubjects);
        for l=1:length(ContactList)
            curr = ContactList(l);
            s1 = curr.subjects(1);
            s2 = curr.subjects(2);
            c(s1, s2) = c(s1, s2) + curr.states(1, 2);
            c(s2, s1) = c(s2, s1) + curr.states(2, 2);
        end
        p = sum(c, 2);
        Profile.Group(end).Val = mean(p ./ TimeOutside);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    
    %% Entropy
    Profile.Group(end+1).Name = sprintf('Entropy');
    try
        p = zeros(obj.nSubjects, obj.ROI.nZones);
        for s=1:obj.nSubjects
            p(s, :) = p(s, :) + histc(obj.Tracking.zones(s, obj.Tracking.valid), 1:obj.ROI.nZones);
        end
        H = zeros(1, obj.nSubjects);
        for s=1:obj.nSubjects
            h = p(s, :) / sum(p(s, :));
            %Profile.Individual(end).Val(s, 1) = -h * log2(h' + (h' == 0)) / VideoDuration;
            H(s) = -h * log2(h' + (h' == 0));
        end
        Profile.Group(end).Val = mean(H);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% entropy outside
    Profile.Group(end+1).Name = 'Entropy outside';
    try
        p = zeros(obj.nSubjects, obj.ROI.nZones);
        for s=1:obj.nSubjects
            p(s, :) = p(s, :) + histc(obj.Tracking.zones(s, obj.Tracking.valid & ~obj.Tracking.sheltered(s, 1:length(obj.Tracking.valid))), 1:obj.ROI.nZones);
        end
        H = zeros(1, obj.nSubjects);
        for s=1:obj.nSubjects
            h = p(s, :) / sum(p(s, :));
            %Profile.Individual(end).Val(s, 1) = -h * log2(h' + (h' == 0)) / VideoDuration;
            H(s) = -h * log2(h' + (h' == 0));
        end
        Profile.Group(end).Val = mean(H);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    
    
    %% In-out Multi-information
    try
        data = double(~obj.Tracking.sheltered);
        iomodel = MaxEntropyGeneral.Train(data', 2, 2, 'Verbose', false);
        iomodel3 = MaxEntropyGeneral.Train(data', 3, 2, 'Verbose', false);
        %iomodel4 = MaxEntropyGeneral.Train(data', 4, 2, 'Verbose', false);
        Profile.Group(end+1).Name = 'InOut Multi-information';
        Profile.Group(end).Val = iomodel.EmpMultiInfo;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% In-out 2nd explained
    Profile.Group(end+1).Name = 'InOut % explained pairwise';
    try
        Profile.Group(end).Val = (iomodel.EmpIndepEntropy - iomodel.Entropy) / iomodel.EmpMultiInfo;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% In-out 3rd explained
    Profile.Group(end+1).Name = 'InOut % explained triplets';
    try
        Profile.Group(end).Val = (iomodel3.EmpIndepEntropy - iomodel3.Entropy) / iomodel3.EmpMultiInfo;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% Time outside
    Profile.Group(end+1).Name = 'Average percentage of time outside of nests';
    try
        Profile.Group(end).Val = mean(sum(~obj.Tracking.sheltered(:, obj.Tracking.valid), 2)) / TotalNumberOfFrames;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% % time high
    Profile.Group(end+1).Name = 'Average percentage of time at high place';
    try
        v = zeros(1, obj.nSubjects);
        for s=1:obj.nSubjects
            z = obj.Tracking.zones(s, obj.Tracking.valid);
            high = sum(ismember(z, HighPlaces));
            total = sum(obj.Tracking.valid);
            
            v(s) = high / total;
        end
        Profile.Group(end).Val = mean(v);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% % time high
    Profile.Group(end+1).Name = 'Average percentage of time at high place per time outside';
    try
        v = zeros(1, obj.nSubjects);
        for s=1:obj.nSubjects
            z = obj.Tracking.zones(s, obj.Tracking.valid);
            high = sum(ismember(z, HighPlaces));
            total = sum(~obj.Tracking.sheltered(s, :));
            
            v(s) = high / total;
        end
        Profile.Group(end).Val = mean(v);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
    %% frequ of entering high-place per time outside
    Profile.Group(end+1).Name = 'Freq of Entering High place per time outside';
    try
        %%
        v = zeros(obj.nSubjects, 1);
        for s=1:obj.nSubjects
            z = obj.Tracking.zones(s, obj.Tracking.valid);
            segs = Segs(ismember(z, HighPlaces));
            segs = segs.Close(round(opt.MinTimeInZone * fps));
            v(s) = length(segs.Events);
        end
        %%
        Profile.Group(end).Val = v ./ TimeOutside;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Group(end).Name ''' - %s\n']);
        Profile.Group(end).Val = nan;
    end
    
end

%% Individual Behaviors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for IndividualSection = 1
    %% Social rank
    Profile.Individual(1).Name = 'Social rank (1 day) [1 dominant - 4 submissive]';
    try
        r = max(obj.Hierarchy.AggressiveChase.rank) - obj.Hierarchy.AggressiveChase.rank + 1;
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = r(s);
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% David Score
    Profile.Individual(end+1).Name = 'Normalized David score (1 day) [higher = more dominant]';
    try
        NormDSonDij = Hierarchy.DavidScore(obj.Hierarchy.AggressiveChase.ChaseEscape);
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = NormDSonDij(s);
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% David Score - NA+Agg
    Profile.Individual(end+1).Name = 'Follow Normalized David score (1 day) [higher = more dominant]';
    try
        NormDSonDij = Hierarchy.DavidScore(obj.Hierarchy.ChaseEscape.ChaseEscape);
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = NormDSonDij(s);
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% David Score - NA+Agg
    Profile.Individual(end+1).Name = 'Approach Normalized David score (1 day) [higher = more dominant]';
    try
        NormDSonDij = Hierarchy.DavidScore(Hierarchy.Approaches(obj));
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = NormDSonDij(s);
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% David Score rank
    Profile.Individual(end+1).Name = 'David score rank (1 day) [1 dominant - 4 submissive]';
    try
        NormDSonDij = Hierarchy.DavidScore(obj.Hierarchy.AggressiveChase.ChaseEscape);
        Profile.Individual(end).Val = Q.rank(NormDSonDij);
        Profile.Individual(end).Val = max(Profile.Individual(end).Val) - Profile.Individual(end).Val + 1;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Classical David Score (not normalized)
    Profile.Individual(end+1).Name = 'Classical David score (1 day) [higher = more dominant, not normalized]';
    try
        [~, ~, DS] = Hierarchy.DavidScore(obj.Hierarchy.AggressiveChase.ChaseEscape);
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = DS(s);
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Time outside
    Profile.Individual(end+1).Name = 'Fraction of time outside';
    try
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = sum(obj.Tracking.valid & ~obj.Tracking.sheltered(s, 1:length(obj.Tracking.valid))) / TotalNumberOfFrames;
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end

    %% Time outside entropy
    Profile.Individual(end+1).Name = 'Time outside entropy';
    try
        nsegs = 12;
        win = obj.Video.FrameRate * 60 * 60;
        c = zeros(obj.nSubjects, nsegs);
        for i=1:nsegs
            r = round((i-1) * win + 1):round(i*win);
            r(r > length(obj.Tracking.valid)) = [];
            c(:, i) = sum(~obj.Tracking.sheltered(:, r), 2) ./ sum(obj.Tracking.valid(r), 2);
        end
        c = bsxfun(@rdivide, c, nansum(c, 2));
        %%
        for s=1:obj.nSubjects
            curr = c(s, :);
            curr(isnan(curr)) = [];
            Profile.Individual(end).Val(s, 1) = InfoTheory.Entropy(curr);
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end

    %% Start of day preference
    Profile.Individual(end+1).Name = 'Start of day preference';
    try
        nsegs = 3;
        win = obj.Video.FrameRate * 60 * 60 * 4;
        c = zeros(obj.nSubjects, nsegs);
        for i=1:nsegs
            r = round((i-1) * win + 1):round(i*win);
            r(r > length(obj.Tracking.valid)) = [];
            c(:, i) = sum(~obj.Tracking.sheltered(:, r), 2) ./ sum(obj.Tracking.valid(r), 2);
        end
        c = bsxfun(@rdivide, c, nansum(c, 2));
        %%
        Profile.Individual(end).Val(:, 1) = c(:, 1);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end

    %% End of day preference
    Profile.Individual(end+1).Name = 'End of day preference';
    try
        nsegs = 3;
        win = obj.Video.FrameRate * 60 * 60 * 4;
        c = zeros(obj.nSubjects, nsegs);
        for i=1:nsegs
            r = round((i-1) * win + 1):round(i*win);
            r(r > length(obj.Tracking.valid)) = [];
            c(:, i) = sum(~obj.Tracking.sheltered(:, r), 2) ./ sum(obj.Tracking.valid(r), 2);
        end
        c = bsxfun(@rdivide, c, nansum(c, 2));
        %%
        Profile.Individual(end).Val(:, 1) = c(:, 3);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Start of day approach preference
    Profile.Individual(end+1).Name = 'Start of day approach preference';
    try
        nsegs = 3;
        win = obj.Video.FrameRate * 60 * 60 * 4;
        c = zeros(obj.nSubjects, nsegs);
        for i=1:nsegs
            r = round((i-1) * win + 1):round(i*win);
            r(r > length(obj.Tracking.valid)) = [];
            c(:, i) = sum(Hierarchy.Approaches(obj, [r(1), r(end)]), 2);
        end
        c = bsxfun(@rdivide, c, nansum(c, 2));
        %%
        Profile.Individual(end).Val(:, 1) = c(:, 1);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% End of day approach preference
    Profile.Individual(end+1).Name = 'End of day approach preference';
    try
        nsegs = 3;
        win = obj.Video.FrameRate * 60 * 60 * 4;
        c = zeros(obj.nSubjects, nsegs);
        for i=1:nsegs
            r = round((i-1) * win + 1):round(i*win);
            r(r > length(obj.Tracking.valid)) = [];
            c(:, i) = sum(Hierarchy.Approaches(obj, [r(1), r(end)]), 2);
        end
        c = bsxfun(@rdivide, c, nansum(c, 2));
        %%
        Profile.Individual(end).Val(:, 1) = c(:, 3);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Number of in\out of nest
    Profile.Individual(end+1).Name = 'Visits outside rate';
    try
        Profile.Individual(end).Val = sum(diff(obj.Tracking.sheltered(:, obj.Tracking.valid), 1, 2) < 0, 2) / VideoDuration;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Foraging correlation
    Profile.Individual(end+1).Name = 'Foraging correlation';
    try
        Profile.Individual(end).Val = sum(corr(obj.Tracking.sheltered(:, obj.Tracking.valid)') .* (1-eye(obj.nSubjects)))' / (obj.nSubjects - 1);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Elo rating
    Profile.Individual(end+1).Name = 'Elo rating (end of day)';
    try
        Elo = Hierarchy.EloRating(obj);
        Profile.Individual(end).Val = Elo.Scores(:, end);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% No. of contacts
    Profile.Individual(end+1).Name = 'Contact rate [1/hour]';
    try
        nc = hist3(sort([ContactList.subjects], 1)', 'Edges', {1:obj.nSubjects, 1:obj.nSubjects});
        nc = nc + nc';
        
        Profile.Individual(end).Val = sum(nc, 2) / VideoDuration;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end

    %% Mean contact rate interval
    Profile.Individual(end+1).Name = 'Mean contact rate interval';
    try
        subj = [ContactList.subjects];
        begf = [ContactList.beg];
        begf = min(begf);
        endf = [ContactList.end];
        endf = max(endf);
        for s=1:obj.nSubjects
            m = any(subj == s);
            b = begf(m);
            e = endf(m);
            len = b(2:end) - e(1:end-1) + 1;
            Profile.Individual(end).Val(s, 1) = mean(len);
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Contact rate
    Profile.Individual(end+1).Name = 'Contact rate outside [1/hour]';
    try
        nc = hist3(sort([ContactList.subjects], 1)', 'Edges', {1:obj.nSubjects, 1:obj.nSubjects});
        nc = nc + nc';
        
        Profile.Individual(end).Val = sum(nc, 2) ./ TimeOutside;
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% % time in contact
    Profile.Individual(end+1).Name = 'Fraction of time in contact outside [1/hour]';
    try
        subj = [ContactList.subjects];
        begf = [ContactList.beg];
        begf = min(begf);
        endf = [ContactList.end];
        endf = max(endf);
        for s=1:obj.nSubjects
            m = any(subj == s);
            len = endf(m) - begf(m) + 1;
            Profile.Individual(end).Val(s, 1) = sum(len) / sum(~obj.Tracking.sheltered(s, obj.Tracking.valid), 2);
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% median duration of contact
    Profile.Individual(end+1).Name = 'Median contact duration [secs]';
    try
        subj = [ContactList.subjects];
        begf = [ContactList.beg];
        begf = min(begf);
        endf = [ContactList.end];
        endf = max(endf);
        for s=1:obj.nSubjects
            m = any(subj == s);
            len = endf(m) - begf(m) + 1;
            Profile.Individual(end).Val(s, 1) = median(len) / fps;
        end
        Profile.Individual(end).Val(Missing | NoContacts) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Mean duration of contact
    Profile.Individual(end+1).Name = 'Mean contact duration [secs]';
    try
        subj = [ContactList.subjects];
        begf = [ContactList.beg];
        begf = min(begf);
        endf = [ContactList.end];
        endf = max(endf);
        for s=1:obj.nSubjects
            m = any(subj == s);
            len = endf(m) - begf(m) + 1;
            Profile.Individual(end).Val(s, 1) = mean(len) / fps;
        end
        Profile.Individual(end).Val(Missing | NoContacts) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end

    %% Contact duration pareto
    Profile.Individual(end+1).Name = 'Contact duration pareto';
    try
        subj = [ContactList.subjects];
        begf = [ContactList.beg];
        begf = min(begf);
        endf = [ContactList.end];
        endf = max(endf);
        for s=1:obj.nSubjects
            m = any(subj == s);
            len = sort(endf(m) - begf(m) + 1, 'descend');
            try
                Profile.Individual(end).Val(s, 1) = (find(cumsum(len) >= .8 * sum(len), 1) - 1) / (length(len) - 1);
            catch
                Profile.Individual(end).Val(s, 1) = 1;
            end
        end
        Profile.Individual(end).Val(Missing | NoContacts) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end

    %% Location duration pareto
    Profile.Individual(end+1).Name = 'Location duration pareto';
    try
        for s=1:obj.nSubjects
            z = obj.Tracking.zones(s, obj.Tracking.valid);
            map = ismember(z, Locations);
            z = z .* uint16(map);
            %%
            ez = [z, 0];
            d = diff([obj.ROI.nZones+1 z obj.ROI.nZones+1]) ~= 0;
            idx = find(d);
            dd = diff(idx);
            durations = sort(dd(ez(idx(1:end-1)) ~= 0), 'descend');
            
            %%
            try
                Profile.Individual(end).Val(s, 1) = (find(cumsum(durations) >= .8 * sum(durations), 1) - 1) / (length(durations) - 1);
            catch
                Profile.Individual(end).Val(s, 1) = 1;
            end
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end    
    
    %% Approach/chase diff
    Profile.Individual(end+1).Name = 'Diff between approaches and chases';
    try
        c = zeros(obj.nSubjects);
        for l=1:length(ContactList)
            curr = ContactList(l);
            s1 = curr.subjects(1);
            s2 = curr.subjects(2);
            c(s1, s2) = c(s1, s2) + curr.states(1, 2);
            c(s2, s1) = c(s2, s1) + curr.states(2, 2);
        end
        Profile.Individual(end).Val(:, 1) = sum(c, 2) - sum(obj.Hierarchy.AggressiveChase.ChaseEscape, 2);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Approach/escape ratio
    Profile.Individual(end+1).Name = 'approach and escape';
    try
        %%
        c = zeros(obj.nSubjects, 1);
        n = zeros(obj.nSubjects, 1);
        for l=1:length(ContactList)
            curr = ContactList(l);
            s1 = curr.subjects(1);
            s2 = curr.subjects(2);
            if l <= length(obj.Hierarchy.Contacts.Behaviors.AggressiveChase.Escaper)
                c(s1) = c(s1) + curr.states(1, 2) & obj.Hierarchy.Contacts.Behaviors.AggressiveChase.Escaper(l) == s1;
                c(s2) = c(s2) + curr.states(2, 2) & obj.Hierarchy.Contacts.Behaviors.AggressiveChase.Escaper(l) == s2;
            end
            n(s1) = n(s1) + curr.states(1, 2);
            n(s2) = n(s2) + curr.states(2, 2);
        end
        %%
        Profile.Individual(end).Val(:, 1) = c./max(n, 1);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Approach/chase ratio
    Profile.Individual(end+1).Name = 'approach and chase';
    try
        %%
        c = zeros(obj.nSubjects, 1);
        n = zeros(obj.nSubjects, 1);
        for l=1:length(ContactList)
            curr = ContactList(l);
            s1 = curr.subjects(1);
            s2 = curr.subjects(2);
            if l <= length(obj.Hierarchy.Contacts.Behaviors.AggressiveChase.Chaser)
                c(s1) = c(s1) + curr.states(1, 2) & obj.Hierarchy.Contacts.Behaviors.AggressiveChase.Chaser(l) == s1;
                c(s2) = c(s2) + curr.states(2, 2) & obj.Hierarchy.Contacts.Behaviors.AggressiveChase.Chaser(l) == s2;
            end
            n(s1) = n(s1) + curr.states(1, 2);
            n(s2) = n(s2) + curr.states(2, 2);
        end
        %%
        Profile.Individual(end).Val(:, 1) = c./max(n, 1);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Fraction of chases per contact
    Profile.Individual(end+1).Name = 'Fraction of chases per contact';
    try
        nc = hist3(sort([ContactList.subjects], 1)', 'Edges', {1:obj.nSubjects, 1:obj.nSubjects});
        nc = nc + nc';
        
        p = sum(obj.Hierarchy.AggressiveChase.ChaseEscape, 2);
        Profile.Individual(end).Val(:, 1) = p ./ sum(nc, 2);
        Profile.Individual(end).Val(Missing | NoContacts) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Number of escapes
    Profile.Individual(end+1).Name = 'Fraction of escapes per contact';
    try
        nc = hist3(sort([ContactList.subjects], 1)', 'Edges', {1:obj.nSubjects, 1:obj.nSubjects});
        nc = nc + nc';
        
        p = sum(obj.Hierarchy.AggressiveChase.ChaseEscape, 1)';
        Profile.Individual(end).Val(:, 1) = p ./ sum(nc, 2);
        Profile.Individual(end).Val(Missing | NoContacts) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Fraction of follows per contact
    Profile.Individual(end+1).Name = 'Fraction of follows per contact';
    try
        nc = hist3(sort([ContactList.subjects], 1)', 'Edges', {1:obj.nSubjects, 1:obj.nSubjects});
        nc = nc + nc';
        
        p = sum(obj.Hierarchy.ChaseEscape.ChaseEscape, 2);
        Profile.Individual(end).Val(:, 1) = p ./ sum(nc, 2);
        Profile.Individual(end).Val(Missing | NoContacts) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Fraction of being followed per contact
    Profile.Individual(end+1).Name = 'Fraction of being followed per contact';
    try
        nc = hist3(sort([ContactList.subjects], 1)', 'Edges', {1:obj.nSubjects, 1:obj.nSubjects});
        nc = nc + nc';
        
        p = sum(obj.Hierarchy.ChaseEscape.ChaseEscape, 1)';
        Profile.Individual(end).Val(:, 1) = p ./ sum(nc, 2);
        Profile.Individual(end).Val(Missing | NoContacts) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Fraction of NA chases per contact
    Profile.Individual(end+1).Name = 'Fraction of NA chases per contact';
    try
        nc = hist3(sort([ContactList.subjects], 1)', 'Edges', {1:obj.nSubjects, 1:obj.nSubjects});
        nc = nc + nc';
        
        p = sum(obj.Hierarchy.ChaseEscape.ChaseEscape - obj.Hierarchy.AggressiveChase.ChaseEscape, 2);
        Profile.Individual(end).Val(:, 1) = p ./ sum(nc, 2);
        Profile.Individual(end).Val(Missing | NoContacts) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Fraction of NA escapes per contact
    Profile.Individual(end+1).Name = 'Fraction of NA escapes per contact';
    try
        nc = hist3(sort([ContactList.subjects], 1)', 'Edges', {1:obj.nSubjects, 1:obj.nSubjects});
        nc = nc + nc';
        
        p = sum(obj.Hierarchy.ChaseEscape.ChaseEscape - obj.Hierarchy.AggressiveChase.ChaseEscape, 1)';
        Profile.Individual(end).Val(:, 1) = p ./ sum(nc, 2);
        Profile.Individual(end).Val(Missing | NoContacts) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Chase rate outside
    Profile.Individual(end+1).Name = 'Aggressive chase rate outside [1/hour]';
    try
        p = sum(obj.Hierarchy.AggressiveChase.ChaseEscape, 2);
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = p(s) / TimeOutside(s);
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Escape rate outside
    Profile.Individual(end+1).Name = 'Aggressive escape rate outside [1/hour]';
    try
        p = sum(obj.Hierarchy.AggressiveChase.ChaseEscape, 1);
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = p(s) / TimeOutside(s);
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Follow rate outside
    Profile.Individual(end+1).Name = 'Follow rate outside [1/hour]';
    try
        p = sum(obj.Hierarchy.ChaseEscape.ChaseEscape, 2);
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = p(s) / TimeOutside(s);
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Being followed rate outside
    Profile.Individual(end+1).Name = 'Being followed rate outside [1/hour]';
    try
        p = sum(obj.Hierarchy.ChaseEscape.ChaseEscape, 1);
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = p(s) / TimeOutside(s);
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% NA chase rate outside
    Profile.Individual(end+1).Name = 'NA chase rate outside [1/hour]';
    try
        p = sum(obj.Hierarchy.ChaseEscape.ChaseEscape - obj.Hierarchy.AggressiveChase.ChaseEscape, 2);
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = p(s) / TimeOutside(s);
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% NA escape rate outside
    Profile.Individual(end+1).Name = 'NA escape rate outside [1/hour]';
    try
        p = sum(obj.Hierarchy.ChaseEscape.ChaseEscape - obj.Hierarchy.AggressiveChase.ChaseEscape, 1);
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = p(s) / TimeOutside(s);
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Aggressive chase rate
    Profile.Individual(end+1).Name = 'Aggressive chase rate [1/hour]';
    try
        p = sum(obj.Hierarchy.AggressiveChase.ChaseEscape, 2);
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = p(s) / VideoDuration;
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Aggressive escape rate
    Profile.Individual(end+1).Name = 'Aggressive escape rate [1/hour]';
    try
        p = sum(obj.Hierarchy.AggressiveChase.ChaseEscape, 1);
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = p(s) / VideoDuration;
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Follow rate
    Profile.Individual(end+1).Name = 'Follow rate [1/hour]';
    try
        p = sum(obj.Hierarchy.ChaseEscape.ChaseEscape, 2);
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = p(s) / VideoDuration;
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Being followed rate
    Profile.Individual(end+1).Name = 'Being followed rate [1/hour]';
    try
        p = sum(obj.Hierarchy.ChaseEscape.ChaseEscape, 1);
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = p(s) / VideoDuration;
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% NA chase rate
    Profile.Individual(end+1).Name = 'NA chase rate [1/hour]';
    try
        p = sum(obj.Hierarchy.ChaseEscape.ChaseEscape - obj.Hierarchy.AggressiveChase.ChaseEscape, 2);
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = p(s) / VideoDuration;
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% NA escape rate
    Profile.Individual(end+1).Name = 'NA escape rate [1/hour]';
    try
        p = sum(obj.Hierarchy.ChaseEscape.ChaseEscape - obj.Hierarchy.AggressiveChase.ChaseEscape, 1);
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = p(s) / VideoDuration;
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% approach
    Profile.Individual(end+1).Name = 'Number of Approachs';
    try
        c = zeros(obj.nSubjects);
        for l=1:length(ContactList)
            curr = ContactList(l);
            s1 = curr.subjects(1);
            s2 = curr.subjects(2);
            c(s1, s2) = c(s1, s2) + curr.states(1, 2);
            c(s2, s1) = c(s2, s1) + curr.states(2, 2);
        end
        p = sum(c, 2);
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = p(s);
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Number of approaches per time outside
    Profile.Individual(end+1).Name = sprintf('Approach rate outside [1/hour]');
    try
        c = zeros(obj.nSubjects);
        for l=1:length(ContactList)
            curr = ContactList(l);
            s1 = curr.subjects(1);
            s2 = curr.subjects(2);
            c(s1, s2) = c(s1, s2) + curr.states(1, 2);
            c(s2, s1) = c(s2, s1) + curr.states(2, 2);
        end
        p = sum(c, 2);
        Profile.Individual(end).Val = p ./ TimeOutside;
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% approach
    Profile.Individual(end+1).Name = 'Number of approaches per couple out';
    try
        c = zeros(obj.nSubjects);
        for l=1:length(ContactList)
            curr = ContactList(l);
            s1 = curr.subjects(1);
            s2 = curr.subjects(2);
            c(s1, s2) = c(s1, s2) + curr.states(1, 2);
            c(s2, s1) = c(s2, s1) + curr.states(2, 2);
            
        end
        p = sum(c, 2);
        for s=1:obj.nSubjects
            t = sum(sum(~obj.Tracking.sheltered(:, obj.Tracking.valid & ~obj.Tracking.sheltered(s, 1:length(obj.Tracking.valid)))) == 2);
            Profile.Individual(end).Val(s, 1) = p(s) / (t / fps / 60 /60);
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% approach
    Profile.Individual(end+1).Name = 'Number of approaches per mice out (two or more)';
    try
        c = zeros(obj.nSubjects);
        for l=1:length(ContactList)
            curr = ContactList(l);
            s1 = curr.subjects(1);
            s2 = curr.subjects(2);
            c(s1, s2) = c(s1, s2) + curr.states(1, 2);
            c(s2, s1) = c(s2, s1) + curr.states(2, 2);
            
        end
        p = sum(c, 2);
        for s=1:obj.nSubjects
            t = sum(sum(~obj.Tracking.sheltered(:, obj.Tracking.valid & ~obj.Tracking.sheltered(s, 1:length(obj.Tracking.valid)))) >= 2);
            Profile.Individual(end).Val(s, 1) = p(s) / (t / fps / 60 /60);
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Number of approaches per contact
    Profile.Individual(end+1).Name = 'Fraction of approaches per contact';
    try
        nc = hist3(sort([ContactList.subjects], 1)', 'Edges', {1:obj.nSubjects, 1:obj.nSubjects});
        nc = nc + nc';
        
        c = zeros(obj.nSubjects);
        for l=1:length(ContactList)
            curr = ContactList(l);
            s1 = curr.subjects(1);
            s2 = curr.subjects(2);
            c(s1, s2) = c(s1, s2) + curr.states(1, 2);
            c(s2, s1) = c(s2, s1) + curr.states(2, 2);
        end
        Profile.Individual(end).Val = sum(c, 2) ./ sum(nc, 2);
        Profile.Individual(end).Val(Missing | NoContacts) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end

    %% Ratio between approaches and being approached
    Profile.Individual(end+1).Name = 'Ratio between approaches and being approached';
    try
        Profile.Individual(end).Val = sum(Hierarchy.Approaches(obj), 2) ./ sum(Hierarchy.Approaches(obj), 1)';
        Profile.Individual(end).Val(Missing | NoContacts) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Number of approaches per hour
    Profile.Individual(end+1).Name = 'Approach rate [1/hour]';
    try
        %nc = hist3(sort([ContactList.subjects], 1)', 'Edges', {1:obj.nSubjects, 1:obj.nSubjects});
        %nc = nc + nc';
        
        c = zeros(obj.nSubjects);
        for l=1:length(ContactList)
            curr = ContactList(l);
            s1 = curr.subjects(1);
            s2 = curr.subjects(2);
            c(s1, s2) = c(s1, s2) + curr.states(1, 2);
            c(s2, s1) = c(s2, s1) + curr.states(2, 2);
        end
        Profile.Individual(end).Val = sum(c, 2) / VideoDuration;
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% approach entropy
    Profile.Individual(end+1).Name = 'Approach entropy';
    try
        %%
        c = Hierarchy.Approaches(obj);
        c = bsxfun(@rdivide, c, sum(c, 2));
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = InfoTheory.Entropy(c(s, Q.exclude(obj.nSubjects, s)));
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Number of being approached per contact
    Profile.Individual(end+1).Name = 'Fraction of being approached per contact';
    try
        nc = hist3(sort([ContactList.subjects], 1)', 'Edges', {1:obj.nSubjects, 1:obj.nSubjects});
        nc = nc + nc';
        
        c = zeros(obj.nSubjects);
        for l=1:length(ContactList)
            curr = ContactList(l);
            s1 = curr.subjects(1);
            s2 = curr.subjects(2);
            c(s1, s2) = c(s1, s2) + curr.states(1, 2);
            c(s2, s1) = c(s2, s1) + curr.states(2, 2);
        end
        Profile.Individual(end).Val = sum(c, 1)' ./ sum(nc, 2);
        Profile.Individual(end).Val(Missing | NoContacts) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% being approach
    Profile.Individual(end+1).Name = 'Being approached rate outside [1/hour]';
    try
        c = zeros(obj.nSubjects);
        for l=1:length(ContactList)
            curr = ContactList(l);
            s1 = curr.subjects(1);
            s2 = curr.subjects(2);
            c(s1, s2) = c(s1, s2) + curr.states(1, 2);
            c(s2, s1) = c(s2, s1) + curr.states(2, 2);
        end
        p = sum(c, 1)';
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = p(s) / TimeOutside(s);
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% being approached rate
    Profile.Individual(end+1).Name = 'Being approached rate [1/hour]';
    try
        c = zeros(obj.nSubjects);
        for l=1:length(ContactList)
            curr = ContactList(l);
            s1 = curr.subjects(1);
            s2 = curr.subjects(2);
            c(s1, s2) = c(s1, s2) + curr.states(1, 2);
            c(s2, s1) = c(s2, s1) + curr.states(2, 2);
        end
        p = sum(c, 1)';
        for s=1:obj.nSubjects
            Profile.Individual(end).Val(s, 1) = p(s) / VideoDuration;
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    %% Number of Approach-Escape per contacts
    Profile.Individual(end+1).Name = 'Fraction of approach escape behavior per aggression';
    try
        c = zeros(obj.nSubjects, 1);
        chaser  = obj.Hierarchy.Contacts.Behaviors.AggressiveChase.Chaser;  chaser(end+1:length(ContactList)) = 0;
        escaper = obj.Hierarchy.Contacts.Behaviors.AggressiveChase.Escaper; escaper(end+1:length(ContactList)) = 0;
        for l=1:length(ContactList)
            curr = ContactList(l);
            s1 = curr.subjects(1);
            s2 = curr.subjects(2);
            %             c(s1, s2) = c(s1, s2) + curr.states(1, 2) & curr.states(1, 6);
            %             c(s2, s1) = c(s2, s1) + curr.states(2, 2) & curr.states(2, 6);
            c(s1) = c(s1) + double(curr.states(1, 2) & (escaper(l) == s1));
            c(s2) = c(s2) + double(curr.states(2, 2) & (escaper(l) == s2));
        end
        agg = zeros(obj.nSubjects, 1);
        for s=1:obj.nSubjects
            agg(s) = sum(chaser == s | escaper == s);
        end
        Profile.Individual(end).Val = c ./ agg;
        Profile.Individual(end).Val(isnan(Profile.Individual(end).Val)) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% entropy
    Profile.Individual(end+1).Name = 'Entropy [bits]';
    try
        p = zeros(obj.nSubjects, obj.ROI.nZones);
        for s=1:obj.nSubjects
            p(s, :) = p(s, :) + histc(obj.Tracking.zones(s, obj.Tracking.valid), 1:obj.ROI.nZones);
        end
        for s=1:obj.nSubjects
            h = p(s, :) / sum(p(s, :));
            %Profile.Individual(end).Val(s, 1) = -h * log2(h' + (h' == 0)) / VideoDuration;
            Profile.Individual(end).Val(s, 1) = -h * log2(h' + (h' == 0));
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% entropy outside
    Profile.Individual(end+1).Name = 'Entropy outside [bits]';
    try
        p = zeros(obj.nSubjects, obj.ROI.nZones);
        for s=1:obj.nSubjects
            p(s, :) = p(s, :) + histc(obj.Tracking.zones(s, obj.Tracking.valid & ~obj.Tracking.sheltered(s, 1:length(obj.Tracking.valid))), 1:obj.ROI.nZones);
        end
        for s=1:obj.nSubjects
            h = p(s, :) / sum(p(s, :));
            %Profile.Individual(end).Val(s, 1) = -h * log2(h' + (h' == 0)) / VideoDuration;
            Profile.Individual(end).Val(s, 1) = -h * log2(h' + (h' == 0));
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% % time near food or water
    Profile.Individual(end+1).Name = 'Fraction of time near food or water [%]';
    try
        for s=1:obj.nSubjects
            z = obj.Tracking.zones(s, obj.Tracking.valid);
            near = sum(ismember(z, FoodWater));
            total = sum(obj.Tracking.valid);
            Profile.Individual(end).Val(s, 1) = near / total;
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% % time near food or water
    Profile.Individual(end+1).Name = 'Food or water per time outside [%]';
    try
        for s=1:obj.nSubjects
            z = obj.Tracking.zones(s, obj.Tracking.valid);
            near = sum(ismember(z, FoodWater));
            %total = sum(obj.Tracking.valid);
            Profile.Individual(end).Val(s, 1) = near / sum(~obj.Tracking.sheltered(s, obj.Tracking.valid), 2);
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% % time near food
    Profile.Individual(end+1).Name = 'Fraction of time in Feeder outside [%]';
    try
        for s=1:obj.nSubjects
            z = obj.Tracking.zones(s, obj.Tracking.valid);
            near = sum(ismember(z, find(ismember(obj.ROI.ZoneNames, {'Feeder', 'Feeder1', 'Feeder2'}))), 2);
            %total = sum(obj.Tracking.valid);
            Profile.Individual(end).Val(s, 1) = near / sum(~obj.Tracking.sheltered(s, obj.Tracking.valid), 2);
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% % time near water
    Profile.Individual(end+1).Name = 'Fraction of time in Water outside [%]';
    try
        for s=1:obj.nSubjects
            z = obj.Tracking.zones(s, obj.Tracking.valid);
            near = sum(ismember(z, find(ismember(obj.ROI.ZoneNames, {'Water', 'Water1', 'Water2'}))), 2);
            %total = sum(obj.Tracking.valid);
            Profile.Individual(end).Val(s, 1) = near / sum(~obj.Tracking.sheltered(s, obj.Tracking.valid), 2);
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Proximate vs distant food
    Profile.Individual(end+1).Name = 'Proximate vs distant food';
    try
        %%
        z = obj.Tracking.zones(:, obj.Tracking.valid);
        f1 = mean(ismember(z, find(ismember(obj.ROI.ZoneNames, {'Feeder', 'Feeder1'}))), 2);
        f2 = mean(ismember(z, find(ismember(obj.ROI.ZoneNames, {'Feeder2'}))), 2);
        %%
        Profile.Individual(end).Val = f1 ./ (f1 + f2);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Proximate vs distant water
    Profile.Individual(end+1).Name = 'Proximate vs distant water';
    try
        %%
        z = obj.Tracking.zones(:, obj.Tracking.valid);
        f1 = mean(ismember(z, find(ismember(obj.ROI.ZoneNames, {'Water', 'Water1'}))), 2);
        f2 = mean(ismember(z, find(ismember(obj.ROI.ZoneNames, {'Water2'}))), 2);
        %%
        Profile.Individual(end).Val = f1 ./ (f1 + f2);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    
    %% % time high
    Profile.Individual(end+1).Name = 'Fraction of time at high place [%]';
    try
        for s=1:obj.nSubjects
            z = obj.Tracking.zones(s, obj.Tracking.valid);
            high = sum(ismember(z, HighPlaces));
            total = sum(obj.Tracking.valid);
            
            Profile.Individual(end).Val(s, 1) = high / total;
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    
    %% % time high
    Profile.Individual(end+1).Name = 'High place per time outside [1/hour]';
    try
        for s=1:obj.nSubjects
            z = obj.Tracking.zones(s, obj.Tracking.valid);
            high = sum(ismember(z, HighPlaces));
            %total = sum(obj.Tracking.valid & ~obj.Tracking.sheltered(s, :));
            
            Profile.Individual(end).Val(s, 1) = high ./ sum(~obj.Tracking.sheltered(s, obj.Tracking.valid), 2);
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Fraction of time in open outside
    Profile.Individual(end+1).Name = 'Fraction of time in the open outside [%]';
    try
        %%
        z = obj.Tracking.zones(:, obj.Tracking.valid);
        openzone = sum(z == find(ismember(obj.ROI.ZoneNames, {'Open'})), 2) / fps / 3600;
        %%
        Profile.Individual(end).Val = openzone ./ TimeOutside;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan;
    end
    
    %% Fraction of time on small nest outside
    Profile.Individual(end+1).Name = 'Fraction of time in SmallNest outside [%]';
    try
        %%
        z = obj.Tracking.zones(:, obj.Tracking.valid);
        openzone = sum(z == find(ismember(obj.ROI.ZoneNames, {'SmallNest'})), 2) / fps / 3600;
        %%
        Profile.Individual(end).Val = openzone ./ TimeOutside;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Fraction of time on ramp outside
    Profile.Individual(end+1).Name = 'Fraction of time on Ramp outside [%]';
    try
        for s=1:obj.nSubjects
            z = obj.Tracking.zones(s, obj.Tracking.valid);
            near = sum(ismember(z, find(ismember(obj.ROI.ZoneNames, {'[Ramp1]', '[Ramp2]', 'Ramp1', 'Ramp2'}))), 2);
            %total = sum(obj.Tracking.valid);
            Profile.Individual(end).Val(s, 1) = near / sum(~obj.Tracking.sheltered(s, obj.Tracking.valid), 2);
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Fraction of time in Labyrinth outside
    Profile.Individual(end+1).Name = 'Fraction of time in Labyrinth outside [%]';
    try
        %%
        z = obj.Tracking.zones(:, obj.Tracking.valid);
        openzone = sum(z == find(ismember(obj.ROI.ZoneNames, {'Labyrinth'})), 2) / fps / 3600;
        %%
        Profile.Individual(end).Val = openzone ./ TimeOutside;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan;
    end
    
    %% Mean distance from walls
    Profile.Individual(end+1).Name = 'Distance from walls in open';
    try
        %%
        ac = obj.Meta.Scale.ArenaCoord;
        wall_n = ac(1):1:ac(1)+ac(3); wall_n(2, :) = ac(2);
        wall_s = ac(1):1:ac(1)+ac(3); wall_s(2, :) = ac(2)+ac(4);
        wall_e(2, :) = ac(2):1:ac(2)+ac(4); wall_e(1, :) = ac(1);
        wall_w(2, :) = ac(2):1:ac(2)+ac(4); wall_w(1, :) = ac(1)+ac(3);
        wall = [wall_n, wall_s, wall_e, wall_w]';
        
        %%
        map = obj.Tracking.zones == find(ismember(obj.ROI.ZoneNames, {'Open'}));
        %%
        for s=1:obj.nSubjects
            %%
            x = obj.Tracking.x(s, map(s, :));
            y = obj.Tracking.y(s, map(s, :));
            Profile.Individual(end).Val(s, 1) = mean(sqrt(sum((wall(knnsearch(wall, [x(:), y(:)]), :) - [x(:), y(:)]).^2, 2))) / PixelsPerCM;
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan;
    end
    
    
    %% Mean distance from nest
    Profile.Individual(end+1).Name = 'Distance from nest';
    try
        %%
        nests = cumsum(cellfun(@(x) ~isempty(x), regexp(obj.ROI.ZoneNames, '\(.*\)')));
        nestid = nests(strcmp(obj.ROI.ZoneNames, '(BigNest)'));
        bounds = bwboundaries(obj.ROI.Hidden{nestid});
        bounds = bounds{1};
        %%
        for s=1:obj.nSubjects
            %%
            x = obj.Tracking.x(s, obj.Tracking.valid & ~obj.Tracking.sheltered(s, :));
            y = obj.Tracking.y(s, obj.Tracking.valid & ~obj.Tracking.sheltered(s, :));
            Profile.Individual(end).Val(s, 1) = mean(sqrt(sum((bounds(knnsearch(bounds, [x(:), y(:)]), :) - [x(:), y(:)]).^2, 2))) / PixelsPerCM;
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% % time alone outside
    Profile.Individual(end+1).Name = 'Fraction of time alone outside';
    try
        for s=1:obj.nSubjects
            other = 1:obj.nSubjects;
            other = other(other ~= s);
            p = sum(all(obj.Tracking.sheltered(other, obj.Tracking.valid)) & ~obj.Tracking.sheltered(s, obj.Tracking.valid));
            Profile.Individual(end).Val(s, 1) = p / sum(~obj.Tracking.sheltered(s, obj.Tracking.valid), 2);
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end

    %% % time alone in shelter
    Profile.Individual(end+1).Name = 'Fraction of time alone in shelter';
    try
        for s=1:obj.nSubjects
            other = 1:obj.nSubjects;
            other = other(other ~= s);
            p = sum(all(~obj.Tracking.sheltered(other, obj.Tracking.valid)) & obj.Tracking.sheltered(s, obj.Tracking.valid));
            Profile.Individual(end).Val(s, 1) = p / sum(obj.Tracking.sheltered(s, obj.Tracking.valid), 2);
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end

    %% mean number of mice outside
    Profile.Individual(end+1).Name = 'Mean number of other mice outside together';
    try
        for s=1:obj.nSubjects
            other = Q.exclude(obj.nSubjects, s);
            ns = ~obj.Tracking.sheltered(:, obj.Tracking.valid);
            p = sum(ns(other, ns(s, :)));
            Profile.Individual(end).Val(s, 1) = mean(p);
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end

    %% mean number of mice outside
    Profile.Individual(end+1).Name = 'Mean number of other mice sheltered together';
    try
        for s=1:obj.nSubjects
            other = Q.exclude(obj.nSubjects, s);
            ins = ~obj.Tracking.sheltered(:, obj.Tracking.valid);
            p = sum(ins(other, ins(s, :)));
            Profile.Individual(end).Val(s, 1) = mean(p);
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% median speed
    Profile.Individual(end+1).Name = 'Median speed outside (rough) [m/sec]';
    Profile.Group(end+1).Name = 'Average median speed (rough) [m/sec]';
    try
        jumpTime = 1; %% sec
        for s=1:obj.nSubjects
            x = Segs(~obj.Tracking.sheltered(s, 1:length(obj.Tracking.valid)) & obj.Tracking.valid, obj.Tracking.x(s, :));
            y = Segs(~obj.Tracking.sheltered(s, 1:length(obj.Tracking.valid)) & obj.Tracking.valid, obj.Tracking.y(s, :));
            dt = round(fps * jumpTime);
            speed = nan(1, x.Length);
            for i=1:length(x.Events)
                cx = x.Events(i).data;
                cx = cx(1:dt:end);
                cy = y.Events(i).data;
                cy = cy(1:dt:end);
                if length(cx) > 2
                    d = sqrt((cx(2:end)-cx(1:end-1)).^2 + (cy(2:end)-cy(1:end-1)).^2) / (PixelsPerCM * 100) / jumpTime;
                    speed(x.Events(i).beg:x.Events(i).beg+length(d)-1) = d;
                end
            end
            Profile.Individual(end).Val(s, 1) =  nanmedian(speed);
        end
        Profile.Group(end).Val = mean(Profile.Individual(end).Val);
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
        Profile.Group(end).Val = nan;
    end
    
    %% mean speed outside
    Profile.Individual(end+1).Name = 'Mean speed outside (rough) [m/sec]';
    try
        if Q.isfield(obj, 'Meta.Scale.PixPerCM')
            PixelsPerCM = obj.Meta.Scale.PixPerCM;
        else
            PixelsPerCM = mean([obj.Meta.Scale.ArenaCoord(3) / obj.Meta.Scale.ArenaWidth, obj.Meta.Scale.ArenaCoord(4) / obj.Meta.Scale.ArenaHeight]);
        end
        jumpTime = opt.DistanceNoiseDuration; %% sec
        for s=1:obj.nSubjects
            x = Segs(~obj.Tracking.sheltered(s, 1:length(obj.Tracking.valid)) & obj.Tracking.valid, obj.Tracking.x(s, :));
            y = Segs(~obj.Tracking.sheltered(s, 1:length(obj.Tracking.valid)) & obj.Tracking.valid, obj.Tracking.y(s, :));
            dt = round(fps * jumpTime);
            speed = nan(1, x.Length);
            for i=1:length(x.Events)
                cx = x.Events(i).data;
                cx = cx(1:dt:end);
                cy = y.Events(i).data;
                cy = cy(1:dt:end);
                if length(cx) > 2
                    d = sqrt((cx(2:end)-cx(1:end-1)).^2 + (cy(2:end)-cy(1:end-1)).^2) / PixelsPerCM;
                    d(d < opt.DistanceNoiseThreshold) = 0;
                    speed(x.Events(i).beg:x.Events(i).beg+length(d)-1) = d / 100 / jumpTime;
                end
            end
            Profile.Individual(end).Val(s, 1) =  nanmean(speed);
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end

    %% mean time in locations
    Profile.Individual(end+1).Name = 'Mean time in locations';
    try
        %%
        for s=1:obj.nSubjects
            z = obj.Tracking.zones(s, obj.Tracking.valid);
            map = ismember(z, Locations);
            z = z .* uint16(map);
            %%
            ez = [z, 0];
            d = diff([obj.ROI.nZones+1 z obj.ROI.nZones+1]) ~= 0;
            idx = find(d);
            dd = diff(idx);
            durations = dd(ez(idx(1:end-1)) ~= 0);
            %%
            Profile.Individual(end).Val(s, 1) = nanmean(durations) / fps;
        end
            
        %%
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end

    %% mean time in locations
    Profile.Individual(end+1).Name = 'std time in locations';
    try
        %%
        for s=1:obj.nSubjects
            z = obj.Tracking.zones(s, obj.Tracking.valid);
            map = ismember(z, Locations);
            z = z .* uint16(map);
            %%
            ez = [z, 0];
            d = diff([obj.ROI.nZones+1 z obj.ROI.nZones+1]) ~= 0;
            idx = find(d);
            dd = diff(idx);
            durations = dd(ez(idx(1:end-1)) ~= 0);
            %%
            Profile.Individual(end).Val(s, 1) = nanstd(durations / fps);
        end
            
        %%
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Tangential velocity
    Profile.Individual(end+1).Name = 'Tangential velocity [m/sec]';
    Profile.Individual(end+1).Name = 'Angular velocity [rad/sec]';
    try
        winsz = round(opt.SpeedWinDuration * fps);
        for s=1:obj.nSubjects
            pos = Segs(~obj.Tracking.sheltered(s, 1:length(obj.Tracking.valid)) & obj.Tracking.valid, [obj.Tracking.x(s, :); obj.Tracking.y(s, :)]');
            TV = nan(1, pos.Length);
            AV = nan(1, pos.Length);
            for i=1:length(pos.Events)
                %%
                d = pos.Events(i).data;
                v = d(1:end-winsz, :) - d(1+winsz:end, :);
                if size(v, 1) < 2
                    continue;
                end
                n  = sqrt(sum(v.^2, 2));
                nv = bsxfun(@rdivide, v, n);
                av = acos(max(min(sum(nv(2:end, :) .* nv(1:end-1, :),2), 1), -1)) / (winsz / fps);
                tv = av .* (n(2:end)  / PixelsPerCM / 100);
                jump = n / PixelsPerCM;
                tv(jump(1:end-1)<opt.SpeedWinThreshold | jump(2:end)<opt.SpeedWinThreshold) = nan;
                av(jump(1:end-1)<opt.SpeedWinThreshold | jump(2:end)<opt.SpeedWinThreshold) = nan;
                TV(pos.Events(i).beg:pos.Events(i).beg+length(tv)-1) = tv;
                AV(pos.Events(i).beg:pos.Events(i).beg+length(av)-1) = av;
            end
            Profile.Individual(end-1).Val(s, 1) = nanmean(TV);
            Profile.Individual(end).Val(s, 1) = nanmean(AV);
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
        Profile.Individual(end-1).Val = nan(obj.nSubjects, 1);
    end
    
    
    %% Distance
    Profile.Individual(end+1).Name = 'Distance outside (rough) [m]';
    try
        Profile.Individual(end).Val = zeros(obj.nSubjects, 1);
        if Q.isfield(obj, 'Meta.Scale.PixPerCM')
            PixelsPerCM = obj.Meta.Scale.PixPerCM;
        else
            PixelsPerCM = mean([obj.Meta.Scale.ArenaCoord(3) / obj.Meta.Scale.ArenaWidth, obj.Meta.Scale.ArenaCoord(4) / obj.Meta.Scale.ArenaHeight]);
        end
        jumpTime = 1; %% sec
        for s=1:obj.nSubjects
            %%
            x = Segs(~obj.Tracking.sheltered(s, 1:length(obj.Tracking.valid)) & obj.Tracking.valid, obj.Tracking.x(s, :));
            y = Segs(~obj.Tracking.sheltered(s, 1:length(obj.Tracking.valid)) & obj.Tracking.valid, obj.Tracking.y(s, :));
            dt = round(fps * jumpTime);
            D = nan(1, x.Length);
            for i=1:length(x.Events)
                cx = x.Events(i).data;
                cx = cx(1:dt:end);
                cy = y.Events(i).data;
                cy = cy(1:dt:end);
                if length(cx) > 2
                    d = sqrt((cx(2:end)-cx(1:end-1)).^2 + (cy(2:end)-cy(1:end-1)).^2) / (PixelsPerCM);
                    d(d < opt.DistanceNoiseThreshold) = 0;
                    D(x.Events(i).beg:x.Events(i).beg+length(d)-1) = d / 100;
                end
            end
            %%
            Profile.Individual(end).Val(s, 1) =  nansum(D);
        end
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Grid Entropy
    griddsize = [6 6];
    Profile.Individual(end+1).Name = sprintf('GridEntropy%dx%d', griddsize(1), griddsize(2));
    try
        
        x = (obj.Tracking.x - obj.Meta.Scale.ArenaCoord(1)) / obj.Meta.Scale.ArenaCoord(3);
        x = min(max(x, 0), 1);
        x = min(floor(x*griddsize(1)) + 1, griddsize(1));
        
        y = (obj.Tracking.y - obj.Meta.Scale.ArenaCoord(2)) / obj.Meta.Scale.ArenaCoord(4);
        y = min(max(y, 0), 1);
        y = min(floor(y*griddsize(2)) + 1, griddsize(2));
        
        ind = reshape(sub2ind(griddsize, x(:), y(:)), size(x));
        p = InfoTheory.IndepProbs(ind', 1:(griddsize(1)*griddsize(2)));
        h = -sum(p .* log2(p));
        Profile.Individual(end).Val = h(:);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Grid MI
    griddsize = [6 6];
    Profile.Individual(end+1).Name = sprintf('GridMI%dx%d', griddsize(1), griddsize(2));
    try
        %%
        x = (obj.Tracking.x - obj.Meta.Scale.ArenaCoord(1)) / obj.Meta.Scale.ArenaCoord(3);
        x = min(max(x, 0), 1);
        x = min(floor(x*griddsize(1)) + 1, griddsize(1));
        
        y = (obj.Tracking.y - obj.Meta.Scale.ArenaCoord(2)) / obj.Meta.Scale.ArenaCoord(4);
        y = min(max(y, 0), 1);
        y = min(floor(y*griddsize(2)) + 1, griddsize(2));
        
        ind = reshape(sub2ind(griddsize, x(:), y(:)), size(x));
        %%
        for s=1:obj.nSubjects
            %%
            uniqs = ind(s, [1 diff(ind(s, :))] ~= 0);
            tr = zeros(prod(griddsize));
            tr(1:prod(griddsize)^2) = histc(sub2ind(size(tr), uniqs(1:end-1)', uniqs(2:end)'), 1:prod(griddsize)^2);
            p = tr / sum(tr(:));
            Profile.Individual(end).Val(s,1) = InfoTheory.MutualInformation(p) / InfoTheory.Entropy(sum(p, 1));
        end
        %%
        
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Grid MultiInfo
    griddsize = [6 6];
    Profile.Individual(end+1).Name = sprintf('GroupGridMultiInfo%dx%d', griddsize(1), griddsize(2));
    try
        
        x = (obj.Tracking.x - obj.Meta.Scale.ArenaCoord(1)) / obj.Meta.Scale.ArenaCoord(3);
        x = min(max(x, 0), 1);
        x = min(floor(x*griddsize(1)) + 1, griddsize(1));
        
        y = (obj.Tracking.y - obj.Meta.Scale.ArenaCoord(2)) / obj.Meta.Scale.ArenaCoord(4);
        y = min(max(y, 0), 1);
        y = min(floor(y*griddsize(2)) + 1, griddsize(2));
        
        ind = reshape(sub2ind(griddsize, x(:), y(:)), size(x));
        mi = InfoTheory.Entropy(InfoTheory.IndepProbsVec(ind')) - InfoTheory.Entropy(InfoTheory.JointProbsVec(ind'));
        Profile.Individual(end).Val = mi * ones(obj.nSubjects, 1);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end

    %% Leaves
    Profile.Individual(end+1).Name = 'Number of Leaves';
    try
        c = Hierarchy.Leaves(obj);
        Profile.Individual(end).Val(:, 1) = sum(c, 2);
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Number of leaves per time outside
    Profile.Individual(end+1).Name = sprintf('Leaves rate outside [1/hour]');
    try
        c = Hierarchy.Leaves(obj);
        p = sum(c, 2);
        Profile.Individual(end).Val = p ./ TimeOutside;
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Number of leaves per contact
    Profile.Individual(end+1).Name = 'Fraction of leaves per contact';
    try
        c = Hierarchy.Leaves(obj);
        nc = Hierarchy.Contacts(obj);
        Profile.Individual(end).Val = sum(c, 2) ./ sum(nc, 2);
        Profile.Individual(end).Val(Missing | NoContacts) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Number of leaves per hour
    Profile.Individual(end+1).Name = 'Leave rate [1/hour]';
    try
        c = Hierarchy.Leaves(obj);
        Profile.Individual(end).Val = sum(c, 2) / VideoDuration;
        Profile.Individual(end).Val(Missing) = 0;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
        Profile.Individual(end).Val = nan(obj.nSubjects, 1);
    end
    
    %% Potts
    if opt.TrainMaxEnt
        Profile.Individual(end+1).Name = 'GroupPotts2nd';
        Profile.Individual(end+1).Name = 'GroupPotts3rd';
        Profile.Individual(end+1).Name = 'GroupPotts4th';
        Profile.Individual(end+1).Name = 'GroupPotts2ndFraction';
        Profile.Individual(end+1).Name = 'GroupPotts3rdFraction';
        Profile.Individual(end+1).Name = 'GroupPotts4thFraction';
        try
            %%
            Console.WriteLine('computing Potts model...');
            Console.Write(1, '2nd order');
            potts2 = MaxEntropyGeneral.Train(double(obj.Tracking.zones(:, obj.Tracking.valid)')-1, 2, obj.ROI.nZones, 'Verbose', false);
            Console.Done();
            Console.Write(1, '3rd order');
            potts3 = MaxEntropyGeneral.Train(double(obj.Tracking.zones(:, obj.Tracking.valid)')-1, 3, obj.ROI.nZones, 'Verbose', false);
            Console.Done();
            %%
            Profile.Individual(end-5).Val = potts2.ConnInfo * ones(obj.nSubjects, 1);
            Profile.Individual(end-4).Val = potts3.ConnInfo * ones(obj.nSubjects, 1);
            Profile.Individual(end-3).Val = potts3.EmpMultiInfo * ones(obj.nSubjects, 1);
            Profile.Individual(end-2).Val = (potts2.ConnInfo ./ potts3.EmpMultiInfo) * ones(obj.nSubjects, 1);
            Profile.Individual(end-1).Val = (potts3.ConnInfo ./ potts3.EmpMultiInfo) * ones(obj.nSubjects, 1);
            Profile.Individual(end-0).Val = ones(obj.nSubjects, 1);
        catch me
            fprintf('\n');
            Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
            Profile.Individual(end-5).Val = nan(obj.nSubjects, 1);
            Profile.Individual(end-4).Val = nan(obj.nSubjects, 1);
            Profile.Individual(end-3).Val = nan(obj.nSubjects, 1);
            Profile.Individual(end-2).Val = nan(obj.nSubjects, 1);
            Profile.Individual(end-1).Val = nan(obj.nSubjects, 1);
            Profile.Individual(end-0).Val = nan(obj.nSubjects, 1);
        end
        %% Foraging MultiInfo
        Profile.Individual(end+1).Name = sprintf('GroupForagingMultiInfo');
        try
            
            foraging = double(obj.Tracking.sheltered(:, obj.Tracking.valid)');
            ForagingMultiInfo = InfoTheory.Entropy(InfoTheory.IndepProbsVec(foraging)) - InfoTheory.Entropy(InfoTheory.JointProbsVec(foraging));
            Profile.Individual(end).Val = ForagingMultiInfo * ones(obj.nSubjects, 1);
        catch me
            Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
            Profile.Individual(end).Val = nan(obj.nSubjects, 1);
        end
        
        %% Foraging
        Profile.Individual(end+1).Name = 'GroupForage2nd';
        Profile.Individual(end+1).Name = 'GroupForage3rd';
        Profile.Individual(end+1).Name = 'GroupForage4th';
        Profile.Individual(end+1).Name = 'GroupForage2ndFraction';
        Profile.Individual(end+1).Name = 'GroupForage3rdFraction';
        Profile.Individual(end+1).Name = 'GroupForage4thFraction';
        try
            %%
            data = double(~obj.Tracking.sheltered(:, obj.Tracking.valid))';
            Console.WriteLine('computing foraging model...');
            Console.Write(1, '2nd order');
            forage2 = MaxEntropyGeneral.Train(data, 2, 2, 'Verbose', false);
            Console.Done();
            Console.Write(1, '3rd order');
            forage3 = MaxEntropyGeneral.Train(data, 3, 2, 'Verbose', false);
            Console.Done();
            Profile.Individual(end-5).Val = forage2.ConnInfo * ones(obj.nSubjects, 1);
            Profile.Individual(end-4).Val = forage3.ConnInfo * ones(obj.nSubjects, 1);
            Profile.Individual(end-3).Val = forage3.EmpMultiInfo * ones(obj.nSubjects, 1);
            Profile.Individual(end-2).Val = (forage2.ConnInfo ./ forage3.EmpMultiInfo) * ones(obj.nSubjects, 1);
            Profile.Individual(end-1).Val = (forage3.ConnInfo ./ forage3.EmpMultiInfo) * ones(obj.nSubjects, 1);
            Profile.Individual(end-0).Val = ones(obj.nSubjects, 1);
        catch me
            fprintf('\n');
            Console.Warning(me, ['# failed for ''' Profile.Individual(end).Name ''' - %s\n']);
            Profile.Individual(end-5).Val = nan(obj.nSubjects, 1);
            Profile.Individual(end-4).Val = nan(obj.nSubjects, 1);
            Profile.Individual(end-3).Val = nan(obj.nSubjects, 1);
            Profile.Individual(end-2).Val = nan(obj.nSubjects, 1);
            Profile.Individual(end-1).Val = nan(obj.nSubjects, 1);
            Profile.Individual(end-0).Val = nan(obj.nSubjects, 1);
        end
    end
end

%% Segmented Behaviors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Individual behaviors in different segments of the day (usually 6 segs of
% 2 hours each)
for SegmentsSection = 1
    nsegs = 6;
    win = fps * 60 * 60 * 2;
    
    %% Contacts
    Profile.Segments(1).Name = 'Contacts';
    try
        for i=1:nsegs
            c = Hierarchy.Contacts(obj, round([(i-1) * win + 1, i*win]));
            Profile.Segments(end).Val(:, i) = sum(c, 1)';
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Segments(end).Name ''' - %s\n']);
        Profile.Segments(end).Val = nan(obj.nSubjects, nsegs);
    end
    
    %% Mean contact duration
    Profile.Segments(1).Name = 'Mean contact duration';
    try
        %%
        for i=1:nsegs
            [c, d] = Hierarchy.Contacts(obj, round([(i-1) * win + 1, i*win]));
            Profile.Segments(end).Val(:, i) = sum(d / fps, 2) ./ sum(c, 2);
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Segments(end).Name ''' - %s\n']);
        Profile.Segments(end).Val = nan(obj.nSubjects, nsegs);
    end
    
    %% Aggressive chase
    Profile.Segments(end+1).Name = 'Aggressive chases';
    try
        for i=1:nsegs
            c = Hierarchy.ChaseEscape(obj, round([(i-1) * win + 1, i*win]));
            Profile.Segments(end).Val(:, i) = sum(c, 2);
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Segments(end).Name ''' - %s\n']);
        Profile.Segments(end).Val = nan(obj.nSubjects, nsegs);
    end
    %% Aggressive escapes
    Profile.Segments(end+1).Name = 'Aggressive escapes';
    try
        for i=1:nsegs
            c = Hierarchy.ChaseEscape(obj, round([(i-1) * win + 1, i*win]));
            Profile.Segments(end).Val(:, i) = sum(c, 1)';
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Segments(end).Name ''' - %s\n']);
        Profile.Segments(end).Val = nan(obj.nSubjects, nsegs);
    end
    
    %% Approaches
    Profile.Segments(end+1).Name = 'Approaches';
    try
        for i=1:nsegs
            c = Hierarchy.Approaches(obj, round([(i-1) * win + 1, i*win]));
            Profile.Segments(end).Val(:, i) = sum(c, 1)';
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Segments(end).Name ''' - %s\n']);
        Profile.Segments(end).Val = nan(obj.nSubjects, nsegs);
    end
    
    %% Being appraoched
    Profile.Segments(end+1).Name = 'Being approached';
    try
        for i=1:nsegs
            c = Hierarchy.Approaches(obj, round([(i-1) * win + 1, i*win]));
            Profile.Segments(end).Val(:, i) = sum(c, 2);
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Segments(end).Name ''' - %s\n']);
        Profile.Segments(end).Val = nan(obj.nSubjects, nsegs);
    end
    
    
    %% Time outside
    Profile.Segments(end+1).Name = 'Percentage of time outside';
    try
        %%
        for i=1:nsegs
            r = round((i-1) * win + 1):round(i*win);
            r(r > length(obj.Tracking.valid)) = [];
            c = sum(~obj.Tracking.sheltered(:, r), 2) / sum(obj.Tracking.valid(r), 2);
            Profile.Segments(end).Val(:, i) = c;
        end
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Segments(end).Name ''' - %s\n']);
        Profile.Segments(end).Val = nan(obj.nSubjects, nsegs);
    end
    
    %% Entropy
    Profile.Segments(end+1).Name = 'Entropy';
    try
        %%
        for i=1:nsegs
            H = zeros(obj.nSubjects, 1);
            r = round((i-1) * win + 1):round(i*win);
            r(r > length(obj.Tracking.valid)) = [];
            for s=1:obj.nSubjects
                z = obj.Tracking.zones(s, r);
                z(~obj.Tracking.valid(r)) = [];
                p = histc(z, 1:obj.ROI.nZones);
                p = p / sum(p);
                H(s) = InfoTheory.Entropy(p);
            end
            Profile.Segments(end).Val(:, i) = H;
        end
        
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Segments(end).Name ''' - %s\n']);
        Profile.Segments(end).Val = nan(obj.nSubjects, nsegs);
    end
    
    %% High place
    Profile.Segments(end+1).Name = 'High place';
    try
        %%
        for i=1:nsegs
            r = round((i-1) * win + 1):round(i*win);
            r(r > length(obj.Tracking.valid)) = [];
            r(~obj.Tracking.valid(r)) = [];
            v = zeros(obj.nSubjects, 1);
            for s=1:obj.nSubjects
                z = obj.Tracking.zones(s, r);
                high = sum(ismember(z, HighPlaces));
                total = size(z, 2);
                v(s) = high / total;
            end
            Profile.Segments(end).Val(:, i) = v;
        end
        
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Segments(end).Name ''' - %s\n']);
        Profile.Segments(end).Val = nan(obj.nSubjects, nsegs);
    end
    
    %% Food or water
    Profile.Segments(end+1).Name = 'Food or water';
    try
        %%
        for i=1:nsegs
            r = round((i-1) * win + 1):round(i*win);
            r(r > length(obj.Tracking.valid)) = [];
            r(~obj.Tracking.valid(r)) = [];
            v = zeros(obj.nSubjects, 1);
            for s=1:obj.nSubjects
                z = obj.Tracking.zones(s, r);
                high = sum(ismember(z, FoodWater));
                total = size(z, 2);
                v(s) = high / total;
            end
            Profile.Segments(end).Val(:, i) = v;
        end
        
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Segments(end).Name ''' - %s\n']);
        Profile.Segments(end).Val = nan(obj.nSubjects, nsegs);
    end
    
    %% Distance from wall
    Profile.Segments(end+1).Name = 'Distance from walls';
    try
        %%
        ac = obj.Meta.Scale.ArenaCoord;
        wall_n = ac(1):1:ac(1)+ac(3); wall_n(2, :) = ac(2);
        wall_s = ac(1):1:ac(1)+ac(3); wall_s(2, :) = ac(2)+ac(4);
        wall_e(2, :) = ac(2):1:ac(2)+ac(4); wall_e(1, :) = ac(1);
        wall_w(2, :) = ac(2):1:ac(2)+ac(4); wall_w(1, :) = ac(1)+ac(3);
        wall = [wall_n, wall_s, wall_e, wall_w]';
        
        for i=1:nsegs
            r = round((i-1) * win + 1):round(i*win);
            r(r > length(obj.Tracking.valid)) = [];
            r(~obj.Tracking.valid(r)) = [];
            v = zeros(obj.nSubjects, 1);
            for s=1:obj.nSubjects
                %%
                R = r;
                R(obj.Tracking.zones(s, r) ~= find(ismember(obj.ROI.ZoneNames, {'Open'}))) = [];
                x = obj.Tracking.x(s, R);
                y = obj.Tracking.y(s, R);
                v(s) = mean(sqrt(sum((wall(knnsearch(wall, [x(:), y(:)]), :) - [x(:), y(:)]).^2, 2))) / PixelsPerCM;
            end
            Profile.Segments(end).Val(:, i) = v;
        end
        
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Segments(end).Name ''' - %s\n']);
        Profile.Segments(end).Val = nan(obj.nSubjects, nsegs);
    end
    
    %% Distance from nest
    Profile.Segments(end+1).Name = 'Distance from nest';
    try
        %%
        nests = cumsum(cellfun(@(x) ~isempty(x), regexp(obj.ROI.ZoneNames, '\(.*\)')));
        nestid = nests(strcmp(obj.ROI.ZoneNames, '(BigNest)'));
        bounds = bwboundaries(obj.ROI.Hidden{nestid});
        bounds = bounds{1};
        
        for i=1:nsegs
            r = round((i-1) * win + 1):round(i*win);
            r(r > length(obj.Tracking.valid)) = [];
            r(~obj.Tracking.valid(r)) = [];
            v = zeros(obj.nSubjects, 1);
            for s=1:obj.nSubjects
                R = r;
                R(obj.Tracking.sheltered(s, r)) = [];
                
                %%
                x = obj.Tracking.x(s, R);
                y = obj.Tracking.y(s, R);
                v(s) = mean(sqrt(sum((bounds(knnsearch(bounds, [x(:), y(:)]), :) - [x(:), y(:)]).^2, 2))) / PixelsPerCM;
            end
            Profile.Segments(end).Val(:, i) = v;
        end
        
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Segments(end).Name ''' - %s\n']);
        Profile.Segments(end).Val = nan(obj.nSubjects, nsegs);
    end
    
    %% Speed
    Profile.Segments(end+1).Name = 'Speed';
    try
        %%
        jumpTime = opt.DistanceNoiseDuration; %% sec
        
        for i=1:nsegs
            r = round((i-1) * win + 1):round(i*win);
            r(r > length(obj.Tracking.valid)) = [];
            v = zeros(obj.nSubjects, 1);
            for s=1:obj.nSubjects
                x = Segs(~obj.Tracking.sheltered(s, r) & obj.Tracking.valid(r), obj.Tracking.x(s, r));
                y = Segs(~obj.Tracking.sheltered(s, r) & obj.Tracking.valid(r), obj.Tracking.y(s, r));
                
                dt = round(fps * jumpTime);
                speed = nan(1, x.Length);
                for e=1:length(x.Events)
                    cx = x.Events(e).data;
                    cx = cx(1:dt:end);
                    cy = y.Events(e).data;
                    cy = cy(1:dt:end);
                    if length(cx) > 2
                        d = sqrt((cx(2:end)-cx(1:end-1)).^2 + (cy(2:end)-cy(1:end-1)).^2) / PixelsPerCM;
                        d(d < opt.DistanceNoiseThreshold) = 0;
                        speed(x.Events(e).beg:x.Events(e).beg+length(d)-1) = d / 100 / jumpTime;
                    end
                end
                v(s) = nanmean(speed);
            end
            Profile.Segments(end).Val(:, i) = v;
        end
        
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Segments(end).Name ''' - %s\n']);
        Profile.Segments(end).Val = nan(obj.nSubjects, nsegs);
    end
    
    %% Distance
    Profile.Segments(end+1).Name = 'Distance';
    try
        %%
        jumpTime = opt.DistanceNoiseDuration; %% sec
        
        for i=1:nsegs
            r = round((i-1) * win + 1):round(i*win);
            r(r > length(obj.Tracking.valid)) = [];
            v = zeros(obj.nSubjects, 1);
            for s=1:obj.nSubjects
                x = Segs(~obj.Tracking.sheltered(s, r) & obj.Tracking.valid(r), obj.Tracking.x(s, r));
                y = Segs(~obj.Tracking.sheltered(s, r) & obj.Tracking.valid(r), obj.Tracking.y(s, r));
                
                dt = round(fps * jumpTime);
                
                D = nan(1, x.Length);
                for e=1:length(x.Events)
                    cx = x.Events(e).data;
                    cx = cx(1:dt:end);
                    cy = y.Events(e).data;
                    cy = cy(1:dt:end);
                    if length(cx) > 2
                        d = sqrt((cx(2:end)-cx(1:end-1)).^2 + (cy(2:end)-cy(1:end-1)).^2) / PixelsPerCM;
                        d(d < opt.DistanceNoiseThreshold) = 0;
                        D(x.Events(e).beg:x.Events(e).beg+length(d)-1) = d / 100;
                    end
                end
                
                v(s) = nansum(D);
            end
            Profile.Segments(end).Val(:, i) = v;
        end
        
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Segments(end).Name ''' - %s\n']);
        Profile.Segments(end).Val = nan(obj.nSubjects, nsegs);
    end
    
end

%% Pairwise Behaviors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Organized in matrices where row represents the 'from-mouse' and column
% the 'to-mouse'. For example, for chases the rows represent the chaser and
% the columns the "escaper".
for PairwiseSection = 1
    %% Aggressive chase-escape
    Profile.Pairwise(1).Name = 'Number of aggressive chase-escapes';
    try
        Profile.Pairwise(end).Val = obj.Hierarchy.AggressiveChase.ChaseEscape;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Pairwise(end).Name ''' - %s\n']);
        Profile.Pairwise(end).Val = nan(obj.nSubjects, obj.nSubjects);
    end
    
    %% Chase-escape
    Profile.Pairwise(end + 1).Name = 'Number of chase-escapes';
    try
        Profile.Pairwise(end).Val = obj.Hierarchy.ChaseEscape.ChaseEscape;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Pairwise(end).Name ''' - %s\n']);
        Profile.Pairwise(end).Val = nan(obj.nSubjects, obj.nSubjects);
    end
    
    %% Contacts
    Profile.Pairwise(end + 1).Name = 'Number of contacts';
    try
        c = zeros(obj.nSubjects);
        for l=1:length(ContactList)
            curr = ContactList(l);
            s1 = curr.subjects(1);
            s2 = curr.subjects(2);
            c(s1, s2) = c(s1, s2) + 1;
            c(s2, s1) = c(s2, s1) + 1;
        end
        Profile.Pairwise(end).Val = c;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Pairwise(end).Name ''' - %s\n']);
        Profile.Pairwise(end).Val = nan(obj.nSubjects, obj.nSubjects);
    end
    
    %% Approaches
    Profile.Pairwise(end + 1).Name = 'Approaches';
    try
        c = zeros(obj.nSubjects);
        for l=1:length(ContactList)
            curr = ContactList(l);
            s1 = curr.subjects(1);
            s2 = curr.subjects(2);
            c(s1, s2) = c(s1, s2) + curr.states(1, 2);
            c(s2, s1) = c(s2, s1) + curr.states(2, 2);
        end
        Profile.Pairwise(end).Val = c;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Pairwise(end).Name ''' - %s\n']);
        Profile.Pairwise(end).Val = nan(obj.nSubjects, obj.nSubjects);
    end
    
    %% Approaches per hour
    Profile.Pairwise(end + 1).Name = 'Approaches per hour';
    try
        c = zeros(obj.nSubjects);
        for l=1:length(ContactList)
            curr = ContactList(l);
            s1 = curr.subjects(1);
            s2 = curr.subjects(2);
            c(s1, s2) = c(s1, s2) + curr.states(1, 2);
            c(s2, s1) = c(s2, s1) + curr.states(2, 2);
        end
        Profile.Pairwise(end).Val = c  / VideoDuration;
    catch me
        Console.Warning(me, ['# failed for ''' Profile.Pairwise(end).Name ''' - %s\n']);
        Profile.Pairwise(end).Val = nan(obj.nSubjects, obj.nSubjects);
    end
end
%% Behavioral Events %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Events: Chase-escape
Profile.Events(1).Name = 'Chase-escape [chaser escaper start_frame end_frame]';
try
    %%
    map = obj.Hierarchy.Contacts.Behaviors.AggressiveChase.Map;
    b = min([ContactList(map).beg]);
    e = max([ContactList(map).end]);
    s = [ContactList(map).subjects];
    aggr = cat(3, ContactList(map).states);
    
    id = squeeze([aggr(:, 5, :)]);
    [~, ids] = max(id);
    chaser = s(sub2ind(size(s), ids, 1:size(s, 2)));
    escaper = s(sub2ind(size(s), 3-ids, 1:size(s, 2)));
    
    Profile.Events(end).Val = [chaser(:) escaper(:) b(:) e(:)];
catch me
    Console.Warning(me, ['# failed for ''' Profile.Events(end).Name ''' - %s\n']);
    Profile.Events(end).Val = nan;
end

%% Create table of behaviors for each individual
% Using data from Profile.Individual
names = cell(1, length(Profile.Individual));
units = cell(1, length(Profile.Individual));
descs = cell(1, length(Profile.Individual));
rows = cell(1, obj.nSubjects);
for i=1:obj.nSubjects
    rows{i} = sprintf('%s.exp%04d.sub%d', obj.Meta.GroupType, obj.Meta.GroupId, i);
end
t = table('rownames', rows);
for i=1:length(Profile.Individual)
    %%
    name = Profile.Individual(i).Name;
    names{i} = regexprep(regexprep(name, '\[.*\]', ''), '\(.*\)', '');
    desc = regexp(name, '\((?<desc>.*)\)','match');
    if isempty(desc); desc = {''}; end
    unit = regexp(name, '\[(?<unit>.*)\]','match');
    if isempty(unit); unit = {''}; end
    
    units{i} = [unit{:}];
    descs{i} = [desc{:}];
    
    vname = regexprep(Q.capitalize(names{i}), ' ', '');
    t = [t, table(Profile.Individual(i).Val, 'VariableNames', {vname})];
    t.Properties.VariableDescriptions{vname} = [desc{:}];
end
t.Properties.VariableUnits = units;
n = CheeseSquare.ParseName(obj.Prefix);
gtype = cell(size(t, 1), 1);
for i=1:size(t, 1)
    gtype{i} = n.GroupType;
end
metatable = table(n.DayId * ones(size(t, 1), 1), n.GroupId * ones(size(t, 1), 1), 1 * ones(size(t, 1), 1), gtype, Q.tocol(1:size(t, 1)), Q.tocol(1:size(t, 1)), 'VariableNames', {'Day', 'GroupID', 'GroupNumber', 'GroupType', 'MouseID', 'MouseNumber'});
metatable.Properties.VariableUnits = {'meta', 'meta', 'meta', 'meta', 'meta', 'meta'};
t = [metatable, t];
t.Properties.UserData = struct('IsMeta', strcmp(t.Properties.VariableUnits, 'meta'));

Profile.Tables.Individual = t;

%% Add pairwise interactions to table 
% The data is stored as a vector inside the table
pairtable =t(:, []);
try
    nApp = zeros(obj.nSubjects);
    nCon = zeros(obj.nSubjects);
    nChas = obj.Hierarchy.AggressiveChase.ChaseEscape;
    for l=1:length(ContactList)
        curr = ContactList(l);
        s1 = curr.subjects(1);
        s2 = curr.subjects(2);
        nApp(s1, s2) = nApp(s1, s2) + curr.states(1, 2);
        nApp(s2, s1) = nApp(s2, s1) + curr.states(2, 2);
        nCon(s1, s2) = nCon(s1, s2) + 1;
        nCon(s2, s1) = nCon(s2, s1) + 1;
        
    end
    
    % contact rate
    curr = array2table(zeros(size(nCon, 1), 1), 'VariableNames', {'PairwiseContactRate'});
    curr.PairwiseContactRate = nCon / VideoDuration;
    pairtable = [pairtable, curr];
    % appraoch rate
    curr = array2table(zeros(size(nApp, 1), 1), 'VariableNames', {'PairwiseApproachRate'});
    curr.PairwiseApproachRate = nApp / VideoDuration;
    pairtable = [pairtable, curr];
    % being appraoched rate
    curr = array2table(zeros(size(nApp, 2), 1), 'VariableNames', {'PairwiseBeingApproachedRate'});
    curr.PairwiseBeingApproachedRate = nApp' / VideoDuration;
    pairtable = [pairtable, curr];
    % chase rate
    curr = array2table(zeros(size(nChas, 1), 1), 'VariableNames', {'PairwiseChaseRate'});
    curr.PairwiseChaseRate = nChas / VideoDuration;
    pairtable = [pairtable, curr];
    % escape rate
    curr = array2table(zeros(size(nChas, 1), 1), 'VariableNames', {'PairwiseEscapeRate'});
    curr.PairwiseEscapeRate = nChas' / VideoDuration;
    pairtable = [pairtable, curr];
    % contact count
    curr = array2table(zeros(size(nCon, 1), 1), 'VariableNames', {'PairwiseContactCount'});
    curr.PairwiseContactCount = nCon;
    pairtable = [pairtable, curr];
    % appraoch count
    curr = array2table(zeros(size(nApp, 1), 1), 'VariableNames', {'PairwiseApproachCount'});
    curr.PairwiseApproachCount = nApp;
    pairtable = [pairtable, curr];
    % being appraoched count
    curr = array2table(zeros(size(nApp, 2), 1), 'VariableNames', {'PairwiseBeingApproachedCount'});
    curr.PairwiseBeingApproachedCount = nApp';
    pairtable = [pairtable, curr];
    % chase count
    curr = array2table(zeros(size(nChas, 1), 1), 'VariableNames', {'PairwiseChaseCount'});
    curr.PairwiseChaseCount = nChas;
    pairtable = [pairtable, curr];
    % escape count
    curr = array2table(zeros(size(nChas, 1), 1), 'VariableNames', {'PairwiseEscapeCount'});
    curr.PairwiseEscapeCount = nChas';
    pairtable = [pairtable, curr];
    
    pOut = zeros(obj.nSubjects);
    out = ~obj.Tracking.sheltered(:, obj.Tracking.valid);
    for i=1:obj.nSubjects
        for j=i+1:obj.nSubjects
            pOut(i,j) = mean(out(i, :) & out(j, :));
            pOut(j,i) = pOut(i,j);
        end
    end
    
    % fraction time outside together
    curr = array2table(zeros(size(pOut, 1), 1), 'VariableNames', {'PairwiseFractionTimeOutside'});
    curr.PairwiseFractionTimeOutside = pOut;
    pairtable = [pairtable, curr];
    
    % foraging correlation
    curr = array2table(zeros(size(nChas, 1), 1), 'VariableNames', {'PairwiseForagingCorr'});
    curr.PairwiseForagingCorr = corr(out');
    pairtable = [pairtable, curr];
    
    full = [t, pairtable];
    full.Properties.UserData.IsMeta(end+1:size(full, 2)) = false;
    full.Properties.UserData.IsPairwise = [false(1, size(t, 2)) true(1, size(pairtable, 2))];
    Profile.Tables.Individual = full;
    
catch me
    Console.Warning(me.message);
end

%% Add behaviors computed on time segments of the video (6 segs) to the table
try
    %%
    t = full;
    for i=1:length(Profile.Segments)
        name = Profile.Segments(i).Name;
        wordstart = regexp([' ' name],'(?<=\s+)\S','start')-1;
        name(wordstart) = upper(name(wordstart));
        name = ['Segs' regexprep(name, ' ', '')];
        t = [t, table(Profile.Segments(i).Val, 'VariableNames', {name})];
    end
    t.Properties.UserData.IsMeta(size(full, 2)+1:size(t, 2)) = false;
    t.Properties.UserData.IsPairwise(size(full, 2)+1:size(t, 2)) = false;
    t.Properties.UserData.IsSegs(size(full, 2)+1:size(t, 2)) = true;
    full = t;
    Profile.Tables.Individual = full;
catch me
    Console.Warning(me.message);
end

%% Create Group behaviors table
foraging = double(obj.Tracking.sheltered(:, obj.Tracking.valid)');
ForagingMultiInfo = InfoTheory.Entropy(InfoTheory.IndepProbsVec(foraging)) - InfoTheory.Entropy(InfoTheory.JointProbsVec(foraging));

Profile.Tables.Group = table(...
    n.DayId,...
    n.GroupId, ...
    1,...
    {n.GroupType}, ...
    ForagingMultiInfo, ...
    'VariableNames', {'Day', 'GroupID', 'GroupNumber', 'GroupType', 'ForagingMultiInfo'}, ...
    'RowNames', {sprintf('%s.exp%04d.day%02d', n.GroupType, n.GroupId, n.DayId)});
%% Add meta data to Profile struct
Profile.HighPlaces = HighPlaces;
Profile.FoodWater = FoodWater;

%% Finalize
obj.Profile = Profile;
