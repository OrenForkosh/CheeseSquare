function obj = CheeseTrackPath(obj, nruns, sourcepath)
% CheeseTrackRun Compute the mice trajectories
%
%   Uses segmentation data files created by CheeseColorSegment:
%       (Prefix).segm.(segid).mat, where (prefix) is the video prefix
%       (SC.exp0001.day01.cam02) and (segid) is the part of the video file
%       that was segmented
%
%   obj = CheeseTrackPath(obj, nruns, sourcepath) Tracks the mice in 'obj'.
%   The algorithm loads segmenation data from 'nruns' parts of the video
%   file. 'sourcepath' specifies where the segmentation files are located.
%
%       Created by OREN FORKOSH
%

if nargin < 3
    sourcepath = obj.OutputPath;
end
%% Default number of segmentation files is 100
if ~exist('nruns', 'var') || isempty(nruns)
    nruns = 100;
end
%%
fprintf('# finding paths\n');
%%
obj.NumOfProbBins = 100; % number of histogram bins for mouse movement statistics
obj.MinJumpProb = 0.0; % deprecated! ignore...
obj.AllowHiddenTransitions = false; % allow transition between hidden states
obj.MinObserved = 4; % ignore if observed less than this
obj.MaxHidden = 100; % if hidden more than this amount of frames, assume in shelter
obj.JumpToShelterDuration = 5 * obj.FrameRate; % deprecated! ignore...
obj.MinEventDuration = 4; % ignore events that are shorter than this
obj.MinShelteredDuration = 6; % ignore if time in shelter is lass than this
obj.MinDistanceToShelterCM = 4; % if distance from shelter is less that this, assume in shelter
obj.MaxIgnore = 12; % deprecated! ignore...

%% Find centers of sheltered areas
HiddenSCenters = zeros(obj.ROI.nHidden, 2);
for r=1:obj.ROI.nHidden
    rp = regionprops(obj.ROI.Hidden{r}, 'Centroid');
    HiddenSCenters(r, :) = rp.Centroid * obj.VideoScale;
end

%% loading segmentation data
% reads the information about detected blobs into the 'cents' structure
% for N frames, S subjects, and at most M blobs 
% the main fields are: 
%   x, y - coordinates (NxM)
%   logprob - probability of the blob to belong to each subject (NxMxS)
fprintf('# - loading segmentations\n');
cents = struct();
cents.x       = [];
idx = 1;
nchars = 0;
for i=1:nruns
    nchars = Console.Reprintf(nchars, '#   . segment no. %3d/%d', i, nruns);
    filename = [sourcepath obj.FilePrefix '.segm.' sprintf('%03d', i) '.mat'];
    currSegm = load(filename);
    currSegm.cents = rmfield(currSegm.cents, 'area');
    currSegm.cents = rmfield(currSegm.cents, 'solidity');
    if size(currSegm.cents.logprob, 3) < obj.nSubjects + 1
        currSegm.cents.logprob = cat(3, currSegm.cents.logprob, zeros(size(currSegm.cents.logprob, 1), size(currSegm.cents.logprob, 2)));
    end
    f = fields(currSegm.cents);
    if i == 1
        for j=1:length(f)
            sz = size(currSegm.cents.(f{j}));
            cents.(f{j}) = zeros([sz(1) * nruns, sz(2:end)]);
        end
    end
    
    for j=1:length(f)
        try
            if ismatrix(currSegm.cents.(f{j}))
                cents.(f{j})(idx:idx+size(currSegm.cents.(f{j}), 1) - 1, :) = currSegm.cents.(f{j})(:, 1:size(cents.(f{j}), 2));
            elseif ndims(currSegm.cents.(f{j})) == 3
                cents.(f{j})(idx:idx+size(currSegm.cents.(f{j}), 1) - 1, :, :) = currSegm.cents.(f{j})(:, 1:size(cents.(f{j}), 2), 1:size(cents.(f{j}), 3));
            end
        catch
        end
    end
    idx = idx + size(currSegm.cents.x, 1);
end
Console.Reprintf(nchars, '#   . segment no. %3d/%d\n', i, nruns);

for j=1:length(f)
    try
    if ismatrix(cents.(f{j}))
        cents.(f{j}) = cents.(f{j})(1:min(idx-1, size(cents.(f{j}), 1)), :);
    else
        cents.(f{j}) = cents.(f{j})(1:idx-1, :, :);
    end
    catch
    end
end
cents.x = single(cents.x);
cents.y = single(cents.y);
cents.logprob = log(bsxfun(@rdivide, exp(cents.logprob), sum(exp(cents.logprob), 3))); % normalize
cents.logprob = single(cents.logprob);

%% 
cents.cumlabel = zeros(obj.nSubjects, size(cents.label, 1));
for i=1:obj.nSubjects
    [b,e,l] = FindEvents(any(cents.label == i, 2)');
    for j=1:length(b)
        cents.cumlabel(i, b(j):e(j)) = 1:l(j);
    end
end
%% Find boundaries of sheltered regions
for i=1:obj.ROI.nHidden
    [x, y] = find((imdilate(obj.ROI.Hidden{i} ~= 0, ones(3,3)) - (obj.ROI.Hidden{i} ~= 0)) ~= 0);
    obj.ROI.HiddenBoundaryCoordinates{i} = [y, x];
    if Q.isfield(obj, 'ROI.Scaling') && ~isempty(obj.ROI.Scaling)
        obj.ROI.HiddenBoundarySCoordinates{i} = obj.ROI.HiddenBoundaryCoordinates{i} * obj.ROI.Scaling;
    else
        obj.ROI.HiddenBoundarySCoordinates{i} = obj.ROI.HiddenBoundaryCoordinates{i} * obj.VideoScale;
    end
    
end

%% Ignore ROIs
if isfield(obj, 'ROI.Ignore')
    ignore = uint8(imread(obj.ROI.Ignore));
    fprintf('#  . removing unwanted cents...');
    for i=1:size(cents.label, 1)
        labels = cents.label(i, :);
        y = cents.y(i,:);
        x = cents.x(i,:);
        y = y(labels~=0) / obj.VideoScale;
        if isempty(y)
            continue;
        end
        x = x(labels~=0) / obj.VideoScale;
        labels(labels~=0) = (1 - ignore( sub2ind(size(ignore), round(y), round(x)) )) .* labels(labels~=0);
        cents.label(i, :) = labels;
    end
    fprintf('[done]\n');
end

%%
StartFrame = obj.StartTime / obj.dt;
EndFrame = obj.EndTime / obj.dt;
if EndFrame < 0
   EndFrame = obj.nFrames + 1;
end

%%
prev = cents.label;
[q, cents.label] = max(cents.logprob, [], 3);
cents.label(prev == 0) = 0;

thresh = obj.ColorMatchThresh; 
ndata = size(cents.x, 1);

%% Initialize trajectories of each mouse
for i=1:obj.nSubjects
    logprobs = cents.logprob(:, :, i);
    logprobs(cents.label ~= i) = -inf;
    [maxlogprobs, centids] = max(logprobs, [], 2);

    maxprobs = exp(maxlogprobs);
    res{i}.x = cents.x(sub2ind(size(cents.x), 1:ndata, centids'));
    res{i}.y = cents.y(sub2ind(size(cents.y), 1:ndata, centids'));
    res{i}.prob = maxprobs';
    res{i}.logprob = maxlogprobs';
    res{i}.centids = centids';
    res{i}.prev = double([true; cents.prev(sub2ind(size(cents.prev), 2:ndata, centids(2:end)'))' == centids(1:end-1)]);
    
    res{i}.valid = maxprobs' > thresh;
    
    %% remove isolated apperances
    d = diff([0 res{i}.valid 0], 1, 2);
    start  = find(d > 0);
    finish = find(d < 0) - 1;
    len = finish - start + 1;
    isolated = start - [1 finish(1:end-1)] > obj.MaxIsolatedDistance & ...
    [start(2:end) size(res{i}.valid, 1)] - finish > obj.MaxIsolatedDistance;
    idx = find(len == 1 & isolated);
    for j=idx
        res{i}.valid(start(j):finish(j)) = false;
    end
    %%
    res{i}.valid(1:StartFrame) = false;
    res{i}.valid(EndFrame:end) = false;
    %%
    res{i}.id = (maxprobs' > thresh) * i;
    res{i}.x(~res{i}.valid) = nan;
    res{i}.y(~res{i}.valid) = nan;
    res{i}.prob(~res{i}.valid) = nan;
    res{i}.logprob(~res{i}.valid) = -inf;
    res{i}.centids(~res{i}.valid) = nan;
    res{i}.prev(~res{i}.valid) = nan;
    %%
end
%% Compute statistics of mouse movements (between frames)
distances_ = [];
dx_ = [];
dy_ = [];
% Compute distances for matching blobs
for i=1:obj.nSubjects
    m = res{i}.prev(2:end) == res{i}.id(1:end-1)';
    dx = diff(res{i}.x);
    dy = diff(res{i}.y);
    dx_ = [dx_, dx(m)];
    dy_ = [dy_, dy(m)];
    distances = sqrt(diff(res{i}.x).^2 + diff(res{i}.y).^2);
    distances = distances(m);
    distances = distances(distances > 0 & ~isnan(distances));
    distances_ = [distances_, distances];
end

% Compute distances between different colors
distances = [];
idistances = [];
for i=1:obj.nSubjects
    for j=1:obj.nSubjects
        if j~=i
            distances = [distances, sqrt((res{i}.x - res{j}.x).^2 +(res{i}.y - res{j}.y).^2)];
        else
            idistances = [idistances, sqrt((res{i}.x(2:end) - res{j}.x(1:end-1)).^2 +(res{i}.y(2:end) - res{j}.y(1:end-1)).^2)];
        end
    end
end

if obj.nSubjects == 1
    x = res{1}.x(~isnan(res{1}.x));
    y = res{1}.y(~isnan(res{1}.y));
    i = randi(length(x), [1, size(res{1}.x, 2)]);
    x = x(i);
    y = y(i);
    distances = sqrt((x(2:end) - x(1:end-1)).^2 +(y(2:end) - y(1:end-1)).^2);
end

for i=1:obj.nSubjects
    res{i}.probBins = linspace(0, max([distances idistances]) + 0.1, obj.NumOfProbBins+1); 
    res{i}.probBinDiff = res{i}.probBins(2) - res{i}.probBins(1);
    
    h1 = histc(distances, res{i}.probBins) + 1;
    h1=h1/sum(h1); 
    
    h2=histc(distances_, res{i}.probBins);
    h2=h2/sum(h2); 
    
    res{i}.jumpProb = h2./(h1(:)'+h2);
    res{i}.jumpProb(res{i}.jumpProb < obj.MinJumpProb) = 0;
    res{i}.jumpProb(res{i}.jumpProb == 0) = min(res{i}.jumpProb(res{i}.jumpProb > 0)) / 2;
    res{i}.jumpLogprob = log(res{i}.jumpProb);
end

%% Compute trajectory of each mouse
nClusters = size(cents.x, 2);
clusters = 1:nClusters;
nFrames = length(res{1}.id);
scoord = cell(1,obj.nSubjects); 
for curr = 1:obj.nSubjects
    fprintf('# - tracking subject %d of %d\n', curr, obj.nSubjects);
    currLogprobs = cents.logprob(:, :, curr);
    currLogprobs(cents.label == 0) = -inf;
    %% build observation probability table
    table = ones(nClusters+obj.ROI.nHidden+1, nFrames) * -inf;
    observed = zeros(nClusters+obj.ROI.nHidden+1, nFrames, 'int16');
    scoord{curr}.x = zeros(nClusters, nFrames, 'single');
    scoord{curr}.y = zeros(nClusters, nFrames, 'single');
    scoord{curr}.hidingNeighbour.x = scoord{curr}.x * 0 + obj.ROI.HiddenBoundarySCoordinates{1}(1, 1);
    scoord{curr}.hidingNeighbour.y = scoord{curr}.y * 0 + obj.ROI.HiddenBoundarySCoordinates{1}(1, 2);
    scoord{curr}.hidingDistance = ones(nClusters, nFrames, obj.ROI.nHidden) * inf;
    path = zeros(size(table), 'uint16');
    nchar = 0;
    fprintf('#   . computing emission probabilities\n');
    for s=1:nClusters
        range = 1:ndata;
        table(s, range) = currLogprobs(:, s);
        scoord{curr}.x(s, range) = cents.x(:, s); % + inf * (cents.label(:, s) == 0);
        scoord{curr}.y(s, range) = cents.y(:, s); % + inf * (cents.label(:, s) == 0);
        vrange = cents.label(:, s) ~= 0;
        for h=1:obj.ROI.nHidden
            nchar = Console.Reprintf(nchar, '#       processing segment %d / %d', ((s - 1) * obj.ROI.nHidden) + h, nClusters*obj.ROI.nHidden);
            hcoord = single(obj.ROI.HiddenBoundarySCoordinates{h});
            d = pdist2([scoord{curr}.x(s, vrange)', scoord{curr}.y(s, vrange)'], hcoord);
            [hvalues, hidx] = min(d, [], 2);
            scoord{curr}.hidingNeighbour.x(s, vrange) = hcoord(hidx, 1);
            scoord{curr}.hidingNeighbour.y(s, vrange) = hcoord(hidx, 2);
            scoord{curr}.hidingDistance(s, vrange, h) = hvalues';
        end
    end
    fprintf('\n');
    
    for s=1:obj.ROI.nHidden+1
        table(nClusters+s, :) = InfoTheory.Log(prod(1 - exp(table(1:nClusters, :)), 1));
    end
    emitlogprobs = table;
    
    valid = sum(isfinite(table(1:nClusters, :)), 1) > 0;
    d = diff([0 valid 0], 1, 2);
    start  = find(d > 0);
    finish = find(d < 0) - 1;
    
    %%
    prev = 1;
    nchar = 0;
    
    scoord{curr}.hidingPos.x = ones(1, nFrames) * inf;
    scoord{curr}.hidingPos.y = ones(1, nFrames) * inf;
    scoord{curr}.hidingTime = zeros(1, nFrames);
    
    zeroprob = InfoTheory.Log(1);
    zeroprobvec = [repmat(zeroprob, nClusters, 1); zeros(obj.ROI.nHidden+1, 1)];
    
    hiddenprobvec = [repmat(zeroprob, nClusters, 1); zeros(obj.ROI.nHidden, 1); zeroprob];
    nohiddenprobvec = [repmat(zeroprob, nClusters, 1); zeros(obj.ROI.nHidden, 1); InfoTheory.Log(0)];
    
    %% Build Viterbi table
    fprintf('#   . building viterbi table\n');
    for n=1:length(start)
        nchar = Console.Reprintf(nchar, '#       processing segment %d / %d', n, length(start));
        r = start(n):finish(n);
        for f=r
            %%
            for s=1:nClusters
                if ~isfinite(table(s, f))
                    continue;
                end
                x = cents.x(f, s);
                y = cents.y(f, s);
                
                distance = scoord{curr}.hidingDistance(s, f, :);
                jump = sqrt((x - scoord{curr}.x(:, prev)).^2 + (y - scoord{curr}.y(:, prev)).^2);
                if prev+1==f && cents.prev(f, s) > 0
                    jump(cents.prev(f, s)) = 0;
                    jump(clusters(clusters ~= cents.prev(f, s))) = inf;
                end
                distance = [...
                    jump; ...
                    distance(:);...
                    sqrt((x - scoord{curr}.hidingPos.x(prev)).^2 + (y - scoord{curr}.hidingPos.y(prev)).^2)];
%%                    sqrt((x - scoord{curr}.hidingPos.x(prev)).^2 + (y - scoord{curr}.hidingPos.y(prev)).^2) / (f - scoord{curr}.hidingTime(prev))];
                idx = floor(distance / res{curr}.probBinDiff) + 1;
                idx(idx < 1) = 1;
                idx(idx > obj.NumOfProbBins) = obj.NumOfProbBins;
                %idx(isnan(idx)) = obj.NumOfProbBins;
                trans = res{curr}.jumpLogprob(idx)';
                
                if f>prev+1
                    trans = trans + zeroprobvec * (f - prev - 1);
                    if f - prev - 1 > obj.MaxHidden
                        trans(nClusters + obj.ROI.nHidden + 1) = InfoTheory.Log(0);
                        trans(1:nClusters) = InfoTheory.Log(0);
                    end
                end
                
                [val, from] = max(table(:, prev) + trans);
                table(s, f) = table(s, f) + val;
                path(s, f) = from;
                observed(s, f) = observed(from, prev) + 1;
            end
            % hidden zones
            for h=1:obj.ROI.nHidden
                if obj.AllowHiddenTransitions
                    distance = [scoord{curr}.hidingDistance(:, prev, h); obj.ROI.HiddenSDistances(:, h); inf];
                else
                    distance = [scoord{curr}.hidingDistance(:, prev, h); ones(obj.ROI.nHidden,1) * inf; inf];
                    distance(nClusters + h) = 0;
                end
                if obj.MinObserved > 0
                    distance(observed(1:nClusters, f) < obj.MinObserved) = inf;
                end
                idx = floor(distance / res{curr}.probBinDiff) + 1;
                idx(idx < 1) = 1;
                idx(idx > obj.NumOfProbBins) = obj.NumOfProbBins;
                %idx(isnan(idx)) = obj.NumOfProbBins;
                trans = res{curr}.jumpLogprob(idx)';
                
                [val, from] = max(table(:, prev) + trans);
                table(nClusters + h, f) = table(nClusters + h, f) + val;
                path(nClusters + h, f) = from;
            end
            %
            
            if f - scoord{curr}.hidingTime(prev) <= obj.MaxHidden
                [val, from] = max(table(:, prev) + hiddenprobvec);
            else
                [val, from] = max(table(:, prev) + nohiddenprobvec);
            end
            
            table(nClusters + obj.ROI.nHidden + 1, f) = table(nClusters + obj.ROI.nHidden + 1, f) + val;
            path(nClusters + obj.ROI.nHidden + 1, f) = from;
            observed(nClusters + obj.ROI.nHidden + 1, f) = observed(from, prev);

            if from <= nClusters
                scoord{curr}.hidingPos.x(f) = scoord{curr}.x(from, prev);
                scoord{curr}.hidingPos.y(f) = scoord{curr}.y(from, prev);
                scoord{curr}.hidingTime(f) = f;
            else
                scoord{curr}.hidingPos.x(f) = scoord{curr}.hidingPos.x(prev);
                scoord{curr}.hidingPos.y(f) = scoord{curr}.hidingPos.y(prev);
                scoord{curr}.hidingTime(f) = scoord{curr}.hidingTime(prev);
            end
            %
            prev = f;
        end
    end
    fprintf('\n');

    %% Backtrack Viterbi table
    [~, maxidx] = max(table);
    maxpath = repmat([1 maxidx(1:end-1)], size(table, 1), 1);
    path(path == 0) = maxpath(path == 0);
    
    %
    fprintf('#   . backtracking\n');
    track{curr}.src = zeros(1, size(table, 2));
    track{curr}.logprob = zeros(1, size(table, 2));
    track{curr}.emitlogprob = zeros(1, size(table, 2));
    track{curr}.x = zeros(1, size(table, 2));
    track{curr}.y = zeros(1, size(table, 2));
    
    for n=length(start):-1:1
        r = start(n):finish(n);
        f=finish(n);
        if n<length(start)
            idx = path(track{curr}.src(start(n+1)), start(n+1));
            
            track{curr}.src(finish(n)) = idx;
            track{curr}.logprob(f) = table(idx, f);
            track{curr}.emitlogprob(f) = emitlogprobs(idx, f);
        else
            [track{curr}.logprob(finish(n)), track{curr}.src(finish(n))] = max(table(:, finish(n)));
            track{curr}.emitlogprob(finish(n)) = emitlogprobs(track{curr}.src(finish(n)), finish(n));
        end
        for f=finish(n)-1:-1:start(n)
            idx = path(track{curr}.src(f+1), f+1);
            
            track{curr}.src(f) = idx;
            track{curr}.logprob(f) = table(idx, f);
            track{curr}.emitlogprob(f) = emitlogprobs(idx, f);
        end
        if n>1
            track{curr}.logprob(finish(n-1)+1:start(n)-1) = track{curr}.logprob(start(n));
            track{curr}.emitlogprob(finish(n-1)+1:start(n)-1) = track{curr}.emitlogprob(start(n));
        end
    end
    
    %%
    fprintf('#   . computing path\n');
    track{curr}.valid = track{curr}.src <= nClusters & track{curr}.src > 0;
    track{curr}.id = track{curr}.valid * curr;
    
    range = 1:length(track{curr}.valid);
    range = range(track{curr}.valid);
    
    src = track{curr}.src(track{curr}.valid);
    
    track{curr}.x = zeros(1, length(track{curr}.valid));
    track{curr}.x(track{curr}.valid) = scoord{curr}.x(sub2ind(size(scoord{curr}.x), src, range));
    track{curr}.x(~track{curr}.valid) = nan;
    
    track{curr}.y = zeros(1, length(track{curr}.valid));
    track{curr}.y(track{curr}.valid) = scoord{curr}.y(sub2ind(size(scoord{curr}.y), src, range));
    track{curr}.y(~track{curr}.valid) = nan;
    
    %if curr == 1; break;end
end
origTrack = track;
if obj.OutputInOldFormat
    filename = [obj.OutputPath obj.FilePrefix '.raw-track.mat'];
    save(filename, 'track');
end

%% Fix trajectories
fprintf('# - setting additional properties & fixes\n');
for i=1:obj.nSubjects;
    % ignore invalid frames
    nValidFrames = length(track{i}.src);
    track{i}.src(nValidFrames+1:obj.nFrames) = nan;
    track{i}.id(nValidFrames+1:obj.nFrames) = nan;
    track{i}.x(nValidFrames+1:obj.nFrames) = nan;
    track{i}.y(nValidFrames+1:obj.nFrames) = nan;
    
    track{i}.hidden = track{i}.src > nClusters;
    
    % remove skipped areas
    d = diff([0 track{i}.src == 0 0], 1, 2);
    start  = find(d > 0);
    finish = find(d < 0) - 1;
    for j=1:length(start)
        if start(j) <= 1
            continue;
        end
        %track{i}.src(start(j)-1)
        if track{i}.src(start(j)-1) > nClusters
            track{i}.src(start(j):finish(j)) = track{i}.src(start(j)-1);
        else
            track{i}.src(start(j):finish(j)) = nClusters + obj.ROI.nHidden + 1;
        end
    end    

    % set zones coordinates
    valid = track{i}.src > nClusters & track{i}.src <= nClusters + obj.ROI.nHidden;
    track{i}.x(valid) = HiddenSCenters(track{i}.src(valid) - nClusters, 1);
    track{i}.y(valid) = HiddenSCenters(track{i}.src(valid) - nClusters, 2);
    
    % interpolate hidden epochs
    d = diff([0 track{i}.src == 0 | track{i}.src == nClusters + obj.ROI.nHidden + 1  0], 1, 2);
    start  = find(d > 0);
    finish = find(d < 0) - 1;
    for j=1:length(start)
        if start(j) <= 1 || finish(j) >= nValidFrames || isnan(track{i}.x(start(j) - 1)) || isnan(track{i}.x(finish(j) + 1))
            continue;
        end
        r = start(j):finish(j);
        [track{i}.x(r), track{i}.y(r)] = MotionInterp(track{i}.x, track{i}.y, r, [obj.VideoWidth * obj.VideoScale, obj.VideoHeight * obj.VideoScale]);
    end
    
    % setting regions
    x = round(min(max(round(track{i}.x), 1), obj.VideoWidth) / obj.VideoScale);
    y = round(min(max(round(track{i}.y), 1), obj.VideoHeight) / obj.VideoScale);
    track{i}.regions = zeros(1, obj.nFrames);
    for r=find(obj.ROI.IsAvail)
        try
            track{i}.regions(obj.ROI.Regions{r}(sub2ind(size(obj.ROI.Regions{r}), y, x))) = r;
        catch
            pausefr = 'franck';
        end
    end
    track{i}.sheltered = track{i}.src > nClusters & track{i}.src <= nClusters + obj.ROI.nHidden;
    
    % setting zones
    track{i}.zones = ones(1, obj.nFrames);
    obj.ROI.ZoneNames{1} = 'Open';
    index = 2;
    for r=find(obj.ROI.IsAvail)
        idx = track{i}.regions == r & ~track{i}.sheltered;
        track{i}.zones(idx) = index;
        if i==1; obj.ROI.ZoneNames{index} = obj.ROI.RegionNames{r}; end
        if obj.ROI.IsSheltered(r)
            % set long invisible periods as sheltered
            invisibleInRegion = track{i}.regions == r & track{i}.hidden & ~track{i}.sheltered;
            d = diff([0 invisibleInRegion  0], 1, 2);
            start  = find(d > 0);
            finish = find(d < 0) - 1;
            len = finish - start + 1;
            for q = find(len > obj.MaxHidden)
                track{i}.sheltered(start(q):finish(q)) = true;
            end
            %
            idx = track{i}.regions == r & track{i}.sheltered;
            track{i}.zones(idx) = index + 1;
            if i==1; obj.ROI.ZoneNames{index + 1} = ['(' obj.ROI.RegionNames{r} ')']; end
            index = index + 2;
        else
            index = index + 1;
        end
    end
    
    % ignore invalid frames
    track{i}.id(isnan(track{i}.x)) = nan;
    track{i}.zones(isnan(track{i}.x)) = nan;
    track{i}.regions(isnan(track{i}.x)) = nan;
    track{i}.hidden(isnan(track{i}.x)) = true;
    track{i}.sheltered(isnan(track{i}.x)) = false;
end
%%
fprintf('# - finding paths\n');
prevtrack = track;
social.options = obj;
social.dt = obj.dt;
social.x = zeros(length(track), length(track{1}.x), 'single');
social.y = zeros(length(track), length(track{1}.x), 'single');
for i=1:length(track)
    social.x(i, :) = track{i}.x / obj.VideoScale;
    social.y(i, :) = track{i}.y / obj.VideoScale;
end
%
fprintf('# - organizing results\n');
social.nData = size(social.x, 2);
social.time = (1:social.nData) * obj.dt;
zones.all = zeros(obj.nSubjects, obj.nFrames);
zones.hidden = false(obj.nSubjects, obj.nFrames);
zones.sheltered = false(obj.nSubjects, obj.nFrames);
zones.regions = zeros(obj.nSubjects, obj.nFrames);
%
for i=1:obj.nSubjects
    x = round(social.x(i, :)); x(x == 0) = 1;
    y = round(social.y(i, :)); y(y == 0) = 1;
    %
    startnan = find(isnan(x), 1);
    if startnan > 1
        x(startnan:end) = x(startnan - 1);
        y(startnan:end) = y(startnan - 1);
    else
        x(startnan:end) = 1;
        y(startnan:end) = 1;
    end
    %
    zones.hidden(i, :) = track{i}.id(1:obj.nFrames) == 0;
    zones.all(i, :) = track{i}.zones(1:obj.nFrames);
    zones.sheltered(i, :) = track{i}.sheltered(1:obj.nFrames);
    zones.regions(i, :) = track{i}.regions(1:obj.nFrames);
end
social.zones = zones;
social.colors = obj.Colors.Centers;
social.nSubjects = obj.nSubjects;
social.zones.labels = obj.ROI.ZoneNames;

%%
if obj.OutputInOldFormat
    fprintf('# - saving data in old format\n');
    filename = [obj.OutputPath obj.FilePrefix '.social.mat'];
    save(filename, 'social', '-v7.3');
end

%% Store tracking data in object
obj.x = social.x(:, 1:obj.nFrames);
obj.y = social.y(:, 1:obj.nFrames);
try
    obj.time = cents.timestamp;
catch
    obj.time = social.time;
end
obj.zones = social.zones.all;
obj.regions = social.zones.regions;
obj.hidden = social.zones.hidden;
obj.sheltered = social.zones.sheltered;
obj.valid = ~isnan(sum(obj.x, 1)); obj.valid = obj.valid(1:obj.nFrames);

%% removing short events
if ~isfield(obj.ROI, 'Scaling')
    obj.ROI.Scaling = 1/4;
end

for s=1:obj.nSubjects
    [b, e, len] = FindEvents(obj.sheltered(s, :));
    
    % remove short epochs in shelter
    idx = find(len < obj.MinShelteredDuration & b > find(obj.valid, 1) & e < find(obj.valid, 1, 'last'));
    for i=idx
        seq = linspace(obj.x(s, b(i) - 1), obj.x(s, e(i) + 1), len(i) + 2); obj.x(s, b(i):e(i)) = seq(2:end-1);
        seq = linspace(obj.y(s, b(i) - 1), obj.y(s, e(i) + 1), len(i) + 2); obj.y(s, b(i):e(i)) = seq(2:end-1);
        obj.zones(s, b(i):e(i)) = 1;
        obj.sheltered(s, b(i):e(i)) = false;
    end
    
    % remove short observations
    [b, e, len] = FindEvents(~obj.hidden(s, :));
    idx = find(len < obj.MinEventDuration & b > find(obj.valid, 1) & e < find(obj.valid, 1, 'last'));
    for i=idx
        obj.hidden(s, b(i):e(i)) = true;
    end
    
    % remove short hidden events
    [b, e, len] = FindEvents(obj.hidden(s, :));
    idx = find(len < obj.MinEventDuration & b > find(obj.valid, 1) & e < find(obj.valid, 1, 'last'));
    for i=idx
        obj.hidden(s, b(i):e(i)) = false;
    end
    
    % move hidden events that end in the shelter to the shelter
    e = find(obj.hidden(s, 1:end-1) & obj.sheltered(s, 2:end) & ~obj.sheltered(s, 1:end-1));
    b = e * 0;
    for i=1:length(e)
        cb = find(~obj.hidden(s, 1:e(i)), 1, 'last');
        if isempty(cb)
            continue;
        end
        b(i) = cb;
        obj.sheltered(s, b(i):e(i)) = true;
        obj.x(s, b(i):e(i)) = obj.x(s, e(i) + 1);
        obj.y(s, b(i):e(i)) = obj.y(s, e(i) + 1);
        obj.zones(s, b(i):e(i)) = obj.zones(s, e(i) + 1);
    end
end

%% if near shelter and hidden, assume in shelter
hiddenIndices = find(cellfun(@length, regexp(obj.ROI.ZoneNames, '^\(.*\)$')));
for s=1:obj.nSubjects
    %%
    [b,e,len] = FindEvents(obj.hidden(s, :) & ~obj.sheltered(s, :));
    idx = find(len >= max(obj.MinEventDuration, obj.MinShelteredDuration));
    %%
    for i=idx
        range = b(i):e(i);
        for h=1:obj.ROI.nHidden
            %%
            hcoord = single(obj.ROI.HiddenBoundaryCoordinates{h});
            d = pdist2([obj.x(s, range)', obj.y(s, range)'], hcoord);
            minDistance = min(d / obj.PixelsPerCM, [], 2);
            if all(minDistance < obj.MinDistanceToShelterCM)
                obj.sheltered(s, range) = true;
                obj.zones(s, range) = hiddenIndices(h);
            end
        end
    end
end

i = 1;
for f=find(cellfun(@length, regexp(obj.ROI.ZoneNames, '^\(.*\)$')))
    obj.x(obj.zones == f) = obj.ROI.HiddenCenters(i, 1);
    obj.y(obj.zones == f) = obj.ROI.HiddenCenters(i, 2);
    i = i + 1;
end
