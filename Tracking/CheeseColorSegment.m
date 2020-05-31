function [obj, cents, segmented] = CheeseColorSegment(obj, varargin)
% CheeseColorSegment Segments frames in video file using color information
%
%   See documentation for algorithm
%
%   [obj, cents] = CheeseColorSegment(obj) Segments all frames in the 
%   video file referenced by 'obj'. Returns the object 'obj' (in classical
%   CheeseSquare format). The 'cents' variable contains the information
%   about each segment. The 'cents' data is also saved to a file named:
%       (Prefix).segm.(segid).mat, where (prefix) is the video prefix
%       (SC.exp0001.day01.cam02) and (segid) is the part of the video file
%       that was segmented
%
%   CheeseColorSegment(obj, range) Instructs the segmentation algorithm to
%   run on frames in the range between range(1) and range(end)
%
%   CheeseColorSegment(obj, nsegs, segid) Instructs the segmentation algorithm 
%   to divide the video file into 'nsegs' segments, and run on segment
%   number 'segid'
%
%       Created by OREN FORKOSH
%

extFeatsVecs = false; % turn ON/OFF to enable/disable the extraction of the features vectors.
contactDistThresh = 4; % [cm]. If set to -1, ALL the features vectors of ALL the frames will be saved.
                        % If set to an integer, only the features vectors
                        % of mice having a bodyEnd closer than this integer
                        % will be saved.

matchCents = true;
if exist('../Basics', 'dir')
    addpath('../Basics');
end
%% Convert CheeseSquare object to classical
segmented = [];
if ischar(obj)
    obj = CheeseSquare(obj, obj);
    obj = obj.ToClassical;
elseif isa(obj, 'CheeseSquare')
    obj = obj.ToClassical;
elseif extFeatsVecs
%     obj = CheeseSquare(obj, obj);
%     obj = obj.ToClassical;
%     obj.VideoScale = 0.25;
end

%% Segmentation parameters
obj.ErodeRadios = 2; % Erosion radius for each blob
% The background subtraction is thereshold such that the number of segments
% does not exceed 'MaxNumOfCents'
obj.MaxNumOfCents = obj.nSubjects * 3;
if ~extFeatsVecs
    obj.RejectBackground = false; % deprecated
end
obj.UseRGB = false; % deprecated
% Store additional properties of each segment:
% ('Solidity', 'Centroid', 'PixelIdxList', 'MajorAxisLength', 'MinorAxisLength', 'Orientation')
SaveFullSegmentData = false;

%% Histogram of background image
chlabels = {'H', 'S', 'V'};
bkghsv = rgb2hsv(obj.BkgImage);
for ch=1:3
    currch = bkghsv(:, :, ch);
    HistBackground.(chlabels{ch}) = histc(currch(:), obj.Colors.Bins)';
end
HistBackground.Count = obj.Colors.Background.Count + 1;

%% 
nFrames = round(obj.nFrames); % number of frames
nColors = obj.nSubjects + 1; % number of colors (#subjects + 1 for background)
%% Parse arguments
snapshot = []; % Deprecated! don't use...
output = false;
if isempty(varargin)
    startframe = 1;
    endframe = nFrames;
    output = true;
elseif length(varargin) == 1
    if isvector(varargin{1})
        startframe = varargin{1}(1);
        endframe = varargin{1}(end);
        output = true;
    else
        startframe = 0;
        endframe = 1;
        output = false;
        snapshot = varargin{1};
    end
    id = 0;
else
    nruns = varargin{1};
    id = varargin{2};
    step = floor(nFrames / nruns);
    curr = 0;
    for i=1:id
        prev = curr + 1;
        if i == nruns
            curr = nFrames;
        else
            curr = prev + step - 1;
        end
    end
    startframe = prev;
    endframe = curr;
end

nFrames = endframe - startframe + 1;
%% Should match segments between video frames
if matchCents
    realStartframe = startframe;
    if nFrames > 1
        startframe = max(1, startframe - 1);
    end
    nFrames = endframe - startframe + 1;
    cents.prev = zeros(nFrames, obj.MaxNumOfCents, 'uint8');
end
%% Initialize the 'cents' struct which stores data regarding all the segments
%  The struct contains:
%       timestamp:  frame number
%       x:          x-coordinate of segment center
%       y:          y-coordinate of segment center
%       label:      most likely mouse label
%       area:       area of segment
%       solidity:   solidity of segment
%       logprob:    logprob of each segment for every mouse (according to its color)

sx = []; sy = [];

cents.timestamp = zeros(nFrames, 1);
cents.x = zeros(nFrames, obj.MaxNumOfCents);
cents.y = zeros(nFrames, obj.MaxNumOfCents);
cents.label = zeros(nFrames, obj.MaxNumOfCents, 'uint8');
cents.area = zeros(nFrames, obj.MaxNumOfCents, 'uint16');
cents.solidity = zeros(nFrames, obj.MaxNumOfCents, 'single');
cents.logprob = zeros(nFrames, obj.MaxNumOfCents, obj.nSubjects);
cents.threshold = zeros(nFrames, 1, 'uint8');
if extFeatsVecs
%     cents.cents2.abMat = zeros(nFrames*obj.nSubjects, 4);
%     cents.cents2.featsVecsMat = zeros(nFrames*obj.nSubjects, 7626, 'uint8');
    cents.abMat = cell(nFrames,1);
    cents.featsVecsMat = cell(nFrames,1);
    cents.ids = cell(nFrames,1);
end
if SaveFullSegmentData
    Properties = {'Solidity', 'Centroid', 'PixelIdxList', 'MajorAxisLength', 'MinorAxisLength', 'Orientation'};
    cents.minorAxisLength = zeros(nFrames, obj.MaxNumOfCents, 'uint16');
    cents.majorAxisLength = zeros(nFrames, obj.MaxNumOfCents, 'uint16');
    cents.orientation = zeros(nFrames, obj.MaxNumOfCents, 'uint8');
else
    Properties = {'Solidity', 'Centroid', 'PixelIdxList'};
    
end

%% 
cmap = [obj.Colors.Centers; 0 0 0]; % main color of each mouse
if ~extFeatsVecs || true % franck scalingPb
    bkgFrame = im2double(imresize(obj.BkgImage, obj.VideoScale)); % background image
else
    bkgFrame = im2double(imresize(obj.BkgImage, 0.25)); % the "obj" passed to myColorSegment in sniff_franck was a "classical" object which derives from a CheeseObject and its VideoScale is set to 0.25 (hard coded)
end
%% 
if isempty(snapshot)
    vid = MyVideoReader(obj.VideoFile); % video file
end
try
    logfile = fopen([obj.OutputPath obj.FilePrefix '.segm.' sprintf('%03d', id) '.log'], 'w');
    vid.FrameNumber = startframe;
    first = true;
    tim = tic;
    stamp = toc(tim);
    for r=1:nFrames % run over all fraomes
        %% Just some progress report...
        RT = toc(tim) / r * obj.FrameRate;
        if toc(tim) - stamp > .5
            if output
                Console.Reprintf('# - segmenting frame %6d [%d-%d] (%3.1f fps) [%5.1f%%]\n', r+startframe-1, startframe, endframe, r / toc(tim), r / (endframe - startframe) * 100);
            else
                if length(varargin) == 2
                    fprintf(logfile, '# - %3d of %3d: segmenting frame %6d [%d-%d] (%6.2fxiRT)\n', varargin{2}, varargin{1}, r+startframe-1, startframe, endframe, RT);
                    fprintf('# - %3d of %3d: segmenting frame %6d [%d-%d] (%3.1f fps) [%5.1f%%]\n', varargin{2}, varargin{1}, r+startframe-1, startframe, endframe, r / toc(tim), r / (endframe - startframe) * 100);
                else
                    fprintf('# - segmenting frame %6d [%d-%d] (%3.1f fps) [%5.1f%%]\n', r+startframe-1, startframe, endframe, r / toc(tim), r / (endframe - startframe) * 100);
                end
            end
            %         fprintf('# - segmenting frame %6d [%d-%d] (%3.1f fps) [%5.1f%%]\n', r+startframe-1, startframe, endframe, Console.Counter / Console.Timer, r / (endframe - startframe) * 100);
            stamp = toc(tim);
            
        end
        currTime = (r+startframe-1) * obj.dt;
        if currTime < obj.StartTime;
            continue;
        end
        if obj.EndTime > 0 && currTime > obj.EndTime
            continue;
        end
        
        if isempty(snapshot)
            if ~first
                [vid, orig] = vid.NextFrame();
            else
                orig = vid.CurrentFrame;
            end
            if isempty(orig)
                continue;
            end
        else
            orig = snapshot;
        end
        if extFeatsVecs
            frame = orig;
        end
        first = false;
        
        %% Video frame
        if extFeatsVecs && false % franck scalingPb
            orig = im2double(imresize(orig, 0.25));
        else
            orig = im2double(imresize(orig, obj.VideoScale));
        end
        m = imsubtract(orig, bkgFrame); % image subtract with background image
        
        if isempty(sx)
            [sx, sy, nc] = size(m);
        end
        
        %% Find a brightness threshold for the background substraction
        %  Set the thereshold such that the number of segments does not
        %  exceed 'MaxNumOfObjects'
        lum = max(m,[],3);
        meanBKG = mean(lum(:));
        stdBKG = std(lum(:));
        if obj.UseAdaptiveThresh
            upper = obj.NoiseThresh;
            lower = 1;
            prev_thresh = round((upper + lower)/2);
            while true
                thresh = round((upper + lower)/2);
                bw = lum > meanBKG + thresh * stdBKG;
                %bw(ignoreRegion) = false;
                cc = bwconncomp(bw);
                if cc.NumObjects < obj.MaxNumOfObjects
                    upper = thresh - 1;
                    prev_thresh = thresh;
                else
                    lower = thresh + 1;
                end
                if lower > upper
                    break
                end
            end
            if thresh ~= prev_thresh
                thresh = prev_thresh;
                bw = lum > meanBKG + thresh * stdBKG;
            end
        else
            thresh = obj.NoiseThresh;
            bw = lum > meanBKG + thresh * stdBKG;
            cc = bwconncomp(bw);
            if cc.NumObjects > obj.MaxNumOfObjects
                continue;
            end
        end
        bw = bwareaopen(bw, obj.MinNumOfPixels); % remove small objects
        cents.threshold(r) = thresh;
        cents.timestamp(r) = vid.Time;
        if ~any(bw(:)) % no blob was detected
            continue;
        end
        
        %% break frame to color components
        if obj.UseRGB
            rm = orig(:, :, 1);
            gm = orig(:, :, 2);
            bm = orig(:, :, 3);
        else
            hsv = rgb2hsv(orig);
            hm = hsv(:,:,1);
            sm = hsv(:,:,2);
            vm = hsv(:,:,3);
        end
        
        %% Compute log-likelihood of each pixel wrt each mouse
        % the colors histograms were computed during the preprocessing step
        if obj.UseRGB
            [b, idx_r] = histc(rm(bw), obj.Colors.Bins);
            [b, idx_g] = histc(gm(bw), obj.Colors.Bins);
            [b, idx_b] = histc(bm(bw), obj.Colors.Bins);
            
            prob_r = zeros(nColors, length(idx_r));
            prob_g = zeros(nColors, length(idx_g));
            prob_b = zeros(nColors, length(idx_b));
            
            for i=1:nColors
                if i <= obj.nSubjects
                    idx = i;
                    source = obj.Colors.Histrogram;
                else
                    idx = 1;
                    source = obj.Colors.Background;
                end
                prob_r(i, :) = max(source.R(idx, idx_r), 1) / sum(source.R(idx, :));
                prob_g(i, :) = max(source.G(idx, idx_g), 1) / sum(source.G(idx, :));
                prob_b(i, :) = max(source.B(idx, idx_b), 1) / sum(source.B(idx, :));
            end
            joint_prob = prob_r .* prob_g .* prob_b;
        else
            [b, idx_h] = histc(hm(bw), obj.Colors.Bins);
            [b, idx_s] = histc(sm(bw), obj.Colors.Bins);
            [b, idx_v] = histc(vm(bw), obj.Colors.Bins);
            
            prob_h = zeros(nColors, length(idx_h));
            prob_s = zeros(nColors, length(idx_s));
            prob_v = zeros(nColors, length(idx_v));
            
            for i=1:nColors
                if i <= obj.nSubjects
                    idx = i;
                    source = obj.Colors.Histrogram;
                else
                    idx = 1;
                    source = HistBackground;
                end
                if source.Count(idx) > 0
                    prob_h(i, :) = max(source.H(idx, idx_h), 1) / sum(source.H(idx, :) + (source.H(idx, :) == 0));% prob_h(i, :) = prob_h(i, :) / sum(prob_h(i, :));
                    prob_s(i, :) = max(source.S(idx, idx_s), 1) / sum(source.S(idx, :) + (source.S(idx, :) == 0));% prob_s(i, :) = prob_s(i, :) / sum(prob_s(i, :));
                    prob_v(i, :) = max(source.V(idx, idx_v), 1) / sum(source.V(idx, :) + (source.V(idx, :) == 0));% prob_v(i, :) = prob_v(i, :) / sum(prob_v(i, :));
                else
                    prob_v(i, :) = 0;
                    prob_s(i, :) = 0;
                    prob_v(i, :) = 0;
                end
                
            end
            joint_prob = prob_h .* prob_s .* prob_v;
        end
        
        %% assign the most likely mouse id to each pixel 
        [m, idx] = max(joint_prob, [], 1);
        nlabels = zeros(sx,sy,'uint8');
        nlabels(bw) = idx;
        
        logprobmap = cell(1, nColors);
        for k=1:nColors
            logprobmap{k} = zeros(sx,sy);
            logprobmap{k}(bw) = InfoTheory.Log(joint_prob(k, :));
        end
        if output || ~isempty(snapshot)
            flabels = zeros(sx,sy,'uint8');
        end
        
        %% filter small regions
        rejected = nlabels == nColors;
        conn = cell(1, obj.nSubjects);
        validmap = cell(1, obj.nSubjects);
        for i=1:obj.nSubjects
            match = nlabels == i;
            map = bwareaopen(match, obj.MinNumOfPixels);
            rejected(match & ~map) = true;
            
            conn{i} = regionprops(map, {'Solidity', 'Area', 'PixelIdxList'});
            validmap{i} = false(size(bw));
            invalid = [conn{i}.Area] .* [conn{i}.Solidity] < obj.MinNumOfPixels;
            rejected   (cat(1, conn{i}( invalid).PixelIdxList)) = true;
            validmap{i}(cat(1, conn{i}(~invalid).PixelIdxList)) = true;
        end
        
        %% reassign rejected regions to clusters
        lrejected = bwlabel(rejected);
        dlrejected = imdilate(lrejected, ones(3,3));
        for i=1:obj.nSubjects
            u_ = unique(dlrejected(validmap{i}))';
            if ~isempty(u_)
                for u=u_
                    if u > 0
                        validmap{i}(lrejected == u) = true;
                        rejected(lrejected == u) = false;
                        changed = true;
                    end
                end
            end
        end
        
        %% compute the new region props
        for i=1:obj.nSubjects
            validmap{i} = AutoErode(validmap{i}, obj.ErodeRadios);
            reg = bwconncomp(validmap{i});
            conn{i} = regionprops(reg, Properties{:});
        end
        rejected = AutoErode(rejected, obj.ErodeRadios);
        reg = bwconncomp(rejected);
        conn{obj.nSubjects+1} = regionprops(reg, Properties{:});
        currConn = cell(1, obj.MaxNumOfCents);
        
        %% save segments into struct
        if matchCents
            centIndex = 1;
        end
        for i=1:obj.nSubjects+1
            %%
            if i <= obj.nSubjects
                currmap = validmap{i};
            else
                currmap = rejected;
            end
            for j=1:length(conn{i})
                if length(conn{i}(j).PixelIdxList) > obj.MinNumOfPixelsAfterErode
                    cents.x(r, centIndex) = conn{i}(j).Centroid(1);
                    cents.y(r, centIndex) = conn{i}(j).Centroid(2);
                    cents.area(r, centIndex) = length(conn{i}(j).PixelIdxList);
                    cents.solidity(r, centIndex) = conn{i}(j).Solidity;
                    if SaveFullSegmentData
                        cents.minorAxisLength(r, centIndex) = conn{i}(j).MinorAxisLength;
                        cents.majorAxisLength(r, centIndex) = conn{i}(j).MajorAxisLength;
                        cents.orientation(r, centIndex) = conn{i}(j).Orientation;
                    end
                    %% log-likelihood of each segment
                    logsum = -inf;
                    for k=1:nColors
                        %cents.logprob(r, centIndex, k) = sum(logprobmap{k}(conn{i}(j).PixelIdxList));
                        logsum = InfoTheory.LogSum(logsum, logprobmap{k}(conn{i}(j).PixelIdxList));
                    end
                    for k=1:nColors
                        cents.logprob(r, centIndex, k) = log(mean(exp(logprobmap{k}(conn{i}(j).PixelIdxList) - logsum)));
                    end
                    
                    %% most likely mouse id for each segment
                    if i <= obj.nSubjects || (~extFeatsVecs && ~obj.RejectBackground)
                        currlabel = i;
                        cents.label(r, centIndex) = i;
                    else
                        [q, currlabel] = max(cents.logprob(r, centIndex, 1:obj.nSubjects));
                        cents.label(r, centIndex) = currlabel;
                    end
                    
                    %% Match segments between frames
                    if matchCents
                        currConn{centIndex} = conn{i}(j).PixelIdxList;
                        if r > 1
                            a = zeros(1, obj.MaxNumOfCents);
                            for k=find(cents.label(r-1, :) == currlabel)
                                a(k) = length(intersect(prevConn{k}, currConn{centIndex}));
                            end
                            [m, k] = max(a);
                            if m > 0
                                cents.prev(r, centIndex) = k;
                            end
                        end
                    end
                    if (output || ~isempty(snapshot)) && i~= cents.label(r, centIndex)
                        validmap{i}(conn{i}(j).PixelIdxList) = false;
                        validmap{cents.label(r, centIndex)}(conn{i}(j).PixelIdxList) = true;
                    end
                    centIndex = centIndex + 1;
                    if centIndex > obj.MaxNumOfCents
                        break;
                    end
                else
                    currmap(conn{i}(j).PixelIdxList) = 0;
                    if output  || ~isempty(snapshot) || extFeatsVecs
                        validmap{i}(conn{i}(j).PixelIdxList) = false;
                    end
                end
            end
            if centIndex > obj.MaxNumOfCents
                cents.label(r, :) = 0;
                break;
            end
            
        end
        if matchCents
            prevConn = currConn;
        end
        %%
        if  ~isempty(snapshot)
            for i=1:obj.nSubjects
                reg = bwconncomp(validmap{i});
                img = labelmatrix(reg);
                flabels(img > 0) = i;
            end
            rgblbls = label2rgb(flabels, cmap, 'k');
            segmented = rgblbls;
        end
        
        %% extFeatsVecs output the nSubjects features vectors for each frame
        if extFeatsVecs
            
            allabMat = []; allfeatsVecsMat = []; allids = [];
            for i = 1 : obj.nSubjects + 1
                for j=1:length(conn{i})
                    if length(conn{i}(j).PixelIdxList) > obj.MinNumOfPixelsAfterErode
                        actBlobImg = zeros(sx,sy,'uint8');
                        actBlobImg(conn{i}(j).PixelIdxList) = i;
%                         [abMat, featsVecsMat, ~] = HandT_tools.extractFeaturesVectors(actBlobImg, obj.nSubjects, 0.25,...
%                             obj.Scale.ArenaWidth, obj.Scale.ArenaHeight, obj.Scale.ArenaCoord, frame, obj.BkgImage, contactDistThresh);
                        oldFVsExtract = false;
                        if oldFVsExtract % works : DO NOT MODIFY but check actBlobImg(actBlobImg ~= 0) = i;
                            actBlobImg = imresize(actBlobImg, 0.25/obj.VideoScale); % franck scalingPb
                            actBlobImg(actBlobImg ~= 0) = i;                        % franck scalingPb
                            [abMat, featsVecsMat, ~] = HandT_tools.extractFeaturesVectors(actBlobImg, obj.nSubjects, 0.25,...
                                obj.Scale.ArenaWidth, obj.Scale.ArenaHeight, obj.Scale.ArenaCoord, frame, ...
                                obj.BkgImage, contactDistThresh);
                        elseif true
                            [abMat, featsVecsMat, ~] = HandT_tools.extractFeaturesVectors(actBlobImg, obj.nSubjects, obj.VideoScale,...
                                obj.Scale.ArenaWidth, obj.Scale.ArenaHeight, obj.Scale.ArenaCoord, frame, ...
                                obj.BkgImage, contactDistThresh);
                            %{
                        elseif false
                            tmpOurScale = size(obj.BkgImage);
                            ourScale = mean([608, 968] ./ tmpOurScale(1:2));
                            actBlobImg = imresize(actBlobImg, ourScale/obj.VideoScale); % franck scalingPb
%                             actBlobImg = imresize(actBlobImg, obj.VideoScale); % franck scalingPb
                            actBlobImg(actBlobImg ~= 0) = i; 
                            [abMat, featsVecsMat, ~] = HandT_tools.extractFeaturesVectors(actBlobImg, obj.nSubjects, 0.25, ... % ourScale,...
                                obj.Scale.ArenaWidth, obj.Scale.ArenaHeight, [obj.Scale.ArenaCoord(1:2) 604 585], imresize(frame, ourScale), ...
                                imresize(obj.BkgImage,ourScale), contactDistThresh);
                            abMat =  abMat ./ ourScale;
                        else
                            tmpOurScale = size(obj.BkgImage);
                            ourScale = mean([608, 968] ./ tmpOurScale(1:2));
                            actBlobImg = imresize(actBlobImg, ourScale/obj.VideoScale); % franck scalingPb
%                             actBlobImg = imresize(actBlobImg, obj.VideoScale); % franck scalingPb
                            actBlobImg(actBlobImg ~= 0) = i; 
                            [abMat, featsVecsMat, ~] = HandT_tools.extractFeaturesVectors(actBlobImg, obj.nSubjects, 0.25, ... % ourScale,...
                                obj.Scale.ArenaWidth, obj.Scale.ArenaHeight, obj.Scale.ArenaCoord, imresize(frame, ourScale), ...
                                imresize(obj.BkgImage,ourScale), contactDistThresh);
                            abMat =  abMat ./ ourScale;
                                %}
                        end
                        allabMat(end+1,:) = abMat; allfeatsVecsMat(end+1,:) = uint8(featsVecsMat); % allids(end+1,:) = uint8(cents.label(r,i));
                    end
                end
            end
            allids = cents.label(r, cents.label(r,:)~=0);
            
            if contactDistThresh ~= -1
                xScale = obj.Scale.ArenaWidth / obj.Scale.ArenaCoord(3);
                yScale = obj.Scale.ArenaHeight / obj.Scale.ArenaCoord(4);
                inContact = zeros(size(allids));
                for i = 1 : size(allabMat,1)
%                     if ismember(i, inContact)
%                         continue
%                     end
                    if inContact(i) == 1
                        continue
                    end
                    for j = 1 : size(allabMat, 1)
                        if allids(i) == allids(j)
                            continue
                        end
                        candidateMin1 = sqrt(((allabMat(i, 1) - allabMat(j, 1)) * xScale)^2 + ((allabMat(i, 2) - allabMat(j, 2)) * yScale)^2); % between leftS bodyEnds
                        candidateMin2 = sqrt(((allabMat(i, 1) - allabMat(j, 3)) * xScale)^2 + ((allabMat(i, 2) - allabMat(j, 4)) * yScale)^2); % between left and right
                        candidateMin3 = sqrt(((allabMat(i, 3) - allabMat(j, 3)) * xScale)^2 + ((allabMat(i, 4) - allabMat(j, 4)) * yScale)^2); % between rightS bodyEnds
                        candidateMin4 = sqrt(((allabMat(i, 3) - allabMat(j, 1)) * xScale)^2 + ((allabMat(i, 4) - allabMat(j, 2)) * yScale)^2); % between right and left
                        actDist = min([candidateMin1, candidateMin2, candidateMin3, candidateMin4]);
%                         actDist = sqrt(((allabMat(i,1) - allabMat(j,1)) * xScale)^2 + ((allabMat(i,2) - allabMat(j,2)) * yScale)^2);
                        if actDist <=  contactDistThresh
%                             inContact = [inContact, i, j];
                            inContact(i) = 1; inContact(j) = 1;
                            break
                        end
                    end
                end
%                 inContact = sort(unique(inContact));
                allfeatsVecsMat = allfeatsVecsMat(inContact~=0, :);
                allabMat = allabMat(inContact~=0, :);
                allids = allids(inContact~=0);
            end
            
            cents.abMat{r} = allabMat;
            cents.featsVecsMat{r} = uint8(allfeatsVecsMat);
            cents.ids{r} = uint8(allids);
        end       
        
        
        %% Plot the segments
        if output 
            %%
            for ii=find(cents.label(r, :) > 0)
                fprintf('   %d: ', cents.label(r, ii));
                fprintf('%.1f ', exp(squeeze(cents.logprob(r, ii, 1:obj.nSubjects+1))));
                fprintf('\n');
            end
            %%
            for i=1:obj.nSubjects
                reg = bwconncomp(validmap{i});
                img = labelmatrix(reg);
                flabels(img > 0) = i;
            end
            subplot(1,2,2);
            cmap = [CheeseSquare.MiceColors('PRBYW'); lines];
            rgblbls = label2rgb(flabels, cmap, 'k');
            imshow(rgblbls);
            subplot(1,2,1);
            boundries = (imdilate(flabels ~= 0, ones(3,3)) - (flabels ~= 0)) ~= 0;
            for i=1:3
                slice = orig(:, :, i);
                slice(boundries) = 1;
                orig(:, :, i) = slice;
            end
            imagesc(orig);
            title(num2str(r));
            drawnow
        end
        %%
        Console.Counter(1);
    end
    
    if logfile >= 0
        fclose(logfile);
    end
    
catch err
    if exist('id', 'var')
        Console.Warning(err, 'segment no. %d/%d;\n%s', id, nruns, err.message);
    else
        Console.Warning(err);
    end
    %error(err.identifier, 'frame %d in segment no. %d/%d in %s (<a href="matlab: opentoline(%s,%d)">line %d</a>);\n%s', r, id, nruns, err.stack(1).name, err.stack(1).file, err.stack(1).line, err.stack(1).line, err.message);
end

Console.NewLine();
fprintf(['# - total time: ' DateTime.SecToString(toc(tim)) '\n']);
%% fixed matched centers information in 'cents;
if matchCents
    if realStartframe > startframe
        startframe = realStartframe;
        cents.x = cents.x(2:end, :);
        cents.y = cents.y(2:end, :);
        cents.label = cents.label(2:end, :);
        cents.area = cents.area(2:end, :);
        cents.solidity = cents.solidity(2:end, :);
        cents.prev = cents.prev(2:end, :);
        cents.logprob = cents.logprob(2:end, :, :);
        if extFeatsVecs
%             cents.cents2.abMat = cents.cents2.abMat(obj.nSubjects+1:end, :);
%             cents.cents2.featsVecsMat = cents.cents2.featsVecsMat(obj.nSubjects+1:end, :);
            cents.abMat = cents.abMat(2:end, :);
            cents.featsVecsMat = cents.featsVecsMat(2:end, :);
            cents.ids = cents.ids(2:end, :);
        end
    end
end
%% same the segments to file
cents.startframe = startframe;
cents.endframe   = endframe;
%
if ~isempty(javachk('desktop')) || (length(varargin) > 1)
    filename = [obj.OutputPath obj.FilePrefix '.segm.' sprintf('%03d', id) '.mat'];
    fprintf(['# - saving segmentation: "' filename '"']);
    save(filename, 'cents');
end

