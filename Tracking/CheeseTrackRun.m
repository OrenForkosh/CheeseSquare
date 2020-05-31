function obj = CheeseTrackRun(obj, varargin)
% CheeseTrackRun Run the tracking algorithm and behavioral analysis
%   obj = CheeseTrackRun(videofile) runs the tracking on a video
%   file named 'videofile'. Preprocessing file should be located in the
%   same folder with a '.obj.mat' extension. The preprocessing files can be
%   created using the Preprocessing module.
%
%   obj = CheeseTrackRun(cheesesquare) runs the tracking on CheeseSquare
%   object. CheeseSquare object should contain preprocessing data.
%
%   obj = CheeseTrackRun(... , options) optional parameters:
%       recover (default: false) - run segmentation only if needed
%       nSegs   (default: 100) - number of segments to divide the video 
%               file for tracking (for parallelization)
%
%       Created by OREN FORKOSH
%

%% parse additional options
p = inputParser;
addOptional(p, 'recover', false, @islogical); % if recover is set to true, doesn't repeat segmentation unless needed
addOptional(p, 'nSegs',   100,   @isnumeric); % number of segments to divide the video file for tracking (for parallelization)
p.parse(varargin{:});
opt = p.Results; % contain all options

%% ensure 'Basics' is in path
if exist('../Basics', 'dir')
    addpath('../Basics');
end

%% create CheeseSquare object
obj = CheeseSquare(obj);

%% convert CheeseSquare object
% convert CheeseSquare object to previous object type (some of the code 
% still runs on the older version of CheeseSquare)
track = obj.ToClassical; % track.VideoFile = 'Z:\Franck_round3_tmp\fromNoaSansPlouf\Fmr1_pilot_WT2.exp0003.day03.cam03.avi'; % franck

%% Set video size scaling 
% Reduces video resolution to increase tracking speed
if isfield(track, 'Scale') && ~isempty(track.Scale)
    wscale = (track.Scale.ArenaCoord(3) / track.Scale.ArenaWidth) / 9;
    hscale = (track.Scale.ArenaCoord(4) / track.Scale.ArenaHeight) / 9;
else
    wscale = 1;
    hscale = 1;
end
% factor by which to reduce video resolution:
track.VideoScale = track.VideoScale / min(wscale, hscale); 

%% Set the Colors.Centers field
% each mouse is assigned a different color (for visualization) by taking
% its average color
track.Colors.Centers = zeros(track.nSubjects, 3);
for i=1:track.nSubjects
    track.Colors.Centers(i, :) = min(mean(obj.Colors.Marks.color(obj.Colors.Marks.id==i, :)/255), 1);
end

%% Copy tracking data 
% Copy tracking data (if exists; i.e. if file was previously tracked)
try
    if ~isempty(obj.Tracking)
        f = fields(obj.Tracking);
        for i=1:length(f)
            track.(f{i}) = obj.Tracking.(f{i});
        end
        track.time = 1/obj.Video.FrameRate * (1:length(track.x));
    end
catch err
    fprintf('%s', Console.StrudelError(err));
end

%% Starting Tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Video segmentation on all frames
% Most time consuming process therefore the video file is split into
% opt.nSegs parts and the segmentation is done seperatly on each part
ns = opt.nSegs;
Console.Message('segmenting video frames (using color information)');
Console.Reprintf(0, '');
h = tic;
segmfile = [track.OutputPath track.FilePrefix '.segm.%03d.mat'];
Console.Timer(tic);
Console.Counter(nan);
parfor i=1:ns % parallalize process
    if opt.recover && exist(sprintf(segmfile, i), 'file') %#ok<PFBNS>
        % if opt.recover is true and segementation file already exists,
        % don't run segmentation again. This should save time.
        try
            m = matfile(sprintf(segmfile, i)); %#ok<NASGU>
            Console.Message(1, 'skipping segment no. %d', i);
        catch
            CheeseColorSegment(track, ns, i);
        end
    else
        CheeseColorSegment(track, ns, i);
%         CheeseColorSegmentBlackFur(track, ns, i);
    end
end
Console.Message(1, 'segmentation took %.1f seconds', toc(h));

%% Creates the movies from the features vectors if those were computed (one movie per mouse).
%  in order to reduce memory costs.
filename = ['Tracking\' track.OutputPath track.FilePrefix '.segm.' sprintf('%03d', 1) '.mat'];
checkExistenceOfFeaturesVectors = load(filename);
checkExistenceOfFeaturesVectors = checkExistenceOfFeaturesVectors.cents;
if isfield(checkExistenceOfFeaturesVectors, 'featsVecsMat') % && false % before CheeseTrackPath invest (2nd condition...)
    outputPath = '';
    splitVideoFilePath = split(track.VideoFile, '\');
    for j = 1 : length(splitVideoFilePath)-1
        outputPath = strcat(outputPath, splitVideoFilePath{j}, '\');
    end
    HandT_tools.makeMoviesFromSeg(track.OutputPath, outputPath, track.FilePrefix, obj.nSubjects, track.nFrames)
end

%% Find the trails of all mice
% uses the segments from the previous step to infer the actual path of each
% mouse
try
%     if exist('outputPath', 'var'); track.VideoScale = 0.25; end % franck scalingPb
    track = CheeseTrackPath(track, opt.nSegs);
    
    trackfields = {'x',      'y',      'zones',  'hidden',  'sheltered', 'valid'};
    trackclass  = {'single', 'single', 'uint16', 'logical', 'logical',   'logical'};
    for i=1:length(trackfields)
        if isfield(track, trackfields{i})
            obj.Tracking.(trackfields{i}) = cast(track.(trackfields{i}), trackclass{i});
        end
    end
    obj.Profile = 'no profile';
    obj.Save;
catch err
    fprintf('\n%s', Console.StrudelError(err));
end

try
    TrackExport(track);
catch err
    fprintf('%s', Console.StrudelError(err));
end

%% Starting Behavioral Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute group structure (or Hierarchy)
% On the way this also classifies all the contacts between the mice and
% infer Chase-Escape interactions
try
    obj = CheeseHierarchyGroup(obj);
catch err
    fprintf('%s', Console.StrudelError(err));
end

%% Compute behavioral profile for each mouse (and group)
% Computes a table which contains a high-dimensional vector of behaviors
% for each mouse (Such as: No. of approaches, Time in the nest, etc.)
try
    obj = CheeseNewProfiler(obj);
    obj.Save;
catch err
    fprintf('%s', Console.StrudelError(err));
end
