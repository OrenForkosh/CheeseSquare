function obj = CheeseTrack(obj, varargin)
% CheeseTrack Run the tracking algorithm
%   obj = CheeseTrack(videofile) runs the tracking on a video
%   file named 'videofile'. Preprocessing file should be located in the
%   same folder with a '.obj.mat' extension. The preprocessing files can be
%   created using the Preprocessing module.
%
%   obj = CheeseTrack(cheesesquare) runs the tracking on CheeseSquare
%   object. CheeseSquare object should contain preprocessing data.
%
%   obj = CheeseTrack(... , options) optional parameters:
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
track = obj.ToClassical;

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
    end
end
Console.Message(1, 'segmentation took %.1f seconds', toc(h));
toc(h)

%% Find the trails of all mice
% uses the segments from the previous step to infer the actual path of each
% mouse
try
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

