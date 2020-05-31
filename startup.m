%% Common environment settings
mystack = dbstack('-completenames'); % use this to determine the basepath
filepath = fileparts(mystack(1).file);
%
addpath(fullfile(filepath, 'Basics'));
addpath(fullfile(filepath, 'Tracking'));
addpath(fullfile(filepath, 'Behavior'));
addpath(fullfile(filepath, 'Preprocessing'));
addpath(fullfile(filepath, 'CheeseSquare'));
addpath(fullfile(filepath, 'MaximumEntropy'));
