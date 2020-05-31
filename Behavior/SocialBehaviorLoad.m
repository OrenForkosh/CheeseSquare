function marks = SocialBehaviorLoad(obj, map)
% Load manually labeled behavioral data. Assumes that data is located under
% the 'aux' folder
if ~exist('map', 'var')
    map = -1;
end

data = dlmread(['data/' obj.FilePrefix '.marks'], '\t', [1 0 500 15]);
if length(map) == 1 && map(1) == -1
else
    data = data(map, :);
end
marks.ids = data(:, 1);
marks.beg = data(:, 2);
marks.end = data(:, 3);
marks.subjects = data(:, 4:5);
marks.approach = data(:, 6:7);
marks.leave = data(:, 8:9);
marks.chase = data(:, 10);
marks.pred = data(:, 11:12);
marks.prey = data(:, 13:14);
marks.artifact = data(:, 15);
marks.agresivness = data(:, 16);
