% Data regarding past experiments. Deprecated!!! Don't use...

experiments = {...
    '5mice.exp0002.day%02d.cam04', ...
    '5mice.exp0003.day%02d.cam04', ...
    '5mice.exp0004.day%02d.cam04', ...
    ...
    };
%'Enriched.exp0006.day%02d.cam01',...

%     'OTRFloxCtrl.exp0001.day%02d.cam04', ...
%     'OTRFloxCtrl.exp0002.day%02d.cam04', ...
%     'OTRFloxCtrl.exp0003.day%02d.cam01', ...
%     'OTRFloxCtrl.exp0007.day%02d.cam01', ...
%     'OTRFloxCtrl.exp0008.day%02d.cam01', ...
%     ...
%     'OTRFloxTreat.exp0001.day%02d.cam01', ...
%     'OTRFloxTreat.exp0002.day%02d.cam01', ...
%     'OTRFloxTreat.exp0003.day%02d.cam04', ...
%     'OTRFloxTreat.exp0007.day%02d.cam04', ...
%     'OTRFloxTreat.exp0010.day%02d.cam04'


%'OTRFloxCtrl.exp0009.day%02d.cam01', ... % recombined
    %'OTRFloxCtrl.exp0008.day%02d.cam04', ... % repeated group
    %'OTRFloxTreat.exp0009.day%02d.cam04', ... % recombined
    %'OTRFloxTreat.exp0008.day%02d.cam01',...  % repeated group

GroupsMeta.Titles = {...
    '5mice',     'Five mice', ...
    'SC',           'Standard', ...
    'OTRFloxCtrl',  'OTR Flox Ctrl', ...
    'OTRFloxTreat', 'OTR Flox Treat', ...
    };

GroupsMeta.Initials = {...
    '5mice',     '5mice', ...
    'SC',           'Std', ...
    'OTRFloxCtrl',  'OTRCtrl', ...
    'OTRFloxTreat', 'OTRTreat', ...
    };

GroupsMeta.Colors = {...
    '5mice',     [80 171 210] / 255, ...
    'SC',           [232 24 46] / 255, ...
    'OTRFloxCtrl',  [120 182 83] / 255, ...
    'OTRFloxTreat', [237 145 61] / 255, ...
    };


nDays = 4;
nSubjects = 5;

%%
seq = 1:length(experiments);
GroupsData = struct();
GroupsData.Types = {};
for i=1:length(experiments)
    type = regexprep(experiments{i}, '^([^\.]*)\..*', '$1');
    if ~any(strcmp(type, GroupsData.Types))
        GroupsData.Types = {GroupsData.Types{:}, type};
    end
    GroupsData.group(i) = find(strcmp(type, GroupsData.Types), 1);
end
GroupsData.nExperiments = length(experiments);
GroupsData.All = 1:GroupsData.nExperiments;
Groups = struct();
for g=unique(GroupsData.group(:))'
    Groups(g).map = false(1, length(experiments));
    Groups(g).map(GroupsData.group == g) = true;
    Groups(g).idx = find(Groups(g).map);
    Groups(g).color = GroupsMeta.Colors{find(strcmp(GroupsData.Types{g}, GroupsMeta.Colors)) + 1};
    Groups(g).title = GroupsMeta.Titles{find(strcmp(GroupsData.Types{g}, GroupsMeta.Titles)) + 1};
    Groups(g).initials = GroupsMeta.Initials{find(strcmp(GroupsData.Types{g}, GroupsMeta.Initials)) + 1};
end
GroupsData.Proportions.Region = 'Block';
GroupsData.Proportions.Size = [7.7 11.4]; % [cm]

%%
% E.map = logical([1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0]);
% E.idx = find(E.map);
% 
% SC.map = ~E.map;
% SC.idx = find(SC.map);

% Groups(1).map = 1 <= seq & seq <= 10;
% Groups(1).idx = find(Groups(1).map);
% Groups(1).title = 'Social Challange';
% Groups(1).initials = 'SC';
% Groups(1).color =  [80 171 210] / 255;
% 
% Groups(2).map = 11 <= seq & seq <= 18;
% Groups(2).idx = find(Groups(2).map);
% Groups(2).title = 'Standard';
% Groups(2).initials = 'STD';
% Groups(2).color =  [232 24 46] / 255;
% 
% Groups(3).map = 19 <= seq & seq <= 25;
% Groups(3).idx = find(Groups(3).map);
% Groups(3).title = 'OTR Flox Ctrl';
% Groups(3).initials = 'OTRFC';
% Groups(3).color =  [120 182 83] / 255;
% 
% Groups(4).map = 26 <= seq & seq <= 31;
% Groups(4).idx = find(Groups(4).map);
% Groups(4).title = 'OTR Flox Treat';
% Groups(4).initials = 'OTRFT';
% Groups(4).color =  [237 145 61] / 255;
% 
% yellow: 120 182 83

% titles = {'Challanged', 'Standard', 'OTR Flox Ctrl', 'OTR Flox Treat'};

GroupsData.Experiments = experiments;
%GroupsData.Initials = {};
GroupsData.Titles = {};
GroupsData.nGroups = length(Groups);
GroupsData.Colors = [];
GroupsData.index = [];
GroupsData.group = [];
for i=1:length(Groups)
    GroupsData.Titles{i} = Groups(i).title;
    GroupsData.Initials{i} = Groups(i).initials;
    GroupsData.Colors(i, :) = Groups(i).color;
    name = regexprep(Groups(i).title, ' ', '');
    GroupsData.(name) = Groups(i);
    GroupsData.index = [GroupsData.index, Groups(i).idx];
    GroupsData.group = [GroupsData.group, Groups(i).idx * 0 + i];
end

GroupsData.nDays = nDays;
GroupsData.nSubjects = nSubjects;
% if isunix
    GroupsData.OutputPath = 'Res/';
    GroupsData.MoviesPath = '../Movies/';
% else
%     GroupsData.OutputPath = 'Z:/tests/SocialMice/base/Res/';
%     GroupsData.MoviesPath = 'Z:/tests/SocialMice/Movies/';
% end
