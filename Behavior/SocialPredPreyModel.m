function obj = SocialPredPreyModel(obj, args)
% Segment and classify contacts between the mice. Each contact is divides
% to what event preceeded it and how it was concluded. Uses Hidden Markov
% model to infer the state of each subject during the different stages of
% the interaction. The topology of the model:
%
%               (pred)        (pred)
%              /      \      /      \
%        (idel) ------ (cont) ------ (idel)
%              \      /      \      /
%              (prey)        (prey)
%
% See documantation.
%
%       Created by OREN FORKOSH
%

fprintf('# Predator-Prey model\n');
%% Hidden Markove model states:
IdelState = 1;      % No interaction
PrePredState = 2;   % Before contact, subject moved in the direction of the other
PrePreyState = 3;   % Before contact, subject moved away from other
ContactState = 4;   % Contact between the mice
PostPredState = 5;  % After contact, subject moved in the direction of the other
PostPreyState = 6;  % After contact, subject moved away from other
%%
obj = CheeseObject.Load(obj);
%%
if ~exist('args', 'var')
    args = struct();
end

% defaults.PredObjectLength = 75;
% defaults.PreyObjectLength = 150;
% defaults.MinSubInteractionDuration = 5;
% defaults.SmoothSpan = 3;
% defaults.ZoneOfProximity = 200;
% defaults.ZoneOfContact = 100;
% defaults.MinNumOfContacts = 3;
% defaults.MinContanctDuration = 6; % minimal duration of a contact
% defaults.MaxShelteredContanctDuration = 25;

defaults.PredObjectLengthCM = 10;
defaults.PreyObjectLengthCM = 10;
defaults.MinSubInteractionDuration = 3;
defaults.ZoneOfProximityCM = 10; % 10?
defaults.ZoneOfContactCM = 10;

defaults.SmoothSpan = 7;
defaults.MinNumOfContacts = 1;
defaults.MinContanctDuration = 10; % minimal duration of a contact
defaults.MaxInterContanctDuration = 10; % minimal duration of a contact

defaults.MaxShelteredContanctDuration = 15;
defaults.ChaseMinOverlap = 6;
defaults.ChaseMinOverlapPercentage = 0;
defaults.RelativeSpeedInChase = 0;
defaults.MismatchProbability = .25;

defaults.MinPreySpeedCMperSec = 0;
defaults.MinPredSpeedCMperSec = 0;
defaults.EndInNestMaxLeg = 1;

defaults.MinShelteredDuration = 6;
defaults.MinEventDuration = 4;
defaults.MinDistanceToShelter = 4;

defaults.AggresivenessPredMaxSpeedThreshold = 20;
defaults.AggresivenessMaxDistanceThreshold = 10;

defaults.MinPrePostInteractionDistanceCM = 5;

args = CheeseObject.AuxParseStructArguments(defaults, args);

f = fieldnames(args);
for i=1:length(f)
    obj.Interactions.(f{i}) = args.(f{i});
end
%% Find all the interactions
[obj, local] = SocialFindInteractions(obj);

%%
converge = zeros(obj.nSubjects, obj.nSubjects, obj.nFrames);
for m1=1:obj.nSubjects
    for m2=1:obj.nSubjects
        if m1 == m2; continue; end;
        q = tan(obj.Interactions.angle{m1, m2}) .* obj.Interactions.distanceCM{min(m1, m2), max(m1, m2)};
        converge(m1, m2, q > 0 & abs(q) < obj.Interactions.PredObjectLengthCM) =  1;
        converge(m1, m2, q < 0 & abs(q) < obj.Interactions.PreyObjectLengthCM) = -1;
    end
end
converge(isnan(converge)) = 0;

%%
speed = obj.speed(:);
minspeed = FindProbabilityCutoffPoint(speed(speed>0), .75);

%%
features = cell(obj.nSubjects);
index = 1;
fprintf('# - processing subjects ');
%ProgressReport();
for m1=1:obj.nSubjects
    for m2=1:obj.nSubjects
        if m1 == m2; continue; end;
        %ProgressReport(index, obj.nSubjects * (obj.nSubjects - 1));
        %%
        others = true(1, obj.nSubjects);
        others(m1) = 0;
        others(m2) = 0;
        sheltered = local.sheltered(m1, :) | local.sheltered(m2, :);
        %%
        features{m1, m2} = ones(1, obj.nFrames) * 2;
        features{m1, m2}(converge(m1, m2, :) == -1) = 1;
        features{m1, m2}(converge(m1, m2, :) ==  0) = 2;
        features{m1, m2}(converge(m1, m2, :) ==  1) = 3;
        features{m1, m2}(local.sheltered(m1, :) | local.sheltered(m2, :)) = 2;
        features{m1, m2}(obj.speed(m1, :) < minspeed) = 2;
        
        features{m1, m2}(obj.Interactions.proximity(m1, m2, :)) = features{m1, m2}(obj.Interactions.proximity(m1, m2, :)) + 3;
        
        %%
        %
        %        (pred)        (pred)
        %       /      \      /      \
        % (idel) ------ (cont) ------ (idel)
        %       \      /      \      /
        %        (prey)        (prey)
        %
        
        x = .5;
        m = defaults.MismatchProbability;
        small = 0.0;
        model.emis = [...
            m   x   m   0.0 0.0 0.0; % idel
            m   m   x   m   m   x; % pred
            x   m   m   x   m   m; % prey
            0.0 0.0 0.0 m   x   m; % cont
            m+small   m+small   x+small   m+small   m+small   x+small; % pred the eps is added for symmetry breaking with the pre contact events
            x+small   m+small   m+small   x+small   m+small   m+small; % prey
            ];
        
        model.trans = [...
            1 1 1 1 0 0;
            0 1 0 1 0 0;
            0 0 1 1 0 0;
            1 0 0 1 1 1;
            1 0 0 0 1 0;
            1 0 0 0 0 1;
            ];
        model.names = {'-', 'pred', 'prey', 'cont', 'pred', 'prey'};
        obj.Interactions.PredPrey.model = model;
        
        states = hmmviterbi(features{m1, m2},model.trans,model.emis);
        
        %% connect proximal events
        fittrans = model.trans;
        fittrans(Q.exclude(size(model.trans, 1), IdelState), IdelState) = 0;

        segs = Segs(states ~= IdelState, states);
        segs = segs.Close(defaults.MaxInterContanctDuration);
        for i=1:length(segs.Events)
            if any(diff(segs.Events(i).data) < 0)
                states(segs.Events(i).beg:segs.Events(i).end) = hmmviterbi(features{m1, m2}(segs.Events(i).beg:segs.Events(i).end), fittrans,model.emis);
            end
        end
        
        %% remove short events
        % remove post-pre short events
        [begF, endF, len] = FindEvents(states, ...
            states ~= IdelState & states ~= ContactState);
        for i=find(len < obj.Interactions.MinSubInteractionDuration)
            states(begF(i):endF(i)) = ContactState;
        end
        % remove post-pre short in distnace events 
        segs = Segs(states ~= IdelState & states ~= ContactState, [local.xCM(m1, :); local.yCM(m1, :)]');
        distance = segs.Run(@(b,e,d) sum(sqrt(sum((d(2:end, :) - d(1:end-1, :)).^2, 2))));
        for i=find(Q.torow(distance) < defaults.MinPrePostInteractionDistanceCM)
            range = segs.Events(i).beg:segs.Events(i).end;
            contact = obj.Interactions.contact(m1, m2, range);
            states(range(contact)) = ContactState;
            states(range(~contact)) = IdelState;
        end
    
        % remove short contacts
        [begF, endF, len] = FindEvents(states, ...
            states == ContactState);
        for i=find(len < obj.Interactions.MinSubInteractionDuration & begF > 1 & endF < obj.nFrames)
            % all chase events
            if ...
                    (states(begF(i)-1) == PrePredState && states(endF(i)+1) == PostPredState) || ...
                    (states(begF(i)-1) == PrePreyState && states(endF(i)+1) == PostPreyState)
                mid = floor((endF(i) + begF(i))/2);
                states(begF(i):mid) = states(begF(i)-1);
                states(mid+1:endF(i)) = states(endF(i)+1);
            % all chase events
%             if ...
%                     (states(begF(i)-1) == PrePredState  || states(begF(i)-1) == PrePreyState) && ...
%                     (states(endF(i)+1) == PostPredState || states(endF(i)+1) == PostPreyState)
                states(find(states(1:begF(i)) == IdelState, 1, 'last')+1:endF(i)) = states(endF(i)+1);
            end
        end
        %%
        [begF, endF, len, events] = FindEvents(states, states > 1);
        valid = true(1, length(events));
        for i=1:length(events)
            %%
            events{i}.desc = [];
            events{i}.bounds = [];
            range = events{i}.BegFrame:events{i}.EndFrame;
            nContacts = sum(obj.Interactions.contact(m1, m2, range));
            valid(i) = nContacts >= obj.Interactions.MinNumOfContacts;
            if ~valid(i)
                continue;
            end
            events{i}.open = any(obj.zones(m1, range) == obj.zones(m2, range) & obj.zones(m1, range) == 1);
            valid(i) = any(obj.zones(m1, range) == obj.zones(m2, range) & obj.zones(m1, range) == 1);
            if ~valid(i)
                continue;
            end
            %             valid(i) = valid(i) && all(sum(obj.Interactions.contact(m2, others, range), 3) < nContacts);
            %             if ~valid(i)
            %                 continue;
            %             end
            if events{i}.data(1) == 4
                f = find(obj.Interactions.contact(m1, m2, range), 1);
                events{i}.data = events{i}.data(f:end);
                events{i}.BegFrame = events{i}.BegFrame + f - 1;
                range = range(f:end);
            end
            if events{i}.data(end) == 4
                f = find(obj.Interactions.contact(m1, m2, range), 1, 'last');
                events{i}.data = events{i}.data(1:f);
                events{i}.EndFrame = events{i}.BegFrame + f - 1;
                range = range(1:f);
            end
            
            for j=1:length(model.names)
                r = GetRange(events{i}.data == j);
                if ~isempty(r)
                    if ~isempty(events{i}.desc)
                        events{i}.desc = [events{i}.desc ' -> '];
                    end
                    %            events{i}.desc = [events{i}.desc, model.names{j} '(' num2str(r(2)-r(1)+1) ')'];
                    events{i}.desc = [events{i}.desc, model.names{j} '(' num2str(events{i}.BegFrame -1 + r(1)) ')'];
                    events{i}.bounds = [events{i}.bounds; events{i}.BegFrame -1 + [r(1), r(2)]];
                end
            end
        end
        events = events(valid);
        obj.Interactions.PredPrey.states{m1, m2} = states;
        obj.Interactions.PredPrey.events{m1, m2} = events;
        index = index + 1;
    end
end

%%
%obj.Contacts.Behaviors = struct();

%%
obj = SocialMatchEvents(obj, obj.Interactions.PredPrey.events);

states = cat(3, obj.Contacts.List.states);
ChaseEscape = squeeze(any(states(:, PostPredState, :)) & any(states(:, PostPreyState, :)));
Chaser = zeros(1, length(ChaseEscape), 'int8');
Escaper = zeros(1, length(ChaseEscape), 'int8');
obj.Contacts.Properties = struct();
for i=find(ChaseEscape(:)')
    curr = obj.Contacts.List(i);
    if all(curr.end < obj.nFrames)
        obj.Contacts.Properties.EndInNest(i) = max(obj.sheltered(curr.subjects(1), curr.end(1) + 1:curr.end(1) + obj.Interactions.EndInNestMaxLeg)) + ...
            max(obj.sheltered(curr.subjects(2), curr.end(2) + 1:curr.end(2) + obj.Interactions.EndInNestMaxLeg));
    else
        obj.Contacts.Properties.EndInNest(i) = 0;
    end
    %%
    pred = find(curr.states(:, PostPredState));
    prey = 3 - pred;
    %%
    [b1, e1] = FindEvents(curr.data{pred} == PostPredState);
    [b2, e2] = FindEvents(curr.data{prey} == PostPreyState);
    b1 = b1 + curr.beg(pred) - 1;
    e1 = e1 + curr.beg(pred) - 1;
    b2 = b2 + curr.beg(prey) - 1;
    e2 = e2 + curr.beg(prey) - 1;
    %%
    begf = max(b1, b2);
    endf = min(e1, e2);
    %%
    predSpeedCMperSec = obj.speedCMperSec(curr.subjects(pred), begf:endf);
    preySpeedCMperSec = obj.speedCMperSec(curr.subjects(prey), begf:endf);
    %%
    InteractionLength = max(e1 - b1 + 1, e2 - b2 + 1);
    
    if endf - begf + 1 >= obj.Interactions.ChaseMinOverlap && (endf - begf + 1) / InteractionLength > obj.Interactions.ChaseMinOverlapPercentage
        ChaseEscape(i) = median(predSpeedCMperSec./preySpeedCMperSec) > obj.Interactions.RelativeSpeedInChase;
    else
        ChaseEscape(i) = false;
    end
    if ChaseEscape(i)
        %%
        Chaser(i) = curr.subjects(pred);
        Escaper(i) = curr.subjects(prey);
        obj.Contacts.Properties.Speed(:, i) = [median(obj.speedCMperSec(curr.subjects(prey), curr.beg(prey):curr.end(prey))) ...
            median(obj.speedCMperSec(curr.subjects(pred), curr.beg(pred):curr.end(pred)))];
        allbeg = max(curr.beg);
        allend = min(curr.end);
        obj.Contacts.Properties.RelativeSpeed(i) = median(obj.speedCMperSec(curr.subjects(prey), allbeg:allend) ./ ...
            obj.speedCMperSec(curr.subjects(pred), allbeg:allend));
%          obj.Contacts.Properties.MaxRelativeSpeed(i) = maxfinite(obj.speedCMperSec(curr.subjects(prey), allbeg:allend) ./ ...
%              obj.speedCMperSec(curr.subjects(pred), allbeg:allend));
        obj.Contacts.Properties.Duration(:, i) = [curr.end(prey)-curr.beg(prey) + 1; curr.end(pred)-curr.beg(pred) + 1];
        obj.Contacts.Properties.PostSpeed(:, i) = [median(preySpeedCMperSec) median(predSpeedCMperSec)];
        obj.Contacts.Properties.PostSpeedJitter(:, i) = [RobustStd(preySpeedCMperSec') RobustStd(predSpeedCMperSec')];
        obj.Contacts.Properties.PostRelativeSpeed(i) = median(predSpeedCMperSec./preySpeedCMperSec);
%          obj.Contacts.Properties.PostMaxRelativeSpeed(i) = maxfinite(predSpeedCMperSec./preySpeedCMperSec);
        obj.Contacts.Properties.PostMaxSpeed(:, i) = [max(preySpeedCMperSec) max(predSpeedCMperSec)];
        obj.Contacts.Properties.PostHighSpeed(:, i) = [quantile(preySpeedCMperSec, .9) quantile(predSpeedCMperSec, .9)];
        obj.Contacts.Properties.PostDuration(i) = endf - begf + 1;
%         obj.Contacts.Properties.PostPathLength(i) = sum(sqrt(diff(obj.x(curr.subjects(prey), begf:endf)).^2 + diff(obj.y(curr.subjects(prey), begf:endf)).^2)) + ...
%             sum(sqrt(diff(obj.x(curr.subjects(pred), begf:endf)).^2 + diff(obj.y(curr.subjects(pred), begf:endf)).^2));
%         obj.Contacts.Properties.PathLength(i) = sum(sqrt(diff(obj.x(curr.subjects(prey), allbeg:allend)).^2 + diff(obj.y(curr.subjects(prey), allbeg:allend)).^2)) + ...
%             sum(sqrt(diff(obj.x(curr.subjects(pred), allbeg:allend)).^2 + diff(obj.y(curr.subjects(pred), allbeg:allend)).^2));
        distanceCM = obj.Interactions.distanceCM{min(curr.subjects(pred), curr.subjects(prey)), max(curr.subjects(pred), curr.subjects(prey))}(begf:endf);
        obj.Contacts.Properties.PostDistance(i) = [median(distanceCM)];
        obj.Contacts.Properties.PostMaxDistance(i) = [max(distanceCM)];
        obj.Contacts.Properties.PostMinDistance(i) = [min(distanceCM)];
        distanceCM = obj.Interactions.distanceCM{min(curr.subjects(pred), curr.subjects(prey)), max(curr.subjects(pred), curr.subjects(prey))}(allbeg:allend);
        obj.Contacts.Properties.Distance(i) = [median(distanceCM)];
        obj.Contacts.Properties.MaxDistance(i) = [max(distanceCM)];
        obj.Contacts.Properties.MinDistance(i) = [min(distanceCM)];
    end
end
%%
%obj.Contacts.Behaviors.ChaseEscape = struct();
obj.Contacts.Behaviors.ChaseEscape.Map = ChaseEscape(:)';
obj.Contacts.Behaviors.ChaseEscape.Chaser = Chaser;
obj.Contacts.Behaviors.ChaseEscape.Escaper = Escaper;

%obj.Contacts.Behaviors.AggressiveChase = struct();
obj.Contacts.Behaviors.AggressiveChase.Map = []; %obj.Contacts.Properties.PostMaxSpeed(2, :) > 20 & obj.Contacts.Properties.PostMaxDistance < 10;

obj.Contacts.Behaviors.AggressiveChase.Chaser = zeros(1, length(obj.Contacts.Behaviors.AggressiveChase.Map));
obj.Contacts.Behaviors.AggressiveChase.Chaser(obj.Contacts.Behaviors.AggressiveChase.Map) = ...
    obj.Contacts.Behaviors.ChaseEscape.Chaser(obj.Contacts.Behaviors.AggressiveChase.Map);
obj.Contacts.Behaviors.AggressiveChase.Escaper = zeros(1, length(obj.Contacts.Behaviors.AggressiveChase.Map));
obj.Contacts.Behaviors.AggressiveChase.Escaper(obj.Contacts.Behaviors.AggressiveChase.Map) = ...
    obj.Contacts.Behaviors.ChaseEscape.Escaper(obj.Contacts.Behaviors.AggressiveChase.Map);
obj.Interactions = rmfield(obj.Interactions, {'contact', 'proximity', 'distanceCM', 'angle', 'valid'});
%% build classifier for Aggressive-Chase
vec = [];
labels = {};
obj.Contacts.Behaviors.AggressiveChase.Classifier.Vec = [];
f = fieldnames(obj.Contacts.Properties);
sz = inf;
for i=1:length(f)
    sz = min(sz, size(obj.Contacts.Properties.(f{i}), 2));
end
lidx = 1;
for i=1:length(f)
    vec = [vec; obj.Contacts.Properties.(f{i})(:, 1:sz)];
    if size(obj.Contacts.Properties.(f{i}), 1) > 1
        for j=1:size(obj.Contacts.Properties.(f{i}), 1)
            labels{lidx} = [f{i} num2str(j)];
            lidx = lidx + 1;
        end
    else
        labels{lidx} = f{i};
        lidx = lidx + 1;
    end        
end
obj.Contacts.Behaviors.AggressiveChase.Classifier.Vec = vec;
obj.Contacts.Behaviors.AggressiveChase.Classifier.LabelNames = labels;
obj.Contacts.Behaviors.AggressiveChase.Classifier.ChaseEscapeMap = obj.Contacts.Behaviors.ChaseEscape.Map(1:size(vec, 2));
%%
fprintf('\n');
try
    % TrackSave(obj);
catch
end

%if isdesktop
    try
        %%
        data = dlmread(['data/' obj.FilePrefix '.marks'], '\t', [1 0 500 15]);
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
        
        amarks = marks;
        amarks.ignore = marks.agresivness == 1;
        amarks.chase = marks.agresivness > 1;
        %
        %[results, match] = SocialBehaviorAnalysisScore(obj, amarks, 'AggressiveChase');
        [results, match] = SocialBehaviorAnalysisScore(obj, marks, 'ChaseEscape');
        
        fprintf('# detection=%.1f%% (%d), false-alarm=%.1f%% (%d)\n', results.chase(1) / (results.chase(1) + results.chase(3)) * 100, results.chase(1), results.chase(2) / results.chase(5) * 100, results.chase(2));
        obj.Contacts.Behaviors.AggressiveChase.Classifier.Labels = match.analysis(6, 1:size(obj.Contacts.Behaviors.AggressiveChase.Classifier.Vec, 2));
        obj.Contacts.Behaviors.AggressiveChase.Classifier.MarkID = match.analysis(8, 1:size(obj.Contacts.Behaviors.AggressiveChase.Classifier.Vec, 2));
    end
%end