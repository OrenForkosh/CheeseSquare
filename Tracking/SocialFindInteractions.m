function [obj, local] = SocialFindInteractions(obj)
% SocialFindInteractions Find contacts and extract interaction properties.
% These properties are then used to classify the different types of
% interactions (see SocialPredPreyModel)
%
%       Created by OREN FORKOSH
%

local = obj;

%% find hidden and sheltered events
if local.Interactions.MinShelteredDuration > 0
   for s=1:local.nSubjects
       [b, e, len] = FindEvents(local.sheltered(s, :));

       idx = find(len < local.Interactions.MinShelteredDuration & b > find(obj.valid, 1) & e < find(obj.valid, 1, 'last'));
       for i=idx
           seq = linspace(local.x(s, b(i) - 1), local.x(s, e(i) + 1), len(i) + 2); local.x(s, b(i):e(i)) = seq(2:end-1);
           seq = linspace(local.y(s, b(i) - 1), local.y(s, e(i) + 1), len(i) + 2); local.y(s, b(i):e(i)) = seq(2:end-1);
           local.zones(s, b(i):e(i)) = 1;
           local.sheltered(s, b(i):e(i)) = false;
       end

       [b, e, len] = FindEvents(~local.hidden(s, :));
       idx = find(len < local.Interactions.MinEventDuration & b > find(obj.valid, 1) & e < find(obj.valid, 1, 'last'));
       for i=idx
           local.hidden(s, b(i):e(i)) = true;
       end
       
       [b, e, len] = FindEvents(local.hidden(s, :));
       idx = find(len < local.Interactions.MinEventDuration & b > find(obj.valid, 1) & e < find(obj.valid, 1, 'last'));
       for i=idx
           local.hidden(s, b(i):e(i)) = false;
       end
       
       e = find(local.hidden(s, 1:end-1) & local.sheltered(s, 2:end) & ~local.sheltered(s, 1:end-1));
       b = e * 0;
       for i=1:length(e)
           pos = find(~local.hidden(s, 1:e(i)), 1, 'last');
           if ~isempty(pos)
               b(i) = pos;
               local.sheltered(s, b(i):e(i)) = true;
           end
       end
   end
end
local.Interactions.valid = ~local.sheltered;

%% convert trajectories from pixels to centimeters
local = TrackSmooth(local, local.Interactions.SmoothSpan);

local.xCM = local.x / local.PixelsPerCM;
local.yCM = local.y / local.PixelsPerCM;

fprintf('# computing mutual quantities\n');

%% find contacts
fprintf('# - finding contacts\n');
local.Interactions.distanceCM = cell(local.nSubjects,local.nSubjects);
local.Interactions.contact  = false(local.nSubjects,local.nSubjects, local.nFrames);
local.Interactions.proximity  = false(local.nSubjects,local.nSubjects, local.nFrames);

for i=1:local.nSubjects-1
    for j=i+1:local.nSubjects
        local.Interactions.distanceCM{i, j} = sqrt((local.xCM(i, :) - local.xCM(j, :)).^2 + (local.yCM(i, :) - local.yCM(j, :)).^2);
        sheltered = local.sheltered(i, :) | local.sheltered(j, :);
        contact = local.Interactions.distanceCM{i, j} < local.Interactions.ZoneOfContactCM & ~sheltered;
        proximity = local.Interactions.distanceCM{i, j} < local.Interactions.ZoneOfProximityCM & ~sheltered;

        % remove short contacts
        [start, finish, len] = FindEvents(contact);
        for r=find(len < local.Interactions.MinContanctDuration)
            contact(start(r):finish(r)) = false;
        end
        
        [start, finish] = FindEvents(contact);
        for r=1:length(start)-1
            if min(sheltered(finish(r) + 1:start(r+1) - 1)) > 0 && start(r+1) - finish(r) - 1 < local.Interactions.MaxShelteredContanctDuration
                contact(finish(r) + 1:start(r+1) - 1) = true;
            end
        end
        
        
        local.Interactions.contact(i, j, :) = contact;
        local.Interactions.proximity(i, j, :) = proximity;
        local.Interactions.contact(j, i, :) = contact;
        local.Interactions.proximity(j, i, :) = proximity;
    end
end

%% relative angles between subjects
fprintf('# - computing relative angles\n');
local.Interactions.angle = cell(local.nSubjects, local.nSubjects);
local.Interactions.jump = cell(local.nSubjects, local.nSubjects);
local.speedCMperSec = zeros(local.nSubjects, local.nFrames);
obj.speed = zeros(obj.nSubjects, obj.nFrames);
for i=1:local.nSubjects
    obj.speed(i, :) = [sqrt((obj.x(i, 2:end) - obj.x(i, 1:end-1)).^2 + (obj.y(i, 2:end) - obj.y(i, 1:end-1)).^2) 0] / obj.dt;
    
    local.speedCMperSec(i, :) = [sqrt((local.xCM(i, 2:end) - local.xCM(i, 1:end-1)).^2 + (local.yCM(i, 2:end) - local.yCM(i, 1:end-1)).^2) 0] / local.dt;
    for j=1:local.nSubjects
        if i == j;
            continue;
        end
        relative = [local.x(j, :) - local.x(i, :); local.y(j, :) - local.y(i, :)];
        relative = relative ./ repmat(sqrt(sum(relative.^2)), 2, 1);
        
        speed = [local.x(i, 2:end) - local.x(i, 1:end-1) 0; local.y(i, 2:end) - local.y(i, 1:end-1) 0];
        direction = speed ./ repmat(sqrt(sum(speed.^2)), 2, 1);
        
        c = sum(direction .* relative);
        c(c >  1) =  1;
        c(c < -1) = -1;
        angle = acos(c);
        local.Interactions.angle{i, j} = angle;
        clear speed;
    end
end

%%
obj.Interactions = local.Interactions;
obj.speedCMperSec = local.speedCMperSec;
