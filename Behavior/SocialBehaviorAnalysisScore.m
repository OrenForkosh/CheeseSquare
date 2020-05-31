function [results, match] = SocialBehaviorAnalysisScore(obj, marks, BehaviorType)
% SocialBehaviorAnalysisScore Auxilary tool for SocialPredPreyModel for
% computing false-alarm and misdetection for labeled chase-escape data.
% Compares the results in 'obj' with labeled 'marks'. Can be used both for
% aggressive and all chase-escape interactions (by setting 'BehaviorType')
results = [];
results.approach = zeros(1, 4);
results.leave = zeros(1, 4);
results.chase = zeros(1, 4);
%%
endf = max([obj.Contacts.List.end]);
begf = min([obj.Contacts.List.beg]);
subjs = [obj.Contacts.List.subjects];
marks.event = {};
output = false;
match.marks = [1:length(marks.beg); marks.agresivness'; marks.beg'; marks.end'; marks.subjects'; ];
b = [obj.Contacts.List.beg]; e = [obj.Contacts.List.end]; s = [obj.Contacts.List.subjects]; 
match.analysis = [1:size(b, 2); b(1, :); e(1, :); s; zeros(4, length(b(1, :)))];
%%
nEvents = sum(begf < max(marks.end));
if ~isfield(marks ,'ignore')
    marks.ignore = false(size(marks.beg));
end
%%
for i=1:length(marks.beg)
    m = ~(marks.beg(i) > endf | marks.end(i) < begf);
    s = (subjs(1, :) == marks.subjects(i, 1)' & subjs(2, :) == marks.subjects(i, 2)') | ...
        (subjs(1, :) == marks.subjects(i, 2)' & subjs(2, :) == marks.subjects(i, 1)');
    m = m & s;
    eidx = find(m);
    marks.event{i} = eidx;
    
    curr = [];
    curr.approach = false(1, 2);
    curr.leave = false(1, 2);
    %curr.pred = false(1, 2);
    %curr.prey = false(1, 2);
    curr.chase = false(1);
    chaseIndex = 0;
    for j=1:length(eidx)
        e = obj.Contacts.List(eidx(j));
        i1 = 2 - (marks.subjects(i, 1) == e.subjects(1));
        i2 = 3 - i1;
        curr.approach([i1, i2]) = curr.approach | [any(e.data{i1} == 2) any(e.data{i2} == 2)];
        
        leave = [any(e.data{i1} == 6) any(e.data{i2} == 6)];
        
        %chase = (any(e.data{i1} == 6) && any(e.data{i2} == 5)) || (any(e.data{i1} == 5) && any(e.data{i2} == 6));
        chase = obj.Contacts.Behaviors.(BehaviorType).Map(eidx(j));
        if chase
            chaseIndex = j;
        end
        
        leave(chase) = false;
        curr.leave([i1, i2]) = curr.leave([i1, i2]) | leave;
        
%         curr.pred([i1, i2]) = curr.pred | [any(e.data{i1} == 5)&&any(e.data{i2} == 6) any(e.data{i2} == 5)&&any(e.data{i1} == 6)];
%         curr.prey([i1, i2]) = curr.prey | prey;

        curr.chase = chase;
    end
    %%
    prevAnalysisIdx = 0;
    matchingIdx = 0;
    if chaseIndex > 0
        prevAnalysisIdx = match.analysis(8, eidx(chaseIndex));
        match.analysis(8, eidx(chaseIndex)) = i;
    end
    if marks.chase(i, :) && curr.chase
        if isempty(eidx)
            match.marks(7, i) = 0;
        else
            match.marks(7, i) = eidx(1);
        end
        match.marks(8, i) = 1;
        chasestr = 'chase  ';
        match.analysis(7, eidx(chaseIndex)) = 1;
        match.analysis(6, eidx(chaseIndex)) = marks.agresivness(i);
    elseif ~marks.chase(i, :) && curr.chase
        if match.analysis(7, eidx) == 0
            match.analysis(6, eidx) = marks.agresivness(i);
            match.analysis(7, eidx) = -1;
        else
            match.analysis(8, eidx(chaseIndex)) = prevAnalysisIdx;
        end
        chasestr = '- chase -';
        if output
            fprintf('%3d (%d-%d) [%3d] %6d %6d %s\n', i-1, marks.subjects(i, 1)', marks.subjects(i, 2)', eidx, marks.beg(i), marks.end(i), chasestr);
        end
    elseif marks.chase(i, :) && ~curr.chase
        if isempty(eidx)
            match.marks(7, i) = 0;
        else
            match.marks(7, i) = eidx(1);
        end
        match.marks(8, i) = -1;
        chasestr = '( chase )';
        if output
            fprintf('%3d (%d-%d) [%3d] %6d %6d %s\n', i-1, marks.subjects(i, 1)', marks.subjects(i, 2)', eidx, marks.beg(i), marks.end(i), chasestr);
        end
    else
        chasestr = 'no-chase';
    end
    
%     curr.predprey = all(curr.pred | curr.prey);
    
    ref = [];
    ref.approach = marks.approach(i, :);
    ref.leave = marks.leave(i, :);
%     ref.predprey = all(marks.pred(i, :) | marks.prey(i, :));
    ref.chase = marks.chase(i, :);

    if ~marks.ignore(i)
        results.approach    = results.approach  + ComputeScore(marks.approach(i, :), curr.approach);
        results.leave       = results.leave     + ComputeScore(marks.leave(i, :), curr.leave);
%     results.predprey    = results.predprey  + ComputeScore(all(marks.pred(i, :) | marks.prey(i, :)), all(curr.pred | curr.prey));
    %     Results.predprey = [Results.predprey; ComputeScore(all(marks.pred(i, :) | marks.prey(i, :)), all(curr.pred | curr.prey))];
        results.chase = results.chase     + ComputeScore(marks.chase(i, :), curr.chase);
    else
    end
   % Results.chase = [Results.chase; ComputeScore(marks.chase(i, :), curr.chase)];
end
results.approach = [results.approach, nEvents];
results.leave = [results.leave, nEvents];
results.chase = [results.chase, nEvents];