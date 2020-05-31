function Score = SocialBehaviorComputeScore(obj, behavior, marks, results, ignore)
% SocialAggrClassifier Auxilary tool for SocialAggrClassifier for
% computing false-alarm and misdetection for labeled aggressive chase-escape data.
% Compares the results in 'obj' with labeled 'marks'

if nargin < 5
    ignore = false(1, length(marks.beg));
end

Contacts.beg = min([obj.Contacts.List.beg]);
Contacts.end = max([obj.Contacts.List.end]);
Contacts.subjects = sort([obj.Contacts.List.subjects]);
Contacts.results = results;
Contacts.results(length(Contacts.results) + 1:length(Contacts.beg)) = false;

Score.Results = struct(...
    'beg', num2cell(Contacts.beg),...
    'end', num2cell(Contacts.end),...
    'output', num2cell(Contacts.results),...
    'subjects', mat2cell(Contacts.subjects, 2, ones(1, length(Contacts.beg))),...
    'match', false,...
    'markIdx', 0, ...
    'ignore', false);
Score.Marks = struct(...
    'beg', num2cell(marks.beg),...
    'end', num2cell(marks.end),...
    'subjects', mat2cell(marks.subjects', 2, ones(1, length(marks.beg)))',...
    'output', num2cell(marks.(behavior)),...
    'match', false,...
    'resultIdx', [], ...
    'ignore', num2cell(ignore));

idx = 1;
for i=1:length(marks.beg)
    %%
    CurrMark.beg = marks.beg(i);
    CurrMark.end = marks.end(i);
    CurrMark.subjects = sort(marks.subjects(i, :))';
    CurrMark.output = marks.(behavior)(i);
    %% find overlap between marks and results
    % overlap times
    overlap = ...
        (CurrMark.end >= Contacts.beg & CurrMark.end <= Contacts.end) | ...
        (CurrMark.end >= Contacts.end & CurrMark.beg <= Contacts.end);
    % overlap subjects
    overlap = overlap & ...
        CurrMark.subjects(1, :) == Contacts.subjects(1, :) & ...
        CurrMark.subjects(2, :) == Contacts.subjects(2, :);
    %% register matches
    for ridx=find(overlap)
        if length(Score.Results) >= ridx && ~isempty(Score.Results(ridx).match) && Score.Results(ridx).match
            continue;
        end
        Score.Results(ridx).match = CurrMark.output == results(ridx);
        Score.Results(ridx).markIdx = i;
        Score.Results(ridx).ignore = ignore(i);
        
        Score.Marks(i).match = Score.Marks(i).match || Score.Results(ridx).match;
    end
    Score.Marks(i).resultIdx = find(overlap);
end
Score.Results = Score.Results([Score.Results.beg] <= max([Score.Marks.end]));
%% compute final score
Score.nFalseAlarm = sum([Score.Results.output] & ~[Score.Results.match] & ~[Score.Results.ignore]);
Score.nDetection = sum([Score.Marks.output] & [Score.Marks.match] & ~[Score.Marks.ignore]);
Score.nTrue = sum(marks.(behavior));
Score.nEvents = length(marks.(behavior));
Score.resultsDetectionIdx = find([Score.Results.output] & [Score.Results.match] & ~[Score.Results.ignore]);
Score.resultsFalseAlarmIdx = find([Score.Results.output] & ~[Score.Results.match] & ~[Score.Results.ignore]);
