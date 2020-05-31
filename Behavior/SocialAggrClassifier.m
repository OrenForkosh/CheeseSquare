function [obj, Score] = SocialAggrClassifier(obj, models)
% Classify Chase-Escape interactions as aggressive (or non-aggressive).
% Usese three different classifiers that were trained on labeled data (see
% documentation)
%
%   obj = SocialAggrClassifier(obj) Classify Chases as aggressive or not
%   and puts the outcome in 'obj.Contacts.Behaviors.AggressiveChase'
%
%   obj = SocialAggrClassifier(obj, models) 'models' is the classifiers
%   that the algorithm will use. By default, the models are loaded from the
%   ./Prototypes/AggrPredPrey.mat file
%
%       Created by OREN FORKOSH
%

if ischar(obj)
    obj = TrackLoad(obj);
    obj = SocialPredPreyModel(obj);
end

if nargin < 2
    [pathstr, name, ext] = fileparts(which(mfilename));
%    models = {'Res/SC.exp0001.day02.cam01.obj.mat'};
    pathstr = regexprep(pathstr, '\\', '/');
    models = {[pathstr '/Prototypes/AggrPredPrey.mat']};
end

if ~iscell(models)
    models = {models};
end
obj.Contacts.Behaviors.AggressiveChase.Map = zeros(1, size(obj.Contacts.Behaviors.AggressiveChase.Classifier.Vec, 2));

nObservations = 0;
for i=1:length(models)
    currModel = models{i};
    if ischar(currModel)
        currModel = TrackLoad(currModel, {'Contacts'});
    end
    
    Classifier = currModel.Contacts.Behaviors.AggressiveChase.Classifier;
    vec = obj.Contacts.Behaviors.AggressiveChase.Classifier.Vec';
    
    train = Classifier.Train.Data;
    label = Classifier.Train.Labels;
    
    for i=1:size(vec, 2)
        vec(~isfinite(vec(:, i)), i) = min(isfinite(vec(:, i)));
    end
    
%     tree = cellfun(@str2num, eval(Classifier.Tree, vec));
    classifierTree = fitctree(train, label);
    tree = predict(classifierTree, vec);
    da = classify(vec(:, Classifier.FeatureSubset), train(:, Classifier.FeatureSubset), label, Classifier.DA.Method, [1.5 ones(1, max(label))*0.1/max(label)]);
    % knn = knnclassify(vec(:, Classifier.FeatureSubset), train(:, Classifier.FeatureSubset), label, Classifier.KNN.nNeighbours, Classifier.KNN.Method);
    knnmdl = fitcknn(train(:, Classifier.FeatureSubset), label, 'NumNeighbors', Classifier.KNN.nNeighbours, 'Distance', Classifier.KNN.Method);
    knn = knnmdl.predict(vec(:, Classifier.FeatureSubset));
    
    obj.Contacts.Behaviors.AggressiveChase.Map = obj.Contacts.Behaviors.AggressiveChase.Map + ...
        ((tree' > 0) + (da'>0) + (knn'>0)) * sum(Classifier.Train.Labels);
    nObservations = nObservations + sum(Classifier.Train.Labels);
end
obj.Contacts.Behaviors.AggressiveChase.Map = (obj.Contacts.Behaviors.AggressiveChase.Map / nObservations) >= 1;
%%
obj.Contacts.Behaviors.AggressiveChase.Chaser = obj.Contacts.Behaviors.AggressiveChase.Map * 0;
obj.Contacts.Behaviors.AggressiveChase.Chaser(obj.Contacts.Behaviors.AggressiveChase.Map) = obj.Contacts.Behaviors.ChaseEscape.Chaser(obj.Contacts.Behaviors.AggressiveChase.Map);
obj.Contacts.Behaviors.AggressiveChase.Escaper = obj.Contacts.Behaviors.AggressiveChase.Map * 0;
obj.Contacts.Behaviors.AggressiveChase.Escaper(obj.Contacts.Behaviors.AggressiveChase.Map) = obj.Contacts.Behaviors.ChaseEscape.Escaper(obj.Contacts.Behaviors.AggressiveChase.Map);

%%
if obj.OutputToFile
    fprintf('# - saving data\n');
    %TrackSave(obj);
end

%%
found = true;
try
    marks = SocialBehaviorLoad(obj);
catch me
    found = false;
end
if found
    %%
    marks.aggressiveChase = marks.agresivness > 1;
    Score = SocialBehaviorComputeScore(obj, 'aggressiveChase', marks, ...
        obj.Contacts.Behaviors.AggressiveChase.Map, marks.agresivness == 1 | marks.artifact == 999 | marks.artifact == 998) % 998 means missed contact
    fprintf('#----------------\n');
    fprintf('# detection=%.1f%% (%d/%d), false-alarm=%.1f%% (%d)\n', Score.nDetection/Score.nTrue*100, Score.nDetection, Score.nTrue, Score.nFalseAlarm / Score.nEvents * 100, Score.nFalseAlarm);
    %     marks.ignore = marks.agresivness == 1;
    %     marks.chase = marks.agresivness > 1;
    %     %
    %     [results, match] = SocialBehaviorAnalysisScore(obj, marks, 'AggressiveChase');
    %     fprintf('#----------------\n');
    %     fprintf('# detection=%.1f%% (%d/%d), false-alarm=%.1f%% (%d)\n', results.chase(1) / (results.chase(1) + results.chase(3)) * 100, results.chase(1), results.chase(1) + results.chase(3), results.chase(2) / results.chase(5) * 100, results.chase(2));
end
