function obj = CheeseHierarchy(obj, force, arg)
% CheeseHierarchy Compute group dominance hierarchy based on chase-escape
% interactions
%   obj = CheeseHierarchy(obj, force, arg) Computes the group structure
%   from the information in 'obj'. If 'force' (optional) is set to true,
%   forces the algorithm to also do the chase-escape analysis (otherwise, 
%   assume it was already done). Set the variable 'arg' if additional args
%   are needed for the SocialPredPreyModel function (see SocialPredPreyModel
%   documantation). 
%
%   Adds the information to 'obj' under the 'Hierarchy' field - see
%   documantation for CheeseSquare for explanation of the fields

obj = CheeseObject.Load(obj);
if ~exist('force', 'var')
    force = true;
end
if ~exist('arg', 'var')
    arg = struct();
end

% if force is true, compute chase-escape interactions, otherwise assume 
% they were previously computed
if force 
    obj = SocialPredPreyModel(obj, arg);
    obj = SocialAggrClassifier(obj);
end

%% Compute hierarchical information 
%  adds the information to 'obj' under the 'Hierarchy' field

% Computes the hiearchy both for all chase-escape (CEs) and just for CEs
% that were classified as aggressive
fields = {'ChaseEscape', 'AggressiveChase'};

for i=1:length(fields)
    edges = cell(1,2); edges{1} = 1:obj.nSubjects; edges{2} = edges{1};
    ce = hist3(...
        [obj.Contacts.Behaviors.(fields{i}).Chaser(obj.Contacts.Behaviors.(fields{i}).Map); ...
        obj.Contacts.Behaviors.(fields{i}).Escaper(obj.Contacts.Behaviors.(fields{i}).Map)]', ...
        'Edges', edges); % number of chase-escapes
    mat = ce - ce';
    mat(binotest(ce, ce + ce')) = 0;
    mat = mat .* (mat > 0);
    
    [rank, removed] = TopoFeedbackArcSetHierarchicalOrder(mat);
    diluted = mat;
    diluted(removed ~= 0) = 0;
    
    obj.Hierarchy.(fields{i}).rank = rank;
    obj.Hierarchy.(fields{i}).nLevels = length(unique(rank));
    obj.Hierarchy.(fields{i}).strength = sum(diluted(:) > 0);
    obj.Hierarchy.(fields{i}).map = diluted;
    obj.Hierarchy.(fields{i}).interactions = obj.Interactions.PredPrey;
    obj.Hierarchy.(fields{i}).removed = removed;
    [q, obj.Hierarchy.(fields{i}).AlphaToDelta] = sort(obj.Hierarchy.ChaseEscape.rank, 'descend');
    obj.Hierarchy.(fields{i}).ChaseEscape = ce;
end
