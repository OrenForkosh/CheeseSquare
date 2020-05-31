function obj = CheeseBehavior(obj, varargin)
% CheeseBehavior Runs the behavioral analysis

%% Compute group structure (or Hierarchy)
% On the way this also classifies all the contacts between the mice and
% infer Chase-Escape interactions
try
    obj = CheeseHierarchyGroup(obj);
catch err
    fprintf('%s', Console.StrudelError(err));
end

%% Compute behavioral profile for each mouse (and group)
% Computes a table which contains a high-dimensional vector of behaviors
% for each mouse (Such as: No. of approaches, Time in the nest, etc.)
try
    obj = CheeseNewProfiler(obj);
    obj.Save;
catch err
    fprintf('%s', Console.StrudelError(err));
end
