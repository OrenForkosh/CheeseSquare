function rank = HierarchicalTopologicalSort(mat)
% HierarchicalTopologicalSort Infers dominance hierarachy from
%   aggressive chases
%
%   [rank, removed, mat] = HierarchicalTopologicalSort(mat)
%   Computes social rank from the number of aggressive chases in matrix
%   'mat' (the rows correspond to the chaser and cols to the "escapers").
%
% The social hierarchy is infered by comparing the number of times 
% mouse A chased B vs. the other way around. If this number was significantly 
% higher for one than the other, then this mouse was ranked higher 
% (more dominant) than the other. Next, the hierarchy map was set as the 
% lowest tree to preserve these ranks. The algorithm is similar to the 
% topological sorting algorithm, but extended in order to allow two or 
% more mice to occupy the same rank. The classification of mice as alpha, 
% beta, gamma, and delta males was done on all four days of experiment. 
% Here we used the Clopper-Pearson test in order to ensure that the 
% relative ranks were significant when constructing the hierarchy maps. The 
% relative ranks were only used when the number of cases, in which one 
% mouse chased another one, was significantly larger than the number of 
% times he escaped from the other mouse to within a 95% percent confidence 
% interval. To quantify the stability of the hierarchy, we measured the 
% fraction of pairs that changed their relative ranking between consecutive 
% days.  

%% if graph is not directed than the algorithm fails 
% (see 'TopoFeedbackArcSetHierarchicalOrder' for a work around)
if ~graphisdag(sparse(mat))
    rank = [];
    return;
end

%%
rank = ones(1, size(mat,1));
r = isValid(mat, rank);
while ~isempty(r)
    rank(r) = rank(r) + 1;
    r = isValid(mat, rank);
end
rank = max(rank) - rank + 1;

function r = isValid(mat, rank)
for i=1:length(rank)
    others = (mat(i, :) > 0) .* rank;
    others(mat(i, :) == 0) = inf;
    if any(others <= rank(i))
        r = find(others <= rank(i), 1);
        return;
    end
end
r = [];

    