function [rank, removed, mat] = TopoFeedbackArcSetHierarchicalOrder(mat)
% TopoFeedbackArcSetHierarchicalOrder Infers dominance hierarachy from
%   aggressive chases. Like 'HierarchicalTopologicalSort' but also copes
%   with undirected graphs.
%
%   [rank, removed, mat] = TopoFeedbackArcSetHierarchicalOrder(mat)
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

%% if this is a directed graph than use 'HierarchicalTopologicalSort'
if graphisdag(sparse(mat))
    removed = mat * 0;
    rank = HierarchicalTopologicalSort(mat);
    return;
end
%% Undirected graphs
%  In this case, look exhustively at all possible choices of a subset of
%  the interactions choosing the permutation that preserves most the the
%  interactions but is directed.
idx = find(mat);
perm = nchoose(idx);
count = zeros(1, length(perm));
for i=1:length(perm)
    s = sum(mat(perm{i})) / sum(mat(:));
    count(i) = length(perm{i});
end
[count, order] = sort(count);
newperm = perm;
for i=1:length(perm)
    newperm{i} = perm{order(i)};
end
perm = newperm;
%
cost = zeros(1, length(perm));
for i=1:length(perm)
    cost(i) = sum(mat(perm{i}));
end
[cost, order] = sort(cost(:), 1, 'ascend');
newperm = perm;
for i=1:length(perm)
    newperm{i} = perm{order(i)};
end
perm = newperm;
%
for i=1:length(perm)
    c = mat;
    c(perm{i}) = 0;
    if graphisdag(sparse(c))
        mat = c;
        rank = myrank(c);
        removed = mat * 0;
        removed(perm{i}) = 1;
        break;
    end
end
