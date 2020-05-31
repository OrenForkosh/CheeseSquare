function [v, m] = RobustStd(x, dim)
% Median absolute deviation
% [v, m] = RobustStd(x, dim) computes the MAD along dimension 'dim' and
% returning both the MAD and median
if nargin == 1
    dim = 1;
end
m = median(x, dim);
sz = size(x);
seq = 1:length(sz);
sz(seq(seq ~= dim)) = 1;
v = median(abs(x - repmat(m, sz)), dim);