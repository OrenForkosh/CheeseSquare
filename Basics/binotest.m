function [issig, pvalue] = binotest(X, N, hypP, alpha)
% BINOTEST Runs a binomial tests
%   [issig, pvalue] = binotest(X, N) runs binomial test for X
%   success out of N trials. Returns whether the trial is significant and
%   the p-value.
%
%   binotest(..., hypP) runs binomial test with a given null
%   hypothesis (the default is 50% success rate or hypP=.5)
%
%   binotest(..., hypP, alpha) runs binomial test with a 'hypP' null
%   hypothesis and 'alpha' confidence interval

if ~exist('hypP', 'var')
    hypP = .5;
end
if ~exist('alpha', 'var')
    alpha = 0.05;
end

side = N * hypP < X;
c = binocdf(X-side, N, hypP);
pvalue = side - side .* c + (1-side) .* c;
issig = pvalue > alpha;
