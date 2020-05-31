function [NormDSonDij, NormDSonPij] = DavidScore(ce)
% Compute Normalized David's score to estimate dominance
% Based upon: Vries, Han de, Jeroen M. G. Stevens, and Hilde 
% Vervaecke. “Measuring and Testing the Steepness of Dominance 
% Hierarchies.” Animal Behaviour 71, no. 3 (March 2006): 585–92. 
% https://doi.org/10.1016/j.anbehav.2005.05.015.
%
%   NormDSonDij = DavidScore(ce) Computes David's score based on chases 
%   matrix ce (the rows represent the chasers, and the columns the escapers) 
%
nij = ce + ce';
Sij = ce;
 
Pij = Sij./nij;
Pij(nij == 0) = 0;
[w, w2, l, l2] = SocialDavidScoreAux(Pij);
DS = w + w2 - l - l2;
NormDSonPij = (DS + size(Sij,1) * (size(Sij,1) - 1)/2) / size(Sij,1);

Dij = (Sij + .5) ./ (nij + 1);
Dij(nij == 0) = 0;
[w, w2, l, l2] = SocialDavidScoreAux(Dij);
DS = w + w2 - l - l2;
NormDSonDij = (DS + size(Sij,1) * (size(Sij,1) - 1)/2) / size(Sij,1);
%plot(sort(NormDS));
