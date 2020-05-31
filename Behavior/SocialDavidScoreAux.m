function [w, w2, l, l2] = SocialDavidScoreAux(Pij)
% Compute the auxiliary tables needed for estimating the dominance
% hierarchy based of David's score. See 
%   De Vries, Han, Jeroen MG Stevens, and Hilde Vervaecke. "Measuring and 
%   testing the steepness of dominance hierarchies." Animal Behaviour 71, 
%   no. 3 (2006): 585-592.
w = sum(Pij, 2);
w2 = Pij * w;
l = sum(Pij, 1)';
l2 = Pij' * l;
