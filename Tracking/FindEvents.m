function [beg, fin, len, events] = FindEvents(data, map)
% Find segments of true valued indices in the 'map' array
%
%       Created by OREN FORKOSH
%

if ~exist('map', 'var')
    map = data;
end
% c = conv([0 map 0], [1 -1]);
% beg = find(c > 0) - 1;
% fin = find(c < 0) - 2;
c = diff([0 map 0]);
beg = find(c > 0);
fin = find(c < 0) - 1;

len = fin - beg + 1;
events = cell(1, length(beg));
if nargout > 3
    for i=1:length(beg)
        events{i}.data = data(beg(i):fin(i));
        events{i}.BegFrame = beg(i);
        events{i}.EndFrame = fin(i);
    end
end