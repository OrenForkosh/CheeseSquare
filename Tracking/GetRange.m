function r = GetRange(map)
% Find start of first segment of true valued indices in the 'map' array, as
% well as the end of the last segment
r1 = find(map == 1, 1);
if isempty(r1)
    r = [];
    return;
end
r2 = find(map == 1, 1,'last');
r = [r1 r2];