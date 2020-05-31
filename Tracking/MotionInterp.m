function [ox, oy] = MotionInterp(x, y, range, ~, ~)
% Interpolation of motion paths in 'x' and 'y'
ox = interp1([range(1)-1 range(end)+1], x([range(1)-1 range(end)+1]), range);
oy = interp1([range(1)-1 range(end)+1], y([range(1)-1 range(end)+1]), range);

