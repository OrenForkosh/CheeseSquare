function obj=CheeseFindDarkPeriod(obj, thresh)
fprintf('# finding dark period\n');
obj.Meta.ExperimentTimes = [BinarySearch(obj, thresh, 1) / obj.Video.FrameRate, BinarySearch(obj, thresh, -1) / obj.Video.FrameRate];
fprintf(['# - from ' DateTime.SecToString(obj.Meta.ExperimentTimes(1)) ' to ' DateTime.SecToString(obj.Meta.ExperimentTimes(2)) '\n']);

function frame = BinarySearch(obj, thresh, dir)
m1 = 1;
m2 = obj.Video.NumberOfFrames;

if dir > 0
    valid = m2;
else
    valid = m1;
end
while (m2 >= m1)
    f = round(m1 + (m2 - m1) / 2);
    obj.Video.FrameNumber = f;
    orig = im2double(obj.Video.CurrentFrame);
    lum = max(orig,[],3);
    lum = mean(lum(:));
    if dir * (lum - thresh) > 0
        m1 = f + 1;
    else
        m2 = f - 1;
    end
    if lum <= thresh
        valid = f;
    end
    %if obj.Output
    imagesc(orig);
    axis off;
                
    if lum <= thresh
        title(sprintf(['frame ' num2str(f) ' - Dark (%.2f)'], lum));
    else
        title(sprintf(['frame ' num2str(f) ' - Light (%.2f)'], lum));
    end
    drawnow;
    %end
end
if lum <= thresh
    frame = f;
else
    frame = valid;
end
