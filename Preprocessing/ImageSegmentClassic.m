function [map, img] = ImageSegmentClassic(obj, frame, useropt)
%%
opt.VideoScale = 1/4;
opt.UseAdaptiveThresh = true;
opt.NoiseThresh = 10;
opt.ValidROI = [];
opt.MaxNumOfObjects = 320/(1/opt.VideoScale)^2;
opt.MinNumOfPixels = 640/(1/opt.VideoScale)^2;
opt.ErodeRadios = 8/(1/opt.VideoScale);

%%
if ischar(obj)
    load(obj);
end

if isscalar(frame)
    obj.Video.FrameNumber = frame;
    img = obj.Video.CurrentFrame;
else
    img = frame;
end
%%
try
    m = imsubtract(img, obj.Background.im);
    m = imresize(m, opt.VideoScale);
    %%
    hsv = rgb2hsv(m);
    vm = hsv(:,:,3);
    %%
    meanBKG = mean(vm(:));
    stdBKG = std(vm(:));
    if opt.UseAdaptiveThresh
        upper = opt.NoiseThresh;
        lower = 1;
        prev_thresh = round((upper + lower)/2);
        while true
            thresh = round((upper + lower)/2);
            bw = vm > meanBKG + thresh * stdBKG;
            if ~isempty(opt.ValidROI)
                bw(~opt.ValidROI) = false;
            end
            cc = bwconncomp(bw);
            if cc.NumObjects < opt.MaxNumOfObjects
                upper = thresh - 1;
                prev_thresh = thresh;
            else
                lower = thresh + 1;
            end
            if lower > upper
                break
            end
        end
        if thresh ~= prev_thresh
            thresh = prev_thresh;
            bw = vm > meanBKG + thresh * stdBKG;
        end
    else
        thresh = obj.NoiseThresh;
        bw = vm > meanBKG + thresh * stdBKG;
        cc = bwconncomp(bw);
        if cc.NumObjects > opt.MaxNumOfObjects
            return;
        end
    end
    boundries = bwareaopen(bw, opt.MinNumOfPixels);
    map = imerode(boundries, strel('disk', opt.ErodeRadios));
    map = imresize(map, 1/opt.VideoScale);
catch
    img = [];
    map = [];
    if isscalar(frame)
        fprintf('# - unable to segmente video frame no. %d', frame);
    else
        fprintf('# - unable to segmente video frame');
    end
end