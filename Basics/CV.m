classdef CV
    % CV Computer vision auxiliary toolbox
    methods (Static = true)
        function lab = RGB2NormLAB(rgb)
            %% convert RGB to CIE 1976 L*a*b* normalized in the range [0, 1]
            lab = rgb2lab(rgb) / 1e4;
            lab(:, :, 2:3) = (1 + lab(:, :, 2:3))/2;
            lab = max(lab, 0);
            lab = min(lab, 1);
        end
        
        function b = Brightness(rgb)
            % compute brightness level of image (converting image to
            % grayscale)
            rgb = im2double(rgb);
            b = max(rgb, [], 3);
        end
        
        function MaskToROI(mask)
            % Converts a mask image (binary image) to a patch and diplay on
            % current image
            %%
            b = bwboundaries(mask);
            Fig.Hon
            for i=1:length(b)
                p = bwtraceboundary(mask, b{i}(1, :), 'N');
                patch(p(:, 2), p(:, 1), 'y', 'EdgeColor', 'w', 'FaceColor', 'y', 'FaceAlpha', .2);
            end
            Fig.Hoff
        end
        
        function img = HistEq(img)
            % histogram equalization for color images
            for i=1:size(img, 3)
                img(:, :, i) = histeq(img(:, :, i));
            end
        end
        
        function img = RegionsToPlot(reg)
            % Plot regions that were extracted using 'regionprops'
            prehold = false;
            if ishold
                prehold = true;
            else
                img = ones([reg.ImageSize 3]);
                imagesc(img);
                Fig.Hon
            end
            cmap = lines(length(reg.PixelIdxList));
            for i=1:length(reg.PixelIdxList)
                [y,x] = ind2sub(reg.ImageSize, reg.PixelIdxList{i});
                %plot(x, y, 'o', 'markerfacecolor', cmap(i, :), 'markeredgecolor', 'none');
                plot(x, y, '.', 'color', cmap(i, :));
                Fig.Hon
            end
            if ~prehold
                Fig.Hoff
            end
            set(gca, 'ydir', 'normal');
        end
        
        function img = RegionsToImage(reg)
            % Convert regions that were extracted using 'regionprops' to
            % image
            img = zeros(reg.ImageSize);
            for i=1:length(reg.PixelIdxList)
                img(reg.PixelIdxList{i}) = i;
            end
            img = ind2rgb(img, [1 1 1; lines]);
        end
        
        function img = RegionsToLabels(reg)
            % Convert regions that were extracted using 'regionprops' to
            % image with different color for each region
            img = zeros(reg.ImageSize);
            for i=1:length(reg.PixelIdxList)
                img(reg.PixelIdxList{i}) = i;
            end
        end
        
        function regions = LabelsToRegions(labels, amplitude)
            % Convert labeled image to ROIs
            nlabels = Q.maxnd(labels);
            regions.ImageSize = size(labels);
            regions.PixelList = cell(1, nlabels);
            regions.PixelIdxList = cell(1, nlabels);
            regions.Amplitude = cell(1, nlabels);
            for i=1:nlabels
                [r,c] = find(labels == i);
                regions.PixelList{i} = [r,c];
                regions.PixelIdxList{i} = sub2ind(regions.ImageSize, r, c);
                if nargin > 1
                    regions.Amplitude{i} = amplitude(regions.PixelIdxList{i});
                end
            end
            ise = cellfun(@isempty, regions.PixelIdxList);
            regions.PixelIdxList(ise) = [];
            regions.PixelList(ise) = [];
            regions.Amplitude(ise) = [];
            regions.NumObjects = length(regions.PixelList);
        end
        
        function [map, score] = Segment(img, background)
            % Simple image segmentation
            opt.nStds = 5;
            opt.MinArea = 550;
            opt.nBins = 10;
            opt.Type = 'classic';
            %%
            switch opt.Type
                case 'color'
                    background = im2double(background);
                    img = im2double(img);
                case {'brighter', 'classic'}
                    background = CV.Brightness(background);
                    img = CV.Brightness(img);
            end
            %%
            nbins = opt.nBins;
            factor = 1/norminv(3/4);
            %thresh = norminv(1-opt.Confidence/2);
            thresh = opt.nStds;
            map = false(size(img, 1), size(img, 2));
            stat.std  = zeros(1, nbins);
            D = img - background;
            for chid=1:size(img, 3)
                bkgch = background(:, :, chid);
                d = D(:, :, chid);
                % NOTE! the following part can (and sould) be computed
                % somewhere else!
                %                 bins = [linspace(0, max(bkgch(:)), nbins) inf];
                %                 [~, ind] = histc(bkgch(:), bins);
                %[bins, ind, count] = Q.binify(bkgch(:), nbins);
                %
                switch opt.Type
                    case {'brighter', 'color'}
                        ind = min(floor(bkgch ./ max(bkgch(:)) * nbins) + 1, nbins);
                        for i=1:nbins
                            stat.std(i) = median(abs(d(ind == i))) * factor;
                            if stat.std(i) == 0
                                stat.std(i) = std(d(ind == i));
                            end
                        end
                        score = d ./ reshape(stat.std(ind), size(d));
                    case 'classic'
                        moments = [mean(d(:)), std(d(:))];
                        score = (d - moments(1)) / moments(2);
                end
                map = map | score > thresh;
            end
            map = bwareaopen(map, opt.MinArea);
        end
        
        function ImageToMice(map)
            % Plot mice as circles
            [label, n] = bwlabel(map, 8);
            imagesc(label);
            Fig.Hon
            for i=1:n
                curr = label == i;
                circles = CV.FitCircles(curr, 2);
                plot(circles(1, 2), circles(1, 1), 'ro')
                plot(circles(2, 2), circles(2, 1), 'g o')
            end
            Fig.Hoff
        end
        
        function FindPath()
            %%
            src = sub2ind(size(im), 1, 1);
            tgt = sub2ind(size(im), 1, size(im, 2));
            %%
            C = im2col(padarray(im, [1 1], inf), [3, 3]);
            L = abs(bsxfun(@minus, C, C((end+1)/2, :)));
            idx = im2col(padarray(reshape(1:numel(im), size(im)), [1 1], 1), [3, 3]);
            %%
            D = inf(size(im, 1), size(im, 2));
            p = nan(size(im, 1), size(im, 2));
            D(src) = 0;
            D = D(:);
            q = true(size(D));
            i = 0;
            while ~isempty(q)
                
                i = i + 1
                [d, u] = min(D);
                if ~isfinite(d)
                    break;
                end
                D(u) = inf;
                alt = d + L(:, u);
                uidx = idx(:, u);
                shorter = alt < D(uidx);
                D(uidx(shorter)) = alt(shorter);
                p(uidx(shorter)) = u;
            end
            
            
            
            %%
        end
        
        function [circles, circmap] = FitCircles(map, n)
            % Fit circles to binary image
            if nargin < 2
                n = 1;
            end
            %%
            circles = zeros(n, 3);
            curr = map;
            distmethod = 'quasi-euclidean';
            for i = 1:n
                d = bwdist(~curr, distmethod);
                [ind, coord] = Q.argmaxnd(d);
                radius = d(ind);
                circles(i,:) = [coord, radius];
                if i < n
                    curr = bwdistgeodesic(curr, coord(2), coord(1), distmethod) > radius;
                end
            end
            %%
            if nargout >= 2 || nargout == 0
                circmap = zeros(size(map), 'uint8');
                for i = 1:n
                    circmap(bwdistgeodesic(map, circles(i, 2), circles(i, 1), distmethod) <= circles(i, 3) & circmap == 0) = i;
                end
            end
            %%
            if nargout == 0
                %%
                imagesc(label2rgb(uint8(map) + circmap, lines));
                Fig.Hon
                for i = 1:n
                    %plot(circles(i, 2), circles(i, 1), 'wo')
                    text(circles(i, 2), circles(i, 1), num2str(i), 'verticalalignment', 'middle', 'horizontalalignment', 'center')
                end
                Fig.Hoff
            end
        end
        
        function c = NumOfNeighbors(a1, sz)
            % Count number of neigbors for binary image (max distance sz)
            if nargin < 2
                sz = 1;
            end
            c = (convn(a1, ones(2*sz+1), 'same') - a1) .* (a1 ~= false);
            
        end
        
        function map = IsNeighbor(a1, a2, conn)
            % Check if two maps are adjacent
            if conn == 4
                zh = false(1, size(a2, 2));
                zv = false(size(a2, 1),1);
                map = (a1 & [zh; a2(1:end-1, :)]) | (a1 & [a2(2:end, :); zh]) | (a1 & [zv, a2(:, 1:end-1)]) | (a1 & [a2(:, 2:end), zv]);
            elseif conn == 8
                [i,j] = ind2sub([3, 3], 1:8);
                map = false(size(a1));
                for idx=Q.exclude(8, 5)
                    map = map | (a1 & Q.shift(a2, [i(idx)-2 j(idx)-2]));
                end
            else
                error
            end
            
        end
        
        
        function lm = NeighMaxima(im, area)
            % maxima in every neighborhood of size 'area'
            sz = (area - 1)/2;
            nim = padarray(im, [sz sz], Q.minval(im));
            lm = col2im(max(im2col(nim, [area area])), [area area], size(nim));
        end
        
        function skel = MedialAxis(im, nscore)
            % Medial axis of a map (Experimental)
            syls = regionprops(im, {'PixelList', 'PixelIdxList'});
            skel = false(size(im));
            for i=1:length(syls)
                r1 = Q.minmax(syls(i).PixelList(:, 1));
                r2 = Q.minmax(syls(i).PixelList(:, 2));
                
                local = zeros(r2(2)-r2(1)+1, r1(2)-r1(1)+1);
                local(sub2ind(size(local), syls(i).PixelList(:, 2)-r2(1)+1, syls(i).PixelList(:, 1)-r1(1)+1)) = nscore(syls(i).PixelIdxList);
                %%
                [mx, idx] = max(local, [], 1);
                local = false(size(local));
                local(sub2ind(size(local), idx, 1:size(local, 2))) = true;
                
                for j=2:size(local, 2)
                    if idx(j) - idx(j-1) > 1
                        mini = nscore(idx(j-1)+1:idx(j)-1, j-1:j);
                        [~, idx2] = max(mini, [], 2);
                        local(sub2ind(size(local), idx(j-1)+1:idx(j)-1, j-2+idx2')) = true;
                    elseif idx(j) - idx(j-1) < -1
                        mini = nscore(idx(j)+1:idx(j-1)-1, j-1:j);
                        [~, idx2] = max(mini, [], 2);
                        local(sub2ind(size(local), idx(j)+1:idx(j-1)-1, j-2+idx2')) = true;
                    end
                end
                %%
                skel(r2(1):r2(2), r1(1):r1(2)) = skel(r2(1):r2(2), r1(1):r1(2)) | local;
%                 imagesc(CV.Join(local>5, skel))
%                 Fig.Hon
%                 plot(1:size(local, 2), idx, 'g.');
%                 
%                 Fig.Hoff
            end
            return;
            
            %%
            skel = bwmorph(im, 'skel', inf);
            d = bwdist(~im);
            while true
            ep = bwmorph(skel, 'endpoints');
            junk = (CV.NeighMaxima(nscore, 3) ~= nscore) & ep;
            if ~any(junk(:))
                break;
            end
            skel(junk) = false;
            end
            
            %%
            %nim = bwmorph(im, 'close');
            nim = im;
            d = bwdist(~nim);
            skel = bwmorph(nim, 'skel', inf);
            ismax = CV.NeighMaxima(d, [3 3]) == d & nim;
%            skel(ismax) = true;
            prev = false(size(skel));
            
            while true
                ep = bwmorph(skel, 'endpoints') | bwmorph(skel, 'branchpoints');
                map = ep & ~prev & ~ismax;
                skel(map) = false;
                if ~any(any(map))
                    break;
                end
                prev = prev | ep;
            end
            d = d .* skel;
        end
        
        function im = Join(varargin)
            im = varargin{1} * 1;
            for i=2:length(varargin)
                im(varargin{i}) = i;
            end
        end
        
        function BlankVideo(filename, nframes, framerate)
            
        end
        
        function SortPoints(x,y)
            %%
            X = x - min(x) + 1;
            Y = y - min(y) + 1;
            m = false(max(Y), max(X));
            m(sub2ind(size(m), Y, X)) = true;
            %%
            d = inf(size(m));
            pos = [Y(1), X(1)];
            d(pos(1), pos(2)) = 0;
            
            imagesc(CV.Join(m, d<inf))
        end
    end
end
