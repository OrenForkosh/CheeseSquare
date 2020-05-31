classdef HandT_tools
    
    methods (Static = true)
        
        %% features vectors extracting
        function [abMat, featsVecsMat, ids] = extractFeaturesVectors(flabels, nSubjects, vidScale, arenaWidth, arenaHeight, arenaCoord, frame, background, contactDistTresh)
            % given 'flabels' (the matrix computed by CheeseColorSegment
            % where each pixel is given an integer value corresponding to
            % the ID of the mouse it belongs to or 0 otherwise) and the
            % other parameters which are pretty self explanatory, returns
            % 'abMat': the matrix of the bodyEnds coordinates of the blobs,
            % 'featsVecsMat': the features vectors corresponding to 'abMat'
            % and 'ids': the IDs corresponding
            featsVecsMat = -ones(nSubjects, 7626);
            abMat = -ones(nSubjects, 4);
            ids = (1:nSubjects+1)';
            subImgSizeBig = 20;
            bodyOrient = zeros(nSubjects); bodyOrientRadian = zeros(nSubjects);
            bodyLength = zeros(nSubjects);
            end_x = -ones(nSubjects, 2); end_y = end_x;
            bodyEnd_A = end_x; bodyEnd_B = end_x;
            getRidOfNonContactFromWithin = [];
            
            if contactDistTresh ~= -1 % we want contact frames only
                xScale = arenaWidth / arenaCoord(3);
                yScale = arenaHeight / arenaCoord(4);
            end
            
            %             vidScale = objClassic.VideoScale;
            scaleX = arenaCoord(3) / 568.8320;  %currently normalizing to 1st group (1st day) scale, on which we determined the thresholds
            scaleY = arenaCoord(4) / 464.9143;  %currently normalizing to 1st group scale, on which we determined the thresholds
            scaleX = scaleX / 70 * arenaWidth;  % obj.Meta.Scale.ArenaWidth;
            scaleY = scaleY / 50 * arenaHeight; % obj.Meta.Scale.ArenaHeight;
            scaleMean = mean([scaleX, scaleY]);
            
            labelStats = regionprops(flabels, 'all');
            labelStats = labelStats([labelStats.Area]>1);
            foreGrndBIG = imsubtract(frame, background);
            foreGrndSmall = imresize(foreGrndBIG, vidScale);
            
            miceInFrame = unique(flabels(:)); % getting the IDs of the mice present in the frame
            missing = 0;
            
            
            for k = 1 : nSubjects + 1
                
                if ~ismember(k, miceInFrame)
                    missing = missing + 1;
                    continue
                end
                bodyIdx = k - missing;
                
                [bodyStatsAdded] = HandT_tools.analyzeBodyStats2(labelStats(bodyIdx), vidScale, false);
                
                bodyOrient(bodyIdx)= bodyStatsAdded.Orientation;
                bodyLength(bodyIdx)= bodyStatsAdded.MajorAxisLength;
                bodyOrientRadian(bodyIdx) = deg2rad(bodyStatsAdded.Orientation);
                [end_x(bodyIdx),end_y(bodyIdx)] = pol2cart(bodyOrientRadian(bodyIdx), bodyLength(bodyIdx) / 2);
                bodyEnd_A(bodyIdx, :) = bodyStatsAdded.Centroid + [end_x(bodyIdx), -end_y(bodyIdx)];
                bodyEnd_B(bodyIdx, :) = bodyStatsAdded.Centroid - [end_x(bodyIdx), -end_y(bodyIdx)];
                
                centRndSmall = round(labelStats(bodyIdx).Centroid);                
                
                mouseSize = round(20 * scaleMean);
                % here we're getting the region which surrounds our actual blob's centroid
                rowIdx = max(1, centRndSmall(2) - mouseSize)     : min(size(foreGrndSmall, 1), centRndSmall(2) + mouseSize);
                colIdx = max(1, centRndSmall(1) - mouseSize - 4) : min(size(foreGrndSmall, 2), centRndSmall(1) + mouseSize + 4);
                mouseFlabels = flabels(rowIdx, colIdx);
                mouseFlabBIG = imresize(mouseFlabels, 1/vidScale, 'nearest');
                flabRotatedBIG = imrotate(mouseFlabBIG, -bodyOrient(bodyIdx), 'loose'); % rotating the image such that the major axis of the ellipse will be parallel to the x-axis
                for i = 1:size(labelStats,1)
                    labelStats(i).label = flabels(labelStats(i).PixelIdxList(1));
                end
                flabRotStatsBIG = regionprops(flabRotatedBIG == labelStats(bodyIdx).label, 'Centroid','BoundingBox','Area','PixelList');
                [~,I] = max([flabRotStatsBIG.Area]); % sometimes the tail gets split into small blobs...
                flabRotStatsBIG = flabRotStatsBIG(I);% only keep the biggest (body) blob stats
                leftEndBig = flabRotStatsBIG.PixelList(1, :);
                rightEndBig = flabRotStatsBIG.PixelList(end, :);
                
                centBig = round(labelStats(bodyIdx).Centroid / vidScale);
                mouseSizeBIG = round(mouseSize/vidScale);
                rowIdxBIG = max(1, centBig(2) - mouseSizeBIG)      : min(size(foreGrndBIG, 1), centBig(2) + mouseSizeBIG);
                colIdxBIG = max(1, centBig(1) - mouseSizeBIG - 16) : min(size(foreGrndBIG, 2), centBig(1) + mouseSizeBIG + 16);
                foreGrndBIG = imsubtract(frame, background);                                                                       % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                mouseFgrndBIG = foreGrndBIG(rowIdxBIG, colIdxBIG, :);
                imgRotatedBIG = imrotate(mouseFgrndBIG, -bodyOrient(bodyIdx), 'loose');
                
                % here we're getting the regions which surround the leftest
                % blob's pixel and the rightest one: our candidates to be
                % head and tail
                [leftSide, rightSide] = HandT_tools.centerEndsConcatZeros2(imgRotatedBIG, leftEndBig, rightEndBig, subImgSizeBig);
                rightSide = imrotate(rightSide, 180);
                imgCutBig = [leftSide, rightSide]; % this is the feature vector as a 2D image (before it's reshaped into a single row vector)
                
                a = bodyEnd_A(bodyIdx, :);
                b = bodyEnd_B(bodyIdx, :);
                if a(1) < b(1)  % we want 'a' to be the rightest extrema
                    tmp = a;
                    a = b;
                    b = tmp;
                end
                if a(1) == b(1)
                    if a(2) > b(2) % since the rotation is counterclockwise (and given the y-axis orientation)
                        tmp = a;
                        a = b;
                        b = tmp;
                    end
                end
                abMat(k, :) = [b, a] .* (1/vidScale); % franck scalingPb
%                 abMat(k, :) = [b, a] .* 4; % abMat(bodyIdx, :) = [b, a] .* 4; % franck scalingPb
                featsVecsMat(k, :) = imgCutBig(:)'; % featsVecsMat(bodyIdx, :) = imgCutBig(:)';
            end
            
            absents = find(featsVecsMat(:,1)==-1);
            featsVecsMat(absents, :) = [];
            abMat(absents, :) = [];
            ids(absents) = [];
            
            for s = getRidOfNonContactFromWithin
                if contactDistTresh ~= -1
                    inContact = [];
                    for i = 1 : size(abMat,1)
                        if ismember(i, inContact)
                            continue
                        end
                        for j = 1 : size(abMat, 1)
                            if i == j
                                continue
                            end
                            actDist = sqrt(((abMat(i,1) - abMat(j,1)) * xScale)^2 + ((abMat(i,2) - abMat(j,2)) * yScale)^2);
                            if actDist <=  contactDistThresh
                                inContact = [inContact, i, j];
                                break
                            end
                        end
                    end
                    inContact = sort(unique(inContact));
                    featsVecsMat = featsVecsMat(inContact);
                    abMat = abMat(inContact);
                    ids = ids(inContact);
                end
            end
            
        end
        
        function [subImgA, subImgB] = centerEndsConcatZeros2(imgRotated, leftEndSmall, rightEndSmall, subImgSizeSmall)
            % given the image 'imgRotated' and the coordinates of two
            % points (our head and tail candidates), returns 2 sub-images
            % around (the "radius" is given by 'subImgSizeSmall') our 2
            % points then paddle them with zeros in case they don't have
            % the wanted size
            rowIdxLeft =    max(1, leftEndSmall(1,2) - subImgSizeSmall)             : min(size(imgRotated, 1), leftEndSmall(1,2) + subImgSizeSmall);
            colIdxLeft =    max(1, leftEndSmall(1,1) - subImgSizeSmall)             : min(size(imgRotated, 2), leftEndSmall(1,1) + round(subImgSizeSmall/2));
            rowIdxRight =   max(1, rightEndSmall(1,2) - subImgSizeSmall)            : min(size(imgRotated, 1), rightEndSmall(1,2) + subImgSizeSmall);
            colIdxRight =   max(1, rightEndSmall(1,1) - round(subImgSizeSmall/2))   : min(size(imgRotated, 2), rightEndSmall(1,1) + subImgSizeSmall);
            
            %% get the sub image
            % TODO: make bigger image for whole mouse
            subImgA = imgRotated(rowIdxLeft, colIdxLeft, :);
            subImgB = imgRotated(rowIdxRight, colIdxRight, :);
            %% concatinate zeros to images cut by the borders of the image
            %% TODO: allocate a zeros(size(image)) mat and put the (maybe partial) image in it...
            if 1 > leftEndSmall(1,1) - subImgSizeSmall
                subImgA = [zeros(size(subImgA,1), -leftEndSmall(1,1) + subImgSizeSmall + 1, 3), subImgA];
            end
            if size(imgRotated, 2) < leftEndSmall(1,1) + subImgSizeSmall
                subImgA = [subImgA, zeros(size(subImgA, 1), leftEndSmall(1,1) + subImgSizeSmall - size(imgRotated, 2), 3)];
            end
            if 1 > leftEndSmall(1,2) - subImgSizeSmall
                subImgA = [zeros(-leftEndSmall(1,2) + subImgSizeSmall + 1, size(subImgA,2), 3); subImgA];
            end
            if size(imgRotated, 1) < leftEndSmall(1,2) + subImgSizeSmall
                subImgA = [subImgA; zeros((leftEndSmall(1,2) + subImgSizeSmall - size(imgRotated, 1)), size(subImgA,2), 3)];
            end
            %%
            if 1 > rightEndSmall(1,1) - subImgSizeSmall
                subImgB = [zeros(size(subImgB,1), -rightEndSmall(1,1) + subImgSizeSmall + 1, 3), subImgB];
            end
            if size(imgRotated, 2) < rightEndSmall(1,1) + subImgSizeSmall
                subImgB = [subImgB, zeros(size(subImgB, 1), rightEndSmall(1,1) + subImgSizeSmall - size(imgRotated, 2), 3)];
            end
            if 1 > rightEndSmall(1,2) - subImgSizeSmall
                subImgB = [zeros(-rightEndSmall(1,2) + subImgSizeSmall + 1, size(subImgB,2), 3); subImgB];
            end
            if size(imgRotated, 1) < rightEndSmall(1,2) + subImgSizeSmall
                subImgB = [subImgB; zeros((rightEndSmall(1,2) + subImgSizeSmall - size(imgRotated, 1)), size(subImgB,2), 3)];
            end
            % end
        end
        
        function [bodyStats] = analyzeBodyStats2(bodyStats, scale, toPlot)
            % add a few properties to those returned by the 'regionprops'
            % function:
            % -the most distant points among the blob's pixels
            % -the angle between the two previous points axis and the
            % x-axis
            % the others properties added are not used anywhere and will be
            % deleted in a future release
            [bodyEndsA, bodyEndsB] = HandT_tools.realFarthestPoints2(bodyStats);
            bodyEndsA = bodyEndsA / scale;
            bodyEndsB = bodyEndsB / scale;
            
            diffs = bodyEndsA - bodyEndsB;
            mids = bodyEndsB + diffs./2;
            bodyAngles = atand(diffs(:,2) ./ diffs(:,1));
            
            % cents = reshape([bodyStats.Centroid],[],2)';
            cents = [bodyStats.Centroid] / scale;
            centroidVector = reshape(cents,[],2) - bodyEndsB;
            midleVector = mids - bodyEndsB;
            lengthCentroid = norm(centroidVector);
            lengthMidle = norm(midleVector);
            for i = 1 : size(centroidVector,1)
                x1 = centroidVector(i,1);
                y1 = centroidVector(i,2);
                x2 = midleVector(i,1);
                y2 = midleVector(i,2);
                centroidMidleAngle = atan2d(x1 * y2 - y1 * x2 , x1 * x2 + y1 * y2);
                
                if toPlot
                    hold on
                    plot(cents(i,1), cents(i,2), 'o')
                    plot(mids(i,1), mids(i,2), 'o')
                    plot([bodyEndsA(i,1),bodyEndsB(i,1)],[bodyEndsA(i,2),bodyEndsB(i,2)],'LineWidth',2)
                    if lengthCentroid > lengthMidle
                        plot(bodyEndsA(i,1),bodyEndsA(i,2),'*w','MarkerSize', 20, 'LineWidth', 1)
                        text(bodyEndsA(i,1)+10,bodyEndsA(i,2)+10, num2str(lengthCentroid / lengthMidle), 'Color', 'w', 'FontSize', 14)
                    else
                        plot(bodyEndsB(i,1),bodyEndsB(i,2),'*w','MarkerSize', 20, 'LineWidth', 1)
                        text(bodyEndsB(i,1)+10,bodyEndsB(i,2)+10, num2str(lengthCentroid / lengthMidle), 'Color', 'w', 'FontSize', 14)
                    end
                end
                bodyStats.centroidMidleAngle = centroidMidleAngle;
                bodyStats.lengthCentroid = lengthCentroid;
                bodyStats.lengthMidle = lengthMidle;
                bodyStats.bodyAngles = bodyAngles;
                bodyStats.bodyEndsA = bodyEndsA;
                bodyStats.bodyEndsB = bodyEndsB;
            end
        end
        
        function [bodyEndsA, bodyEndsB] = realFarthestPoints2(bodyStats)
            % given a list of blobs, returns the most distant two pixels
            % for every blob (these will be our head and tail candidates)
            for i = 1:length(bodyStats)
                pixels = bodyStats(i).PixelList;
                maxDist = 0;
                for j = 1:size(pixels,1)-1
                    for k = j+1:size(pixels,1)
                        far = sum((pixels(j,:) - pixels(k,:)).^2);
                        if far > maxDist
                            maxDist = far;
                            bodyEndsA(i,:) = pixels(j,:);
                            bodyEndsB(i,:) = pixels(k,:);
                        end
                    end
                end
            end
        end
        
        
        %% handling of extracted features vectors
        function [] = makeMovieFromFV2(pathOfFV, outputPath, expName, nSubjects)
            % create movies composed of features vectors (one movie per mouse) in order
            % to save memory since the movies are compressed.
            % WARNING: lossy
            
            % pathOfFV = 'X:\sniffing library\bigFeatsVecsVar\';
            fv = load([pathOfFV expName '.FeatsVecs.mat']); % fv = load([pathOfFV 'bigFeatsVecsVarMat2.mat']);
            fv = fv.featsVecs;
            
            for j = 1 : nSubjects
                eval(['movieForMouse' num2str(j) ' = {};']);
            end
            
            for i = 1 : size(fv,1)/nSubjects % 400
                for j = 1 : nSubjects
                    newFrame = reshape(fv((i-1)*nSubjects+j, :), [41, 62, 3]);
                    newFrame = im2frame(newFrame);
                    eval(['movieForMouse' num2str(j) '{end+1} = newFrame.cdata;' ]); % im2frame(newFrame);']);
                end
            end
            
            for j = 1 : nSubjects
                %     eval(['v = VideoWriter(''X:\sniffing library\bigFeatsVecsVar\movieFor' num2str(j) '.avi'');']);
                eval(['v = VideoWriter(' outputPath '\' expName '.ForMouse' num2str(j) '.avi'');']);
                open(v);
                for k = 1 : length(movieForMouse1)
                    eval(['writeVideo(v, movieForMouse' num2str(j) '{' num2str(k) '}' ');']);
                end
                close(v);
            end
            
            % for j = 1 : nSubjects
            %     eval(['movieForMouse' num2str(j) ' = {};']);
            % end
            %
            % for i = 1 : size(fv,1)/nSubjects % 400
            %     for j = 1 : 1 % nSubjects
            %         newFrame = reshape(fv((i-1)*nSubjects+j, :), [41, 62, 3]);
            %         eval(['movieForMouse' num2str(j) '{end+1} = newFrame;' ]); % im2frame(newFrame);']);
            %     end
            % end
            %
            % for j = 1 : 1 % nSubjects
            %     eval(['v = VideoWriter(''X:\sniffing library\bigFeatsVecsVar\movieFor' num2str(j) 'Uncompressed.avi'', ''Uncompressed AVI'');']);
            %     v.FrameRate = 25;
            %     open(v);
            %
            %     for k = 1 : length(movieForMouse1)
            %         eval(['writeVideo(v, movieForMouse' num2str(j) '{' num2str(k) '}' ');']);
            %     end
            %     close(v);
            % end
            
        end
        
        function [] = concatFeatsVecsData(resPath, expName, numOfFrames, savePath)
            
            % savePath = 'X:\sniffing library\bigFeatsVecsVar\';
            
            allFeatsVecs = struct();
            allFeatsVecs.featsVecs = cell(numOfFrames,1);
            allFeatsVecs.ab = cell(numOfFrames, 1);
            allFeatsVecs.ids = cell(numOfFrames, 1);
            
            filenames = dir(fullfile(resPath, ['*' expName '.segm*.mat' ]));
            ns = length(filenames);
            
            tic
            for i = 1 : ns-1
                act = load([filenames(i).folder '\' filenames(i).name]); act = act.cents;
                allFeatsVecs.featsVecs((i-1)*floor(numOfFrames/ns)+1:i*floor(numOfFrames/ns), :) = act.featsVecsMat;
                allFeatsVecs.ab((i-1)*floor(numOfFrames/ns)+1:i*floor(numOfFrames/ns), :) = act.abMat;
                allFeatsVecs.ids((i-1)*floor(numOfFrames/ns)+1:i*floor(numOfFrames/ns), :) = act.ids;
                act = rmfield(act, 'featsVecsMat');
                act = rmfield(act, 'abMat');
                act = rmfield(act,'ids');
                cents = act;
                save([filenames(i).folder '\' filenames(i).name], 'cents');
            end
            
            act = load([filenames(i).folder '\' filenames(ns).name]); act = act.cents;
            allFeatsVecs.featsVecs(i*floor(numOfFrames/ns)+1:end, :) = act.featsVecsMat;
            allFeatsVecs.ab(i*floor(numOfFrames/ns)+1:end, :) = act.abMat;
            allFeatsVecs.ids(i*floor(numOfFrames/ns)+1:end, :) = act.ids;
            act = rmfield(act, 'featsVecsMat');
            act = rmfield(act, 'abMat');
            act = rmfield(act,'ids');
            cents = act;
            save([filenames(i).folder '\' filenames(i).name], 'cents');
            
            
            
            % save([savePath 'bigFeatsVecsVar.mat'], '-struct', 'allFeatsVecs', '-v7.3');
            save([savePath expName '.FeatsVecs.mat'], '-struct', 'allFeatsVecs', '-v7.3');
            toc
        end
        
        function [] = makeMovieFromFV(pathOfFV, outputPath, expName, nSubjects)
            % create movies composed of features vectors (one movie per mouse) in order
            % to save memory since the movies are compressed.
            % WARNING: lossy
            
            % pathOfFV = 'X:\sniffing library\bigFeatsVecsVar\';
            fvFile = load([pathOfFV expName '.FeatsVecs.mat']); % fv = load([pathOfFV 'bigFeatsVecsVarMat2.mat']);
            fv = fvFile.featsVecs;
            ids = fvFile.ids;
            
            for j = 1 : nSubjects
                eval(['movieForMouse' num2str(j) ' = {};']);
            end
            
            for i = 1 : size(fv,1) % 400
                for j = 1 : size(fv{i}, 1)
                    newFrame = reshape(fv{i}(j, :), [41, 62, 3]);
                    newFrame = im2frame(newFrame);
                    eval(['movieForMouse' num2str(ids{i}(j)) '{end+1} = newFrame.cdata;' ]); % im2frame(newFrame);']);
                end
            end
            
            for j = 1 : nSubjects
                %     eval(['v = VideoWriter(''X:\sniffing library\bigFeatsVecsVar\movieFor' num2str(j) '.avi'');']);
                eval(['v = VideoWriter(''' outputPath expName '.ForMouse' num2str(j) '.avi'');']); % eval(['v = VideoWriter(' outputPath '\' expName '.ForMouse' num2str(j) '.avi'');']);
                open(v);
                eval(['actMovieForMouse = movieForMouse' num2str(j) ';']);
                for k = 1 : length(actMovieForMouse)
                    eval(['writeVideo(v, movieForMouse' num2str(j) '{' num2str(k) '}' ');']);
                end
                close(v);
            end
        end
        
        function [] = makeMoviesFromSeg(pathOfSegs, outputPath, expName, nSubjects, numOfFrames)
            % iterates over all the segments from the CheeseColorSegment
            % step ('Res' directory) and make movies from the features
            % vectors (one by subject). This is supposed to be temporary
            % and was implemented in order to save memory since these
            % videos are less heavy than the matrices of the features
            % vectors themselves.
            % It also concatenates all the coordinates of the heads and
            % tails into a single matrix and all the corresponding IDs
            % into another one and save both these matrices into a .mat
            % file (with the 'Featsvecs' suffix)
            allFeatsVecs = struct();    % will contain the bodyEnds coordinates (ab/abMat) and their respective IDs
            
            allFeatsVecs.ab = cell(numOfFrames, 1);
            allFeatsVecs.ids = cell(numOfFrames, 1);
            fvDim = [41, 62, 3]; % that's the size of our features vectors images (before they were reshaped to a single row vector)
            
            filenames = dir(fullfile(pathOfSegs, ['*' expName '.segm*.mat' ]));
            ns = length(filenames);
            
            for j = 1 : nSubjects   % creates a VideoWriter for each mouse
                eval(['v' num2str(j) ' = VideoWriter(''' outputPath expName '.FVMouse' num2str(j) '.avi'');']);
                eval(['open(v' num2str(j) ');']);
            end
            
            for i = 1 : ns-1        % iterates over all the segments of the actual experiment (in the 'Res' directory) besides the last one
                act = load([filenames(i).folder '\' filenames(i).name]); act = act.cents;
                fv = act.featsVecsMat;
                actIds = act.ids;
                allFeatsVecs.ab((i-1)*floor(numOfFrames/ns)+1:i*floor(numOfFrames/ns), :) = act.abMat;
                allFeatsVecs.ids((i-1)*floor(numOfFrames/ns)+1:i*floor(numOfFrames/ns), :) = act.ids;
                act = rmfield(act, 'featsVecsMat');
                act = rmfield(act, 'abMat');
                act = rmfield(act,'ids');
                cents = act; %#ok<NASGU>
                save([filenames(i).folder '\' filenames(i).name], 'cents');
                
                % iterates over all the features vectors which were saved
                % for this very movie segment then write each one of them
                % to the good video 
                for k = 1 : size(fv,1)
                    for j = 1 : size(fv{k}, 1)
                        newFrame = reshape(fv{k}(j, :), fvDim);
                        newFrame = im2frame(newFrame); %#ok<NASGU>
                        eval(['writeVideo(v' num2str(actIds{k}(j)) ', newFrame);']);
                    end
                end
            end
            
            % the last segment isn't in the loop since it may contain less
            % frames than the previous segments
            act = load([filenames(ns).folder '\' filenames(ns).name]); act = act.cents;
            fv = act.featsVecsMat;
            actIds = act.ids;
            allFeatsVecs.ab((ns-1)*floor(numOfFrames/ns)+1:end, :) = act.abMat;
            allFeatsVecs.ids((ns-1)*floor(numOfFrames/ns)+1:end, :) = act.ids;
            act = rmfield(act, 'featsVecsMat');
            act = rmfield(act, 'abMat');
            act = rmfield(act,'ids');
            cents = act; %#ok<NASGU>
            save([filenames(ns).folder '\' filenames(ns).name], 'cents');
            
            for k = 1 : size(fv,1)
                for j = 1 : size(fv{k}, 1)
                    newFrame = reshape(fv{k}(j, :), fvDim);
                    newFrame = im2frame(newFrame); %#ok<NASGU>
                    eval(['writeVideo(v' num2str(actIds{k}(j)) ', newFrame);']);
                end
            end
            
            for j = 1 : nSubjects
                eval(['close(v' num2str(j) ');']);
            end
            
            save([pathOfSegs expName '.FeatsVecs.mat'], '-struct', 'allFeatsVecs', '-v7.3');
        end
        
    end
end