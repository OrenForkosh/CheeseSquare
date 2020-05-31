classdef Hierarchy
    % Hierarchy related tools
    methods (Static = true)
        
        function r = num2rank(n)
            % Convert number to string rank ('alpha', 'beta', etc)
            names = {'\alpha', '\beta', '\gamma', '\delta'};
            r = names(n);
        end
        
        function [NormDSonDij, NormDSonPij, DS] = DavidScore(varargin)
            % Compute Normalized David's score to estimate dominance
            % Based upon: Vries, Han de, Jeroen M. G. Stevens, and Hilde 
            % Vervaecke. “Measuring and Testing the Steepness of Dominance 
            % Hierarchies.” Animal Behaviour 71, no. 3 (March 2006): 585–92. 
            % https://doi.org/10.1016/j.anbehav.2005.05.015.
            %
            %   NormDSonDij = DavidScore(cheese) Computes David's score 
            %   based on behavior in the 'cheese' CheeseSquare object.
            %
            %   NormDSonDij = DavidScore(Sij) Computes David's score 
            %   based on the chase-escape matrix Sij (where each row represents
            %   the chaser, and the cols the "escaper")
            %
            if isa(varargin{1}, 'CheeseSquare')
                obj = varargin{1};
                nij = obj.Hierarchy.AggressiveChase.ChaseEscape + obj.Hierarchy.AggressiveChase.ChaseEscape';
                Sij = obj.Hierarchy.AggressiveChase.ChaseEscape;
            else
                Sij = varargin{1};
                nij = Sij + Sij';
            end
            
            aux = @(Pij) {sum(Pij, 2), Pij * sum(Pij, 2), sum(Pij, 1)', Pij' * sum(Pij, 1)'};
            
            Pij = Sij./nij;
            Pij(nij == 0) = 0;
            WL = aux(Pij);
            [w, w2, l, l2] = WL{:};
            DS = w + w2 - l - l2;
            NormDSonPij = (DS + size(Sij,1) * (size(Sij,1) - 1)/2) / size(Sij,1);
            
            Dij = (Sij + .5) ./ (nij + 1);
            Dij(nij == 0) = 0;
            WL = aux(Dij);
            [w, w2, l, l2] = WL{:};
            DS = w + w2 - l - l2;
            NormDSonDij = (DS + size(Sij,1) * (size(Sij,1) - 1)/2) / size(Sij,1);
        end
        
        function [s, ds] = DavidStatus(varargin)
            % Gives the mouse rank according to David's score. Same
            % parameters as 'DavidScore'
            ds = Hierarchy.DavidScore(varargin{:});
            [~, o] = sort(ds, 'descend');
            seq = 1:length(ds);
            s(o) = seq;
        end

        function [rank, removed]  = Rank(varargin)
            % Compute social rank based on binomial test. 
            % 
            %   rank = Rank(src) Computes the social rank where src is
            %   either a CheeseSquare object or the chase-escape matrix
            if isa(varargin{1}, 'CheeseSquare')
                obj = varargin{1};
                if nargin < 2
                    ce = obj.Hierarchy.AggressiveChase.ChaseEscape;
                else
                    ce = varargin{2};
                end
            else
                ce = varargin{1};
            end
            mat = ce - ce';
            mat(binotest(ce, ce + ce')) = 0;
            mat = mat .* (mat > 0);
            
            [rank, removed] = TopoFeedbackArcSetHierarchicalOrder(mat);
        end
        
        function c = Approaches(obj, range)
            % Compute the approaches matrix
            %
            %   c = Approaches(obj, range) Where obj is a CheeseSquare object.
            %   Range is optional and represent on which frames to computer
            %   the approaches
            if nargin > 1
                map = Q.torow(Q.overlap(Q.tocol(min(cat(2, obj.Hierarchy.Contacts.List.beg))), Q.tocol(max(cat(2, obj.Hierarchy.Contacts.List.end))), range(1), range(2)));
            else
                map = true(1, length(obj.Hierarchy.Contacts.List));
            end
            c = zeros(obj.nSubjects);
            for l=find(map)
                curr = obj.Hierarchy.Contacts.List(l);
                s1 = curr.subjects(1);
                s2 = curr.subjects(2);
                c(s1, s2) = c(s1, s2) + curr.states(1, 2);
                c(s2, s1) = c(s2, s1) + curr.states(2, 2);
            end
        end

        function c = Leaves(obj, range)
            % Compute the number of times the mouse went away from an
            % interaction (produces a matrix)
            %
            %   c = Leaves(obj, range) Where obj is a CheeseSquare object.
            %   Range is optional and represent on which frames to computer
            %   the approaches
            if nargin > 1
                map = Q.torow(Q.overlap(Q.tocol(min(cat(2, obj.Hierarchy.Contacts.List.beg))), Q.tocol(max(cat(2, obj.Hierarchy.Contacts.List.end))), range(1), range(2)));
            else
                map = true(1, length(obj.Hierarchy.Contacts.List));
            end
            c = zeros(obj.nSubjects);
            for l=find(map)
                curr = obj.Hierarchy.Contacts.List(l);
                s1 = curr.subjects(1);
                s2 = curr.subjects(2);
                c(s1, s2) = c(s1, s2) + (curr.states(1, 6) == 1 &  curr.states(2, 5) == 0);
                c(s2, s1) = c(s2, s1) + (curr.states(2, 6) == 1 &  curr.states(1, 5) == 0);
            end
        end
        
        function c = ChaseEscape(obj, range)
            % Computes the chase-escape matrix
            % See 'Hierarchy.Approaches' for parameters            
            if nargin > 1
                map = Q.torow(Q.overlap(Q.tocol(min(cat(2, obj.Hierarchy.Contacts.List.beg))), Q.tocol(max(cat(2, obj.Hierarchy.Contacts.List.end))), range(1), range(2)));
            else
                map = true(1, length(obj.Hierarchy.Contacts.List));
            end
            agg = obj.Hierarchy.Contacts.Behaviors.AggressiveChase.Map;
            agg = [agg, false(1, length(map) - length(agg) + 1)];
            c = zeros(obj.nSubjects);
            for l=find(map)
                curr = obj.Hierarchy.Contacts.List(l);
                s1 = curr.subjects(1);
                s2 = curr.subjects(2);
                c(s1, s2) = c(s1, s2) + (curr.states(1, 5) & agg(l)) * 1;
                c(s2, s1) = c(s2, s1) + (curr.states(2, 5) & agg(l)) * 1;
            end
        end

        function [c, duration] = Contacts(obj, range)
            % Computes the contact matrix
            % See 'Hierarchy.Approaches' for parameters            
            if nargin > 1
                map = Q.torow(Q.overlap(Q.tocol(min(cat(2, obj.Hierarchy.Contacts.List.beg))), Q.tocol(max(cat(2, obj.Hierarchy.Contacts.List.end))), range(1), range(2)));
            else
                map = true(1, length(obj.Hierarchy.Contacts.List));
            end
            c = zeros(obj.nSubjects);
            duration = zeros(obj.nSubjects);
            for l=find(map)
                curr = obj.Hierarchy.Contacts.List(l);
                s1 = curr.subjects(1);
                s2 = curr.subjects(2);
                c(s1, s2) = c(s1, s2) + 1;
                c(s2, s1) = c(s2, s1) + 1;
                duration(s1, s2) = duration(s1, s2) + (curr.end(1) - curr.beg(1) + 1);
                duration(s2, s1) = duration(s2, s1) + (curr.end(2) - curr.beg(2) + 1);
            end
        end
        
        function [timeout, pout] = TimeOutside(obj)
            % Amount of time each pair of mice spend outside together
            fps = obj.Video.FrameRate;
            if fps == 0
                try
                    fps = median(1./diff(obj.Tracking.time));
                catch
                end
                if fps == 0
                    fps = 25;
                end
            end
            
            %%
            pout = zeros(obj.nSubjects);
            timeout = zeros(obj.nSubjects);
            out = ~obj.Tracking.sheltered(:, obj.Tracking.valid);
            for i=1:obj.nSubjects
                for j=i+1:obj.nSubjects
                    pout(i,j) = mean(out(i, :) & out(j, :));
                    pout(j,i) = pout(i,j);
                    timeout(i, j) = sum(out(i, :) & out(j, :)) / fps;
                    timeout(j, i) = timeout(i, j);
                end
            end

        end
        
        function Elo = EloRating(objs)
            % Elo rating dominance score - a time continous dominance rank
            % Based upon :
            %    Christof Neumann, Julie Duboscq, Constance Dubuc, Andri Ginting,
            %    Ade Maulana Irwan, Muhammad Agil, Anja Widdig, Antje Engelhardt,
            %    Assessing dominance hierarchies: validation and advantages of
            %    progressive evaluation with Elo-rating, Animal Behaviour,
            %    Volume 82, Issue 4, October 2011, Pages 911-921, ISSN 0003-3472,
            %    10.1016/j.anbehav.2011.07.016.
            
            if isa(objs, 'CheeseFarm')
                objs = objs.Sources;
            end
            
            if ~isa(objs, 'cell')
                temp = objs;
                objs = {};
                objs{1} = temp;
            end
            obj = objs{1};
            
            Elo.WindowSize = 2; % hours
            Elo.DayDuration = 12;
            Elo.k = 100;
            Elo.Scores = ones(obj.nSubjects, 1) * 1000;
            Elo.Times = find(obj.Tracking.valid, 1);
            Elo.Day = 1;
            Elo.Map = false(obj.nSubjects, 1);
            starttime = 0;
            allStartTimes = zeros(1, length(objs)+1);
            %%
            idx = 1;
            for oid=1:length(objs)
                obj = objs{oid};
                map = obj.Hierarchy.Contacts.Behaviors.AggressiveChase.Map;
                Elo.Scores = [Elo.Scores, zeros(obj.nSubjects, sum(map))];
                Elo.Times = [Elo.Times, zeros(1, sum(map))];
                Elo.Day = [Elo.Day, zeros(1, sum(map))];
                Elo.Map = [Elo.Map, false(obj.nSubjects, sum(map))];
                for i=find(map)
                    %   obj.Contacts.List.beg
                    %    obj.Contacts.Behaviors.AggressiveChase.Chaser(i)
                    Pred = obj.Hierarchy.Contacts.Behaviors.AggressiveChase.Chaser(i);
                    Prey = obj.Hierarchy.Contacts.Behaviors.AggressiveChase.Escaper(i);
                    
                    d = Elo.Scores(Pred, idx) - Elo.Scores(Prey, idx);
                    ploss = 1 ./ (1 + 10 .^ (d/400));
                    Elo.Scores(Pred, idx + 1) = Elo.Scores(Pred, idx) + ploss * Elo.k;
                    Elo.Scores(Prey, idx + 1) = Elo.Scores(Prey, idx) - ploss * Elo.k;
                    other = 1:obj.nSubjects;
                    other = other(other ~= Pred);
                    other = other(other ~= Prey);
                    Elo.Scores(other, idx + 1) = Elo.Scores(other, idx);
                    Elo.Times(idx + 1) = min(obj.Hierarchy.Contacts.List(i).beg) + starttime;
                    Elo.Day(idx + 1) = oid;
                    Elo.Map(Pred, idx + 1) = true;
                    Elo.Map(Prey, idx + 1) = true;
                    idx = idx + 1;
                end
                %%
                % for s=obj.nSubjects
                %     plot(Elo.Times, Elo.Scores', '.-')
                % end
                starttime = starttime + obj.Video.NumberOfFrames;
                allStartTimes(oid+1) = starttime;
            end

            last = zeros(1, objs{1}.nSubjects);
            Elo.Segments = struct();
            for day=Q.torow(unique(Elo.Day))
                realt = Elo.Times(day == Elo.Day) - allStartTimes(day);
                scores = Elo.Scores(:, day == Elo.Day);
                Elo.Segments(day).Times = ((Elo.WindowSize:Elo.WindowSize:Elo.DayDuration) + (day-1) * Elo.DayDuration)*obj.Video.FrameRate*3600;
                t = ((Elo.WindowSize:Elo.WindowSize:Elo.DayDuration))*obj.Video.FrameRate*3600;
                for i=1:length(t)
                    idx = find(t(i) > realt, 1, 'last');
                    if isempty(idx)
                        Elo.Segments(day).RealTimes(i) = nan;
                        Elo.Segments(day).Scores(i, :) = last;
                    else
                        Elo.Segments(day).RealTimes(i) = realt(idx) + allStartTimes(day);
                        Elo.Segments(day).Scores(i, :) = scores(:, idx);
                    end
                end
                if size(scores, 2) > 0
                    last = scores(:, end);
                end
            end
            
            %%
            Elo.Rank = zeros(size(Elo.Scores));
            seq = 1:obj.nSubjects;
            r = zeros(1, obj.nSubjects);
            for i=1:size(Elo.Scores, 2)
                [~, o] = sort(Elo.Scores(:, i), 'descend');
                r(o) = seq;
                Elo.Rank(:, i) = r;
            end
            %%
            if nargout == 0
                %%
                %     plot(Elo.Times, Elo.Scores, 'color');
                %     ax
                
                cmap = CheeseSquare.MiceColors('PRBYW');
                scores = cat(2, Elo.Segments.Scores);
                for s=1:obj.nSubjects
                    subplot(2,1,1);
                    plot(Elo.Times/obj.Video.FrameRate/60/60, Elo.Scores(s, :), 'color', cmap(s, :));
                    hon;
                    plot(cat(2, Elo.Segments.Times)/obj.Video.FrameRate/60/60, scores(s, :), 'o-', 'MarkerFaceColor', cmap(s, :), 'MarkerEdgeColor', 'k', 'color', 'k');
                    subplot(2,1,2);
                    %plot(Elo.Median.Times/obj.Video.FrameRate/60/60, Elo.Median.Scores(s, :), 'o-', 'MarkerFaceColor', cmap(s, :), 'MarkerEdgeColor', 'k', 'color', 'k');
                    plot(cat(2, Elo.Segments.Times)/obj.Video.FrameRate/60/60, scores(s, :), 'o-', 'MarkerFaceColor', cmap(s, :), 'MarkerEdgeColor', 'k', 'color', 'k');
                    hon
                end
                for i=[true find(diff(Elo.Day) ~= 0)]
                    for s=1:2
                        subplot(2,1,s);
                        t = Elo.Times(i)/obj.Video.FrameRate/60/60;
                        a = axis;
                        Fig.VLine(t, 'color', Colors.LightGray);
                        text(t, a(4), sprintf(' day %d', Elo.Day(i+1)), 'verticalalignment', 'top', 'horizontalalignment', 'left', 'color', Colors.DarkGray);
                    end
                end
                
                subplot(2,1,1);
                Fig.Labels('time [hours]', 'score')
                Fig.Fix
                hoff;
                subplot(2,1,2);
                Fig.Labels('time [hours]', 'score')
                Fig.Fix
                hoff
            end
        end
        
        function ce = FromTable(tt)
            % Infers chase-escape matrix from profile table
            warning('using PairwiseChaseRate');
            %%
            ce = nan(max(tt.MouseID), max(tt.MouseID), max(tt.GroupNumber));
            for g=Q.torow(unique(tt.GroupNumber))
                ce(:, :, g) = Q.accumrows(tt.MouseID(tt.GroupNumber == g), tt.PairwiseChaseRate(tt.GroupNumber == g, :), @sum);
            end
            %%
        end
        
        function Show(obj)
            % Display Elo rating data for group
            for type=1:2
                dayid = obj.Hierarchy.Group.AggressiveChase.DayID;
                cds = zeros(obj.nSubjects, length(dayid)); % cummulative David-score
                dds = zeros(obj.nSubjects, length(dayid)); % daily David-score
                if type == 1
                    CE = obj.Hierarchy.Group.AggressiveChase.ChaseEscape;
                else
                    CE = obj.Hierarchy.Group.ChaseEscape.ChaseEscape;
                end
                
                for i=1:length(dayid)
                    ce = sum(CE(:, :, 1:i), 3);
                    cds(:, i) = Hierarchy.DavidScore(ce);
                    dds(:, i) = Hierarchy.DavidScore(CE(:, :, i));
                    %%
                    %                 subplot(2,2,3);
                    %                 Fig.Square(length(dayid), i);
                    %                 Plot.ImagePatch(CE(:, :, i), Colormaps.BlueRed, [0, 20]);
                    %                 title(sprintf('day %d', dayid(i)));
                end
                %%
                subplot(2,2,((type-1)*2) + 1);
                cmap = CheeseSquare.MiceColors('PRBYW');
                for i=1:obj.nSubjects
                    plot(dayid, cds(i, :), 'o-', 'markeredgecolor', 'none', 'markerfacecolor', cmap(i, :), 'color', cmap(i, :));
                    Fig.Hon
                end
                Fig.Hoff
                Fig.Fix;
                Fig.Labels('day', 'cummulative David-score');
                if type == 1
                    Fig.Title('Aggressive');
                else
                    Fig.Title('Aggressive & non-aggressive');
                end

                %%
                subplot(2,2,((type-1)*2) + 2);
                cmap = CheeseSquare.MiceColors('PRBYW');
                for i=1:obj.nSubjects
                    plot(dayid, dds(i, :), 'o-', 'markeredgecolor', 'none', 'markerfacecolor', cmap(i, :), 'color', cmap(i, :));
                    Fig.Hon
                end
                Fig.Hoff
                Fig.Fix;
                Fig.Labels('day', 'daily David-score');
            end
        end
        
        function Line(varargin)
            % Show group structure based on David's score (experimental. Avoid!)
            %
            %   Line(src) where src is either the profile table of
            %   chase-escape matrix
            if isa(varargin{1}, 'CheeseSquare')
            elseif istable(varargin{1})
                t = varargin{1};
                %%
                [~, ~, mn] = unique(t.MouseNumber);
                CE = t.PairwiseChaseCount;
                CE = Q.accumrows(mn, CE, @sum);
                %%
                gn = Q.accumrows(mn, t.GroupNumber, @mode);
                ugn = unique(gn);
                [~, ~, gt] = unique(t.GroupType);
                %if length(ugn) > 1
                gmap = lines;
                cmap = [CheeseSquare.MiceColors('PRBY'); lines];
                label = {};
                    idx = 1;
                    for g=ugn(:)'
                        %Fig.Square(length(ugn), idx);
                        %
                        ce = CE(gn == g, :);
                        ds = Hierarchy.DavidScore(ce);
                        grouptype = unique(t.GroupType(t.GroupNumber == g));
                        plot([idx idx], Q.minmax(ds), 'k', 'color', gmap(unique(gt(t.GroupNumber == g)), :));
                        hold on;
                        for i=1:length(ds)
                            %Patches.Circle(idx, ds(i), .25, cmap(i, :));
                            plot(idx, ds(i), 'o', 'MarkerFaceColor', cmap(i, :), 'MarkerEdgeColor', 'none', 'MarkerSize', 10);
                        end
                        %
                        
                        label{idx} = sprintf('%s exp%04d', grouptype{1}, unique(t.GroupID(t.GroupNumber == g)));
                        idx = idx + 1;
                    end
                    hold off;
                    %axis equal;
                    set(gca, 'XTick', 1:length(ugn), 'XTickLabel', label, 'XTickLabelRotation', -30)
                    set(gca, 'XMinorTick', 'off')
                    Fig.Fix
                    grid on
                %end
            elseif ismatrix(varargin{1})
                ce = varargin{1};
                ds = Hierarchy.DavidScore(ce);
                cmap = [CheeseSquare.MiceColors('PRBY'); lines];
                plot([1 1], Q.minmax(ds), 'k-');
                hold on;
                for i=1:length(ds)
                    Patches.Circle(1, ds(i), .25, cmap(i, :));
                end
                hold off;
                axis equal;
                axis off
            end
        end
            
        function Tree(varargin)
            % Tree Show group structure based on binomial test
            %
            %   Tree(src) Where src is either a CheeseSquare object,
            %   profile table or chase-escape matrix
            %
            if nargin < 1
                error;
            end
            if isa(varargin{1}, 'CheeseSquare')
                Hierarchy.Graph(varargin{:});
            elseif istable(varargin{1})
                t = varargin{1};
                %%
                [~, ~, mn] = unique(t.MouseNumber);
                ce = t.PairwiseChaseCount;
                ce = Q.accumrows(mn, ce, @sum);
                %%
                gn = Q.accumrows(mn, t.GroupNumber, @mode);
                ugn = unique(gn);
                if length(ugn) > 1
                    idx = 1;
                    for g=ugn(:)'
                        Fig.Square(length(ugn), idx);
                        Hierarchy.Graph([], ce(gn == g, :), varargin{2:end});
                        grouptype = unique(t.GroupType(t.GroupNumber == g));
                        Fig.Title(sprintf('%s exp%04d', grouptype{1}, unique(t.GroupID(t.GroupNumber == g))));
                        idx = idx + 1;
                    end
                end
                
                %%
            elseif ismatrix(varargin{1})
                Hierarchy.Graph([], varargin{:})
            else
                error;
            end
        end

        function Classic(obj, ce)
            % Show group structure based on binomial test
            %
            %   Classic(obj, ce) Where obj is a CheeseSquare object and ce
            %   is the chase-escape matrix
            if nargin < 2
                ce = obj.Hierarchy.AggressiveChase.ChaseEscape;
            end
            Hierarchy.Graph(obj, ce, [], [], false);
        end
        
        function DavidScoreGraph(varargin)
            % Show group structure based on David's score
            %   DavidScoreGraph(src) where src is either a CheeseSquare object,
            %   vector of David's scores or chase-escape matrix
            opt.Radius = .3;
            opt.Linear = true;
            opt.LineWidth = .1;
            %%
            if length(varargin) == 1 && isa(varargin{1}, 'CheeseSquare')
                obj = varargin{1};
                ce = obj.Hierarchy.AggressiveChase.ChaseEscape;
                ds = Hierarchy.DavidScore(obj, ce);
                varargin = varargin(2:end);
            end
            if length(varargin) == 1
                if isvector(varargin{1})
                    ds = varargin{1};
                    ce = [];
                else
                    ce = varargin{1};
                    ds = Hierarchy.DavidScore(obj, ce);
                end
            end
            nsubj = length(ds);

            %%
            cmap = [CheeseSquare.MiceColors('PRBY'); lines];
            rank = Q.rank(ds);
            lcolor = ones(1,3) *.95;
            Patches.Line([1 1], [0 nsubj-1], opt.LineWidth, lcolor);
            hold on
            Patches.Circle(1, 0, opt.LineWidth, lcolor);
            Patches.Circle(1, nsubj-1, opt.LineWidth, lcolor);
            for s=1:nsubj
                if ~opt.Linear
                    r = mod(rank(s)-1, 4);
                    switch r
                        case {0, 3}
                            Patches.Circle(1, ds(s), opt.Radius, cmap(s, :));
                        case 1
                            Patches.Circle(1-2*opt.Radius, ds(s), opt.Radius, cmap(s, :));
                        case 2
                            Patches.Circle(1+2*opt.Radius, ds(s), opt.Radius, cmap(s, :));
                    end
                else
                    Patches.Circle(1, ds(s), opt.Radius, cmap(s, :), 'FaceAlpha', .8);
                end
                hold on;
            end
            ylim([-.5, nsubj-.5]);
            axis equal
            hold off;
        end
        
        function Graph(obj, ce, interactions, scale, useDS)
            % Show group structure based on binomial test
            %
            %   Graph(obj, ce, interactions, scale, useDS) Where obj is a
            %   CheeseSquare object, ce (optional) is the chase-escape
            %   matrix, interactions (optional) the width of interactions
            %   lines (equals ce by default), scale (optional) scaling for
            %   the line widths, and useDS (optional, default false) use 
            %   David's score to determine rank
            radius = .3;
            dw = 1;
            if nargin < 2 || isempty(ce)
                ce = obj.Hierarchy.AggressiveChase.ChaseEscape;
            end
            if ~exist('interactions', 'var') || isempty(interactions)
                interactions = ce;
            end
            if nargin < 4 || isempty(scale)
                scale = 1500;
            end
            if nargin < 5
                useDS = false;
            end
            %%
            if useDS
                DS = Hierarchy.DavidScore(ce);
                rank = Q.rank(DS);
            else
                rank = Hierarchy.Rank(obj, ce);
                DS = rank;
            end
             %%
            ds = Hierarchy.DavidScore(ce);
            
            SO = max(rank) - rank;
            cmap = [CheeseSquare.MiceColors('PRBY'); lines];
            loc = zeros(length(rank), 2);
            for i=1:length(rank)
                so = SO(i);
                count = sum(SO == so);
                %offset = -(2*mod(so,2) - 1) * (count - 1 + (so > 0 & so < max(SO))) / 2 + sum(SO(1:i-1) == so);
                if count > 1
                    offset = -(count-1)/2 + sum(SO(1:i-1) == so);
                else
                    if so==0 || so == max(SO)
                        offset = 0;
                    else
                        offset = mod(so, 2)*2-1;
                    end
                end
                loc(i, :) = [offset * dw, DS(i)];
                %loc(i, :) = [offset * dw, mean(DS(SO == os))];
            end
            for i=1:size(loc, 1)
                for j=i+1:size(loc, 1)
                    if interactions(i, j) == 0 && interactions(j, i) == 0
                        continue;
                    end
                    if interactions(i, j) > interactions(j, i)
                        if interactions(i, j) > 0 
                            Patches.Line(loc([i,j], 1), loc([i,j], 2), interactions(i, j)/scale, cmap(i, :), radius/5);
                            Fig.Hon
                        end
                        if interactions(j, i) > 0 
                            Patches.Line(loc([j,i], 1), loc([j,i], 2), interactions(j, i)/scale, cmap(j, :), radius/5);
                            Fig.Hon
                        end
                    else
                        if interactions(j, i) > 0 
                            Patches.Line(loc([j,i], 1), loc([j,i], 2), interactions(j, i)/scale, cmap(j, :), radius/5);
                            Fig.Hon
                        end
                        if interactions(i, j) > 0 
                            Patches.Line(loc([i,j], 1), loc([i,j], 2), interactions(i, j)/scale, cmap(i, :), radius/5);
                            Fig.Hon
                        end
                    end
                    Fig.Hon
                end
            end
            for i=1:size(loc, 1)
               Patches.Circle(loc(i, 1), loc(i, 2), radius, cmap(i, :));
               text(loc(i, 1), loc(i, 2), sprintf('%.1f', ds(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'Color', Colors.CompLight(cmap(i, :)))
            end
            Fig.Hoff
            axis equal;
            axis off;
        end
        
        
    end
    
end