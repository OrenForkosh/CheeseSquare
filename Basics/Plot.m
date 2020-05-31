classdef Plot
    methods (Static = true)
        function Graph(d, varargin)
            %% Note! assuming symmetry
            p = inputParser;
            addOptional(p, 'AttractorStrength', 1); 
            addOptional(p, 'ColorScale', 2); 
            addOptional(p, 'LineScale', 2); 
            addOptional(p, 'CenterColors', Colors.PrettyBlue); 
            addOptional(p, 'CenterLabels', {}); 
            addOptional(p, 'MaxWidth', 5); 
            addOptional(p, 'Gap', 0); 
            p.parse(varargin{:});
            opt = p.Results;

            %%
            dim = size(d, 1);
            radius = 2/dim;
            maxwidth = opt.MaxWidth;
            annotate = true;
            maxval = max(abs(d(:)));
            
            centers = [cos(linspace(0, 2*pi - opt.Gap, dim+1)); -sin(linspace(0, 2*pi - opt.Gap, dim+1))]';
            
            % draw edges
            linemap = Colormaps.BlueWhiteRed;
            midmap = floor(size(linemap, 1)/2 - 1);
            
            scalfun = @(w) sign(w).*(exp(opt.ColorScale*abs(w))-1)/(exp(opt.ColorScale)-1);
            colorfun = @(w) linemap(round(midmap + 1 + scalfun(w) * midmap), :);

            linescalfun = @(w) sign(w).*(exp(opt.LineScale*abs(w))-1)/(exp(opt.LineScale)-1);
            widthfun = @(w) abs(linescalfun(w)) * maxwidth;
            
            [idx2, idx1] = meshgrid(1:size(d, 2), 1:size(d, 1));
            uptri = triu(true(size(d)),1);
            idx = [idx1(uptri), idx2(uptri)];
            [~, order] = sort(abs(d(uptri)));
            
            for i=order'
                %%
                fr = centers(idx(i, 1), :);
                to = centers(idx(i, 2), :);
                %%
                dfr = [sin(atan2(fr(1),fr(2)))/4, cos(atan2(fr(1),fr(2)))/4];
                dto = [sin(atan2(to(1),to(2)))/4, cos(atan2(to(1),to(2)))/4];
                %bezier = Shapes.Bezier(fr, fr - dfr, to - dto, to);
                
                bezier = Plot.QuadraticBezier(fr, to, [0, 0], 101, opt.AttractorStrength);
                w = max(min(d(idx(i, 1), idx(i, 2)) / maxval, 1), -1);
                if w ~= 0
                    plot(bezier(:, 1), bezier (:, 2), '-', 'Color',colorfun(w), 'LineWidth',  widthfun(w));
                    hold on
                end
            end
            
            % draw nodes
            for i=1:size(d, 1)
                color = opt.CenterColors(mod(i-1, size(opt.CenterColors, 1)) + 1, :);
                Patches.Circle(centers(i, 1), centers(i, 2), radius, color);
                if annotate 
                    if isempty(opt.CenterLabels)
                        text(centers(i, 1), centers(i, 2), num2str(i), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', Q.ifthen(Colors.Value(color) > .92, 'k', 'w'));
                    else
                        text(centers(i, 1), centers(i, 2), opt.CenterLabels{i}, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', Q.ifthen(Colors.Value(color) > .92, 'k', 'w'));
                    end
                    
                end
                hold on;
            end
            hold off;
            axis square
        end
        
        function Ternary(varargin)
            if size(varargin{1}, 2) == 3
                a = varargin{1}(:, 1);
                b = varargin{1}(:, 2);
                c = varargin{1}(:, 3);
                opt = varargin(2:end);
            elseif nargin >= 3 && isnumeric(varargin{1}) && isnumeric(varargin{2}) && isnumeric(varargin{3})
                a = varargin{1};
                b = varargin{2};
                c = varargin{3};
                opt = varargin(4:end);
            else
                error;
            end
            p = inputParser;
            addOptional(p, 'MarkerSize', 150, @isnumeric); 
            addOptional(p, 'EdgeColor', 'none'); 
            addOptional(p, 'FaceColor', [000 177 229] / 255); 
            addOptional(p, 'Alpha', .8); 
            p.parse(opt{:});
            opt = p.Results;
            
            x = 1/2*(2*b+c)./(a+b+c);
            y = sqrt(3)/2*c./(a+b+c);
            
            ish = ishold;
            if ~ish
                cla
                patch([0 1 1/2 0], [0 0 sqrt(3)/2 0], [.9 .9 .9], 'EdgeColor', [.8 .8 .8]);
                hold on;
            end
            scatter(x, y, opt.MarkerSize, opt.FaceColor, 'filled', 'MarkerFaceAlpha', opt.Alpha);
            if ~ish
                hold off;
            end
            axis off
        end
        
        function Scatter(x,y,varargin)
            % Scatter   Scatter plot.
            %   Scatter(X,Y) plots vector Y versus vector X. 
            %
            %   Scatter(X,Y, label) plots vector Y versus vector X. Points
            %   that belong to the same class have similar style. labels can
            %   be empty and then it is ignored.
            %
            %   Scatter(X,Y, label, errb) also shows errorbars on points. If
            %   errb is a vector, errors are assumed to be on the Y axis.
            %
            %   Scatter(X,Y...style) use style on plot (similar to plot
            %   options.
            
            p = inputParser;
            
            offset = find(cellfun(@ischar, varargin), 1);
            offset = Q.ifthen(isempty(offset), length(varargin)+1, offset);
            style = varargin(offset:end);
            if offset > 1
                label = varargin{1};
            else
                label = [];
            end
            if offset > 2
                errb = varargin{2};
            else
                errb = [];
            end
            %%
            ish = ishold;
            if ~ish
                set(gca, 'ColorOrderIndex',  1);
            end
            if isempty(label)
                label = ones(size(x));
            end
            [label, ~, id] = unique(label);
            nlabels = label(end);
            %cmap = Colormaps.SubCategorical(nlabels);
            for l=1:nlabels
                i = find(id == l);
                colorindex = get(gca, 'ColorOrderIndex');
                a = gca;
                if colorindex > size(a.ColorOrder, 1)
                    colorindex = mod(colorindex, size(a.ColorOrder, 1));
                end
                c = a.ColorOrder(colorindex, :);
                if ~isempty(errb)
                    if size(errb, 2) == 1
                        x1 = x(i) - errb(i, 1);
                        x2 = x(i) + errb(i, 1);
                        plot([x1 x2]', [y(i) y(i)]', '.-', 'Color', c, style{:}, 'HandleVisibility', 'off');
                        hold on;
                    end
                    if size(errb, 2) >= 2
                        y1 = y(i) - errb(i, 2);
                        y2 = y(i) + errb(i, 2);
                        %plot([x(i) x(i)]', [y1 y2]', '.-', 'Color', c,  style{:}, 'HandleVisibility', 'off');
                        for j=1:length(i)
                            Patches.Circle(x(i(j)), y(i(j)), [errb(i(j), 1), errb(i(j), 2)],c,  'HandleVisibility', 'off', 'FaceAlpha', .2);
                            hold on;
                        end
                    end
                end
                plot(x(i(1)), y(i(1)), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', c,  style{:});
                hold on;
                plot(x(i(2:end)), y(i(2:end)), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', c,  style{:}, 'HandleVisibility', 'off');
            end
            set(gca, 'ColorOrderIndex',  colorindex + 1);
            if ish
                hold on;
            else
                hold off;
            end
            Fig.Fix
        end
        
        function Line(x,y,width,c,varargin)
            plot(x,y,'color', c, 'linewidth', width, varargin{:});
        end

        function Hinton(m, varargin)
            %%
            p = inputParser;
            addOptional(p, 'annotate', false, @islogical); 
            addOptional(p, 'abs', true, @islogical); 
            addOptional(p, 'minmax', false, @islogical); 
            addOptional(p, 'mark', false(0), @islogical); % logical matrix showing which entries to highlight
            addOptional(p, 'colormap', [Colors.PrettyBlue; Colors.PrettyRed], @(x) isnumeric(x) && size(x, 2) == 3);
            addOptional(p, 'range', [], @(x) isnumeric(x) && length(x) == 2);
            p.parse(varargin{:});
            
            opt = p.Results;
            
            %%
            cla
            %mm = Q.minmax(m(:));
                csz = size(opt.colormap, 1);
            if isempty(opt.range)
                opt.range = quantile(m(:), [.025, .975]);
                if opt.abs
                    opt.range = max(abs(opt.range)) * [-1 1];
                end
            end
            r = max(min(m / (opt.range(2)-opt.range(1)), 1), -1);
            %if opt.abs
            %    cidx = floor(r * (csz/2) + (csz/2)) + 1;
            %else
                cidx = floor((m - opt.range(1)) / (opt.range(2)-opt.range(1)) * (csz)) + 1;
            %end
            cidx = max(min(cidx, csz), 1);
            
            WH = max(min((m - opt.range(1)) / (opt.range(2)-opt.range(1)), 1), -1);
            
            [mx, mxij] = Q.argmaxnd(m);
            [mn, mnij] = Q.argminnd(m);
            mx = m(mx);
            mn = m(mn);
            for i=1:size(m, 1)
                for j=1:size(m, 2)
                    wh = abs(r(i,j));
                    % color = Q.ifthen(m(i,j) < 0, Colors.PrettyBlue, Colors.PrettyRed);
                    color = opt.colormap(cidx(i, j), :);
                    if isnan(m(i, j))
                        continue;
                    end
                        if ~isempty(opt.mark) && opt.mark(i, j)
                            Patches.Rect(j-wh/2, i-wh/2, wh, wh, color, 'EdgeColor', 'k');
                        else
                            Patches.Rect(j-wh/2, i-wh/2, wh, wh, color);
                        end
                    hold on;
                    if opt.annotate
                        c = Colors.RGB2HSV(color);
                        c(2) = 1 - c(2);
                        c(3) = 1 - c(3);
                        c = Colors.HSV2RGB(c);
                        text(j, i, sprintf('%3.2g', m(i,j)), 'horizontalalignment', 'center', 'Color', c);
                    end
                end
            end
                
            if opt.minmax
                for i=1:size(m, 1)
                    for j=1:size(m, 2)
                        if all([i,j] == mxij)
                            text(j, i, sprintf('%3.1g', mx), 'horizontalalignment', 'center');
                        elseif all([i,j] == mnij)
                            text(j, i, sprintf('%3.1g', mn), 'horizontalalignment', 'center');
                        end
                    end
                end
            end
            
            
            hold off;
            axis equal
            xlim([0, size(m,2)+1]);
            ylim([0, size(m,1)+1]);
%             set(gca, 'XTick', 1:min(diff(get(gca, 'XTick'))):max(get(gca, 'XLim')));
%             set(gca, 'YTick', 1:min(diff(get(gca, 'YTick'))):max(get(gca, 'YLim')));
            set(gca, 'XTick', 1:size(m,2));
            set(gca, 'YTick', 1:size(m,1));

            % set(gca, 'YTick', get(gca, 'YTick') + (1 - min(get(gca, 'YTick'))));
            set(gca, 'YDir', 'reverse')
            Fig.Fix
            box on;
            
            grid on;
            return
            f = factor(size(m, 1)); f = [f prod(f)]; f = min(f(f >= 5));
            if ~isempty(f)
                for i=f:f:size(m, 1)-1
                    Fig.HLine(i+.5);
                end
            end
            f = factor(size(m, 2)); f = [f prod(f)]; f = min(f(f >= 5));
            if ~isempty(f)
                for i=f:f:size(m, 2)-1
                    Fig.VLine(i+.5);
                end
            end
        end
        
        function Lines(m)
            thresh = quantile(Q.tocol(diff(m)), 0.01);
            a = log10(abs(thresh));
            a = sign(a) * ceil(abs(a));
            jump = ceil(abs(thresh) / 10 ^ a) * 10 ^ a;
            plot([m + repmat(jump * (0:size(m, 1)-1)', 1, size(m, 2))]');
            xlim([0, size(m, 2)]);
        end
        
        function Clock(t, y, c, varargin)
            opt = struct(varargin{:});
            opt = Q.defaultargs(false, opt, ...
                'start', 10/24 ...
                );
            if nargin >= 3 && ~isempty(c)
                polar(2*pi-(t-opt.start)*2*pi + pi/2, y, c);
            else
                polar(2*pi-(t-opt.start)*2*pi + pi/2, y);
            end
            %%
            src = linspace(0, 330, 12);
            tgt = [opt.start*24 + 6:-2:2 24:-2:opt.start*24 + 4];
            valid = mod(tgt - opt.start*24, 12) == 0;
            for i = 1:length(src)
                r = src(i);
                if valid(i)
                    set(findall(gcf, 'string', num2str(r)), 'String', sprintf('%02d:00', tgt(i)));
                else
                    set(findall(gcf, 'string', num2str(r)), 'String', '');
                end
                
            end
        end
        
        function Network(weights, varargin)
            if ~isempty(varargin) && isnumeric(varargin{1})
                nodes = varargin{1};
                varargin = varargin(2:end);
            else
                r = 10;
                n = max(size(weights, 1), size(weights, 2));
                nodes = zeros(n, 2);
                for i=1:n
                    nodes(i, 1) = r*cos((i-1)/n*2*pi);
                    nodes(i, 2) = r*sin((i-1)/n*2*pi);
                end
            end
            opt = Q.defaultargs(false, varargin, ...
                'Colormap', [] ...
                , 'NodeColormap', lines(size(nodes, 1)) ...
                , 'NodeRadius', 1 ...
                , 'MaxEdgeWidth', .4 ...
                , 'CenterAttraction', 2 ...
                , 'Annotate', true ...
                , 'Range', [] ...
                );
            issym = issymmetric(weights) || all(Q.tocol(triu(weights) == weights));
            if isempty(opt.Range)
                opt.Range = quantile(abs(weights(:)), .95);
            end
            weights(eye(size(weights)) > 0) = 0;
            nweights = max(min(weights/opt.Range, 1), -1);
            weights = abs(nweights) * opt.MaxEdgeWidth;
            cc = mean(nodes);
            if isempty(opt.Colormap)
                if issym
                    opt.Colormap = flip(Colormaps.BlueWhiteRed, 1);
                else
                    opt.Colormap = lines(size(nodes, 1));
                end
            end
      
            for i=1:size(nodes, 1)
                Patches.Circle(nodes(i, 1), nodes(i, 2), opt.NodeRadius, opt.NodeColormap(i, :));
                Fig.Hon              
            end            
            [~, sorted] = sort(abs(weights(:)));
            %[i,j] = ind2sub(size(d), sorted);
            
            for idx=sorted(:)'
                [i, j] = ind2sub(size(weights), idx);
                if weights(i,j) == 0
                    continue;
                end
                %for j=1:size(nodes, 1)
                    x = nodes([i,j], 1);
                    y = nodes([i,j], 2);
                    a = atan2(y(2)-y(1), x(2)-x(1));
                    shift = opt.CenterAttraction;
                    mx = (x(2)+x(1))*1/2;
                    my = (y(2)+y(1))*1/2;
                    
                    mx1 = mx + shift * sin(a);
                    my1 = my - shift * cos(a);
                    
                    mx2 = mx - shift * sin(a);
                    my2 = my + shift * cos(a);
                    
                    if issym 
                        color = opt.Colormap(min(max(round((1+nweights(i, j)) * size(opt.Colormap,1)), 1), size(opt.Colormap, 1)), :);
                        
                        if pdist2([mx, my], cc) < min(pdist2([mx1, my1], cc), pdist2([mx2, my2], cc))
                            mx1 = mx;
                            my1 = my;
                        elseif pdist2([mx1, my1], cc) > pdist2([mx2, my2], cc) 
                            continue;
                        elseif pdist2([mx1, my1], cc) == pdist2([mx2, my2], cc)
                            if i<j
                                continue;
                            end
                            mx1 = mx;
                            my1 = my;
                        end
                    else
                        color = opt.Colormap(i, :);
                    end
                    
                    l = Plot.Bezier([x(1) y(1)], [mx1, my1], [mx1 my1],  [x(2) y(2)], linspace(0, 1, 10));
                    Patches.Curve(l(:, 1), l(:, 2), weights(i,j), color);
                    Fig.Hon
                %end
            end
            for i=1:size(nodes, 1)
                Patches.Circle(nodes(i, 1), nodes(i, 2), opt.NodeRadius, opt.NodeColormap(i, :));
                Fig.Hon              
            end
            if opt.Annotate
                height = diff(ylim);
                for i=1:size(nodes, 1)
                    text(nodes(i, 1), nodes(i, 2), num2str(i), 'Color', 'w', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized', 'FontSize', 1.5*opt.NodeRadius/height);
                end
            end
            Fig.Hoff
            axis equal;
            axis off;
        end

        
        function pts = QuadraticBezier(pt1, pt2, attractor, n, strength)
            t = linspace(0,1,n);
            if nargin < 5
                strength = 1;
            end
            alpha = strength + 1;
            pts = (kron((1-t).^alpha,pt1(:)) + kron(alpha*(1-t).*t,attractor(:)) + kron(t.^alpha,pt2(:)))';
        end
        
        function [Q]=Bezier(P0,P1,P2,P3,varargin)
            % Bezier interpolation for given four control points.
            % Each control point can be in N-Dimensional vector space.
            % Input:
            % P0,P1,P2,P3 : four control points of bezier curve
            %               control points can have any number of coordinates
            % t : vector values of t paramter at which bezier curve is evaluated
            %   (optional argument) default 101 values between 0 and 1.
            
            % Output:
            % Q evaluated values of bezier curves. Number of columns of Q are equal to
            % number of coordinates in control point. For example for 2-D, Q has two
            % columns. Column 1 for x value and column 2 for y values. Similarly for
            % 3-D, Q will have three columns

            % % % --------------------------------
            % % % Author: Dr. Murtaza Khan
            % % % Email : drkhanmurtaza@gmail.com
            % % % --------------------------------

            
            if (nargin<4)
                disp('Atleast four input arguments (four control points) are required');
                return
            end
            
            [r0 c0]=size(P0); [r1 c1]=size(P1); [r2 c2]=size(P2); [r3 c3]=size(P3);
            if (r0~=r1 || r0~=r2 || r0~=r3 || c0~=c1 || c0~=c2 || c0~=c3)
                disp('arg1,arg2,arg3,arg4 must be of equal size');
                return
            end
            
            %%% Default Values %%%
            t=linspace(0,1,101); % uniform parameterization
            defaultValues = {t};
            %%% Default Values %%%
            %%% Assign Valus %%%
            nonemptyIdx = ~cellfun('isempty',varargin);
            defaultValues(nonemptyIdx) = varargin(nonemptyIdx);
            [t] = deal(defaultValues{:});
            %%% Assign Valus %%%
            
            [r c]=size(t);
            if(r>1 && c>1)
                disp('arg5 must be a vector');
                return
            end
            
            % Equation of Bezier Curve, utilizes Horner's rule for efficient computation.
            % Q(t)=(-P0 + 3*(P1-P2) + P3)*t^3 + 3*(P0-2*P1+P2)*t^2 + 3*(P1-P0)*t + Px0
            
            c3 = -P0 + 3*(P1-P2) + P3;
            c2 = 3*(P0 - (2*P1)+P2);
            c1 = 3*(P1 - P0);
            c0 = P0;
            for k=1:length(t)
                Q(k,:)=((c3*t(k)+c2)*t(k)+c1)*t(k) + c0;
            end
        end            
            
        
        
        function Errorbar(varargin)
            styleidx = find(cellfun(@ischar, varargin), 1);
            if ~isempty(styleidx)
                style = varargin(styleidx:end);
                varargin = varargin(1:styleidx-1);
            else
                style = {};
            end
            
            p = inputParser;
            addOptional(p, 'FaceColor', 'w'); 
            addOptional(p, 'EdgeColor', 'k'); 
            p.parse(style{:});
            opt = p.Results;

            
            color = Colors.ParseHex('314e54');
            
            
            if length(varargin) < 3
                y = varargin{1};
                e = varargin{2};
                x = 1:length(y);
            else
                x = varargin{1};
                y = varargin{2};
                e = varargin{3};
            end
            valid = ~isnan(x(:)) & ~isnan(y(:)) & ~isnan(e(:));
            x = x(valid)';
            y = y(valid)';
            e = e(valid)';
            y1 = y+e;
            y2 = y-e;
            plot(x(:), y(:));
            Fig.Hon
            %Patches.Polygon([x; x(end:-1:1)], [y1; y2(end:-1:1)], 'k', 'FaceColor', 'w', 'facealpha', .7);
            Patches.Polygon([x; x(end:-1:1); x(1)], [y1; y2(end:-1:1); y1(1)], 'k', 'FaceColor', opt.FaceColor, 'FaceAlpha', .3, 'EdgeColor', opt.EdgeColor, 'edgealpha', .5);
            plot(x(:), y(:), '.-', 'markeredgecolor',opt.EdgeColor, 'markerfacecolor', opt.EdgeColor, 'color', opt.EdgeColor);
            Fig.Hoff
        end
        
        function Multiprob(x, h, offset, varargin)
            %%
            if ~ishold
                cla reset
                ish = false;
            else
                ish = true;
            end
            if nargin < 3
                offset = 1;
            end
            h = h / sum(h);
            maxwidth = 3;
            dx = median(diff(x));
            cmap = lines(size(h, 1));
            for j=1:size(h, 1)
                for i=1:length(x)
                    width = maxwidth * h(j, i);
                    Patches.Rect(x(i)-dx/2, j-width/2+offset-1, dx, width, cmap(j, :), 'FaceAlpha', h(j, i) / max(h(j, :)), varargin{:});
                    hold on
                end
            end
            if ~ish
                hold off;
            end
        end
        
        function P = Statistics(data, varargin)
            opt = Q.defaultargs(true, varargin,  ...
                'colormap', lines, ...
                'test', 'ttest2',... % 'ttest2', or 'ttest'
                'names', {}, ...
                'nbins', 20, ...
                'type', 'bar',...
                'ShowPoints', false...
                );

            if ~ishold
                cla reset
                ish = false;
            else
                ish = true;
            end
            w = .4;
            %
            cmap = opt.colormap;
            if ~iscell(data)
                temp = cell(1, size(data, 2));
                for i=1:length(data)
                    temp{i} = data(:, i);
                end
                data = temp;
                temp = {};
            end
            m = zeros(1, length(data));
            e = zeros(1, length(data));
            for i=1:length(data)
                m(i) = nanmean(data{i});
                e(i) = Q.nanstderr(data{i});
            end
            switch opt.type
                case 'bar'
                    for i=1:length(m)
                        patch([i-w i+w i+w i-w], [0 0 m(i) m(i)], 1, 'FaceColor', cmap(i, :), 'EdgeColor', 'none');
                    end
                    lw = 0.02; %diff(lim) / 200;
                    yrange = [inf -inf];
                    for i=1:length(m)
                        patch([i-lw i+lw i+lw i-lw], m(i) + [0 0 e(i) e(i)], 1, 'FaceColor', cmap(i, :), 'EdgeColor', 'none');
                        patch([i-lw i+lw i+lw i-lw], m(i) + [-e(i) -e(i) 0 0 ], 1, 'FaceColor', cmap(i, :), 'EdgeColor', 'none');
                        if m(i) < 0
                            patch([i-lw i+lw i+lw i-lw], [m(i) * [1 1] min(m(i) + e(i), 0) * [1 1]], 1, 'FaceColor', 'w', 'EdgeColor', 'none');
                        else
                            patch([i-lw i+lw i+lw i-lw], [max(m(i) - e(i), 0) * [1 1] m(i) * [1 1]], 1, 'FaceColor', 'w', 'EdgeColor', 'none');
                        end
                        hold on;
                        %%
                        %plot(i + randn(size(data{i})) * .1, data{i}, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', cmap(i, :));
                        yrange = [min(yrange(1), m(i) - e(i)) max(yrange(2), m(i) + e(i))];
                    end
                    %%
                    Fig.HLine(0, 'LineStyle', '-', 'color', 'k');
                case 'candle'
                    res = [];
                    group = [];
                    for i=1:length(data)
                        res = [res; data{i}(:)];
                        group = [group; ones(length(data{i}), 1) * i];
                    end
                    boxplot(res, group, 'Colors',lines, 'Notch','marker', 'Symbol','o', 'Widths',0.2);
                    yrange = ylim;
            end
            %%
            if opt.ShowPoints
            hold on;
            w = .02;
            for i=1:length(m)
                [h, bin] = histc(data{i}, linspace(min(data{i}), max(data{i}), opt.nbins));
                count = zeros(1, opt.nbins);
                for j=1:length(data{i})
                    plot(i + w * count(bin(j)) - w * (h(bin(j)) - 1) / 2, data{i}(j), 'o', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', cmap(i, :), 'MarkerSize', 5);
                    count(bin(j)) = count(bin(j)) + 1;
                end
            end
            hold off;
            end
            %%
            lim = ylim;
            pos = max(lim(2), yrange(2)) + (lim(2) - lim(1))/10;
            P = zeros(length(data));
            if strcmpi(opt.test, 'ttest2')
                lh = diff(lim) / 200;
                hold on;
                for i=1:length(data)
                    for j=i+1:length(data)
                        [h, P(i,j)] = ttest2(data{i}, data{j}, 'Alpha',0.05);
                        if ~isnan(h) && h
                            fprintf('# (%d, %d)*\n', i, j);
                            %patch([i j j i], [pos-lh, pos-lh, pos+lh, pos + lh], 'k');
                            plot([i j], [pos-lh, pos-lh], 'k-');
                            %text((i+j)/2, pos, '*');
                            pos = pos + diff(lim) / 10;
                        end
                    end
                    [h, P(i,i)] = ttest(data{i});
                end
                hold off;
                P = P + P';
                t = get(gca, 'xtick');
                set(gca, 'xtick', t(t == round(t)));
            elseif isnumeric(opt.test)
                for i=1:length(data)
                    h = opt.test(i);
                    if ~isnan(h) && h
                        %text(i, pos, '*', 'FontSize', 25);
                        %pos = pos + diff(lim) / 10;
                        fprintf('# (%d)*\n', i);
                    end
                end
            else
                for i=1:length(data)
                    h = ttest(data{i});
                    if ~isnan(h) && h
                        %text(i, pos, '*', 'FontSize', 25);
                        %pos = pos + diff(lim) / 10;
                        fprintf('# (%d)*\n', i);
                    end
                end
            end
            %%
            x = xlim;
            x(1) = 1-1.5*w;
            x(2) = length(m)+1.5*w;
            xlim(x);
            y = ylim;
            y(1) = min(y(1), yrange(1));
            y(2) = max(y(2), pos);
            ylim(y);
            %%
            Fig.HLine(0, 'linestyle', ':', 'color', 'k');
            %%
            set(gca, 'xtick', 1:length(data));
%            if ~isempty(opt.names)
                set(gca, 'xticklabel', {});
 %           end
            %%
            set(gca, 'Xcolor', 'w');
            y = ylim;
            style = Fig.DefaultLabelFont;
            for i=1:length(data)
                if ~isempty(opt.names)
                    s = opt.names{i};
                else
                    s = num2str(i);
                end
                text(i, y(1) - (y(2)-y(1))/50, s, 'horizontalalignment', 'center', 'verticalalignment', 'top', style{:}, 'color', cmap(i, :));
            end
            text();
            if ~ish
                hold off;
            end
            Fig.Fix
            set(gca, 'XMinorTick', 'off')
        end
        
        function Hist(varargin)
            %%
            p = find(cellfun(@ischar, varargin), 1);
            if ~isempty(p)
                opt = struct(varargin{p:end});
                varargin = varargin(1:p-1);
            else
                opt = struct();
            end
            opt = Q.defaultargs(false, opt, ...
                'n', 25 ...
                );
            
            m1 = min(cellfun(@min, varargin));
            m2 = max(cellfun(@max, varargin));
            x = linspace(m1, m2, opt.n);
            cmap = lines;
            X = [x(1:end-1); x(2:end)]; X = X(:);
            for i=1:length(varargin)
                h = histc(varargin{i}, x);
                h = h(1:end-1);
                h = h / sum(h);
                H = repmat(h(:)', 2, 1); H = H(:);
                plot(X, H, 'color', cmap(i, :));
                Fig.Hon
            end
            
            for i=1:length(varargin)
                h = histc(varargin{i}, x);
                h = h(1:end-1);
                h = h / sum(h);
                H = repmat(h(:)', 2, 1); H = H(:);
                %area(x, h, 'facecolor', cmap(i, :), 'edgecolor', 'none');
                patch([X(:)', X(end) X(1)], [H(:)' 0 0], 1, 'facecolor', cmap(i, :), 'edgecolor', 'none', 'facealpha', .2);
                %patch([x(:)', x(end) x(1)], [h(:)' 0 0], 1, 'facecolor', cmap(i, :), 'edgecolor', 'none', 'facealpha', .2);
                %Plot.LineBar(h);
            end
            xlim([X(1), X(end)]);
            Fig.Fix
            Fig.Hoff
        end

        function Area(varargin)
            if ~ishold
                cla reset
            end
            offset = 2;
            x = [];
            if nargin >= 2
                if isnumeric(varargin{2})
                    x = varargin{1};
                    y = varargin{2};
                    offset = 3;
                else
                    y = varargin{1};
                end
            else
                y = varargin{1};
            end
            %%
            p = inputParser;
            p.KeepUnmatched = true;
            addOptional(p, 'Offset', 0); 
            p.parse(varargin{offset:end});
            opt = p.Results;
            style = Q.struct2cell(p.Unmatched);
            %%
            if isempty(x)
                patch([1:length(y), length(y), 1], [y 0 0] + opt.Offset, 1, style{:});
                %patch([1:length(y), length(y)], [y 0], 1, style{:});
            else
                patch([x, x(end), x(1)], [y 0 0] + opt.Offset, 1, style{:});
                %patch([x x(end)], [y 0], 1, style{:});
            end
        end
        
        function LineBar(varargin)
            % similar to bar only using lines
            ish = ishold;
            if ~ish
                cla reset
            end

            style = {};
            for i=1:length(varargin)
                if ischar(varargin{i})
                    style = varargin(i:end);
                    varargin = varargin(1:i-1);
                    break;
                end
            end
%             if mod(length(varargin), 2) == 0
%                 for i=1:2:length(varargin)
%                     x = varargin{i};
%                     y = varargin{i+1};
%                 end
%             end
            if length(varargin) > 1
                x = varargin{1};
                y = varargin{2};
                if isvector(y) == 1
                    y = y(:);
                end
                idx = floor((0:2*size(y, 1)-1) / 2) + 1;
                fidx = repmat(idx(:), 1, size(y, 2));
                Y = y(fidx);
                if any(diff(x) < 0)
                    error;
                end
                dx = diff(x);
                dxr = dx; dxr(end+1) = dxr(end);
                dxl = dx; dxl = [dxl(1) dxl];
                X = sort([x + dxr/2 x - dxl /2]);
            else
                y = varargin{1};
                if isvector(y) == 1
                    y = y(:);
                end
                idx = floor((0:2*size(y, 1)-1) / 2) + 1;
                fidx = repmat(idx(:), 1, size(y, 2));
                Y = y(fidx);
                X = repmat(floor((1:2*size(y, 1))' / 2) + .5, 1, size(y, 2));
            end
            if true
                patch([X(:)' X(end) X(1) X(1)], [Y(:)' 0, 0 Y(1)], 'b', style{:});
            else
                plot(X, Y, style{:});
            end
            if ~ish
                Fig.XAxis(Q.minnd(X), Q.maxnd(X));
                %set(gca, 'XTick', 1:min(diff(get(gca, 'XTick'))):max(get(gca, 'XLim')));
            end
            
        end

        function LineBar2(varargin)
            % similar to bar only using lines
            ish = ishold;
            if ~ish
                cla reset
            end

            style = {};
            for i=1:length(varargin)
                if ischar(varargin{i})
                    style = varargin(i:end);
                    varargin = varargin(1:i-1);
                    break;
                end
            end
%             if mod(length(varargin), 2) == 0
%                 for i=1:2:length(varargin)
%                     x = varargin{i};
%                     y = varargin{i+1};
%                 end
%             end
            y = varargin{1};
            if isvector(y) == 1
                y = y(:);
            end
            idx = floor((0:2*size(y, 1)-1) / 2) + 1;
            fidx = repmat(idx(:), 1, size(y, 2));
            Y = y(fidx);
            X = repmat(floor((1:2*size(y, 1))' / 2) + .5, 1, size(y, 2));
            if true
                patch([X(:)' X(end) X(1) X(1)], [Y(:)' 0, 0 Y(1)], 'b', style{:});
            else
                plot(X, Y, style{:});
            end
            if ~ish
                Fig.XAxis(Q.minnd(X), Q.maxnd(X));
                set(gca, 'XTick', 1:min(diff(get(gca, 'XTick'))):max(get(gca, 'XLim')));
            end
            
        end
        
        function XLess(varargin)
            style = {};
            out = {};
            for i=1:length(varargin)
                if isnumeric(varargin{i})
                    curr = varargin{i};
                    out{end+1} = Q.index(curr);
                    out{end+1} = curr;
                else
                    style = {varargin(i:end)};
                    break;
                end
            end
            plot(out{:}, style{:});
        end
        
        function Equality(X, Y, varargin)
            % plot a scatter graph assuming that x and y should be equal
            %style = {'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', Colormaps.Red, 'MarkerSize', 2};
            if isempty(varargin) 
                varargin = {10, Colormaps.Red};
            end
            scatter(X, Y, varargin{:});
            ish = ishold;
            Fig.Hon
            a = axis;
            mn = min(a(1), a(3));
            mx = max(a(2), a(4));
            Fig.XAxis(mn, mx);
            Fig.YAxis(mn, mx);
            plot([mn mx], [mn mx], ':', 'color', ones(1,3)*.2, 'HandleVisibility','off');
            axis square
            if ~ish
                Fig.Hoff
            end
            Fig.Fix
        end
        
        function HexaPlot(bins, range, varargin)
            p = inputParser;
            addOptional(p, 'Colormap', Colormaps.FromTo(Colors.White, Colors.PrettyRed, 256)); 
            %addOptional(p, 'Colormap', Colormaps.UglyBlueWhiteRed); 
            %addOptional(p, 'Colormap', Colormaps.FromTo(Colors.White, Colors.Red, 256)); 
            %addOptional(p, 'Colormap', flip(hot(256))); 
            addOptional(p, 'ShowBar', false); 
            addOptional(p, 'Range', []); 
            addOptional(p, 'MaxColor', [0 0 0]); 
            p.parse(varargin{:});
            opt = p.Results;
            
            %%
            gap = 0.001;
            
            orig = bins;
            bins(isnan(orig)) = 0;
            cla
            mx1 = range(1);
            mx2 = range(2);
            my1 = range(3);
            my2 = range(4);
            dxs = (mx2 - mx1) / (size(bins, 2) * 3 - 1);
            dx = dxs * 3;
            dys = (my2 - my1) / (size(bins, 1) * 2 - 1);
            dy = dys * 2;
            %%
            % xoff = dx / 4;
            % yoff = dy / 2;
            cmap = opt.Colormap;
            tau = mean(bins(:));
            bins(bins > tau * 3) = tau * 3;
            %bins = log(bins);
            
            if isempty(opt.Range)
                m1 = min(bins(isfinite(bins)));
                m2 = max(bins(isfinite(bins)));
            else
                m1 = opt.Range(1);
                m2 = opt.Range(2);
                
            end
            
%            m1 = 0;
%            m2 = ceil(max(bins(isfinite(bins))) / sum(bins(:)) * 1000) / 1000 * sum(bins(:));
            
            if opt.ShowBar;
                colormap(cmap);
                h = colorbar('EastOutside');
                set(h, 'YTick', [0 0.5 1])
                ticks = (get(h, 'YTick') - 1) / (size(cmap, 1) - 1) * (m2 - m1) / sum(bins(:)) + m1 / sum(bins(:));
                labels = {};
                for i=1:length(ticks)
                    labels{i} = sprintf('%.3f', ticks(i));
                end
                set(h, 'YTickLabel', labels);
            end
            
            for i=1:size(bins, 2)
                for j=1:size(bins, 1)
                    x = [(i - 1) * dx - dxs, (i - 1) * dx, i * dx - dxs, i * dx, i * dx - dxs, (i - 1) * dx, (i - 1) * dx - dxs];
                    b = mod(i, 2) - 1;
                    y = [(j - 1) * dy + dys + b * dys, (j - 1) * dy + b * dys, (j - 1) * dy + b * dys, (j - 1) * dy + dys + b * dys, (j - 1) * dy + 2 * dys + b * dys, (j - 1) * dy + 2 * dys + b * dys, (j - 1) * dy + dys + b * dys];
                    if ~isfinite(orig(j, i))
                        patch(x + mx1, y + my1, 'w', 'EdgeColor', 'none'); hold on;
                    else
                        cidx = floor((bins(j, i) - m1) / (m2 - m1) * (size(cmap, 1) - 1)) + 1;
                        cidx = max(cidx, 1);
                        if cidx > size(cmap, 1)
                            color = opt.MaxColor;
                        else
                            color = cmap(cidx, :);
                        end
                        
                        patch(x + mx1, y + my1, color, 'EdgeColor', 'none'); hold on;
                        if gap > 0
                            patch(x + mx1, y + my1, color, 'FaceColor', 'none', 'EdgeColor', [1 1 1]); hold on;
                        end
                    end
                end
            end
            hold off;
            axis([mx1 - dxs, mx2 + dxs, my1 - dys, my2 + dys]);
            axis off;
            axis equal
        end

        function [bins, range] = HexaHist(data, nbins, varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addOptional(p, 'Bounds', []); 
            p.parse(varargin{:});
            opt = p.Results;
            
            %%
            x = data(:, 1);
            y = data(:, 2);
            %%
            if isempty(opt.Bounds)
                mx1 = min(x);
                mx2 = max(x);
                my1 = min(y);
                my2 = max(y);
            else
                range = opt.Bounds;
                mx1 = range(1);
                mx2 = range(2);
                my1 = range(3);
                my2 = range(4);
            end
            
            %%
            %mx1 = 0;
            %mx2 = 1;
            dxs = (mx2 - mx1) / (nbins(2) * 3 - 1);
            %%
            % my1 = 0;
            % my2 = 20;
            dys = (my2 - my1) / (nbins(1) * 2 - 1);
            %%
            bx = floor((x - mx1) / dxs) + 1;
            by = floor((y - my1) / dys) + 1;
            
            ix = x - mx1 - (bx - 1) * dxs;
            iy = y - my1 - (by - 1) * dys;
            
            underSlash = (ix / dxs) > (iy / dys);
            underBSlash = (ix / dxs) < ((dys - iy) / dys);
            
            %%
            border1 = mod(bx, 3) == 0 & mod(bx, 6)/3 == mod(by, 2);
            bx(border1) = bx(border1) + (2 * underSlash(border1) - 1);
            
            border2 = mod(bx, 3) == 0 & mod(bx, 6)/3 == 1-mod(by, 2);
            bx(border2) = bx(border2) - (2 * underBSlash(border2) - 1);
            
            %%
            binx = floor(bx / 3) + 1;
            biny = floor((by - mod(binx, 2)) / 2) + 1;
            %%
            bins = zeros(nbins(1), nbins(2));
            biny(biny > nbins(1)) = nbins(1);
            binx(binx > nbins(2)) = nbins(2);
            binx(binx < 1) = 1;
            biny(biny < 1) = 1;
            %idx = sub2ind([count count], biny(biny <= count & binx <= count), binx(biny <= count & binx <= count));
            idx = sub2ind(nbins, biny, binx);
            bins(1:length(bins(:))) = histc(idx, 1:length(bins(:)));
            range = [mx1, mx2, my1, my2];
            %%
            if nargout == 0
                Plot.HexaPlot(bins, range, p.Unmatched);
            end
        end
        
        function Labeled(varargin)
            h = plot(varargin{:});
            x = get(h, 'XData');
            y = get(h, 'YData');
            for i=1:length(x)
                text(x(i), y(i), sprintf('(%g,%g)\n', x(i), y(i)), 'HorizontalAlignment', 'center', 'verticalalignment', 'bottom');
            end
        end
        
        function Stat(data)
            nbins = max(min(2 * round(sqrt(length(data))), 100), 10);
            [h,x] = hist(data, nbins);
            il = Fig.IgnoreInLegend;
            area(x, h, il{:});
            Fig.Hon
            m = mean(data);
            Fig.VLine(m, 'linestyle', '-', 'Color', 'k');
            md = median(data);
            Fig.VLine(md, 'linestyle', ':', 'Color', 'g');
            s = std(data);
            Fig.VLine(m + s, 'linestyle', ':', 'Color', 'k');
            sd = Q.stdR(data);
            Fig.VLine(md + sd, 'linestyle', ':', 'Color', 'g');
            Fig.Hoff
            Fig.Legend({sprintf('mean (%g)', m), sprintf('median (%g)', md), sprintf('std (%g)', s), sprintf('stdR (%g)', sd)});
        end
        
        function ImagePatch(mat, cmap, range, iscolorbar, istext)
            valid = ~isnan(mat);
            mat(~valid) = min(mat(:));
            if nargin < 5
                istext = false;
            end
            if nargin < 4
                iscolorbar = true;
            end
            if nargin < 3 || isempty(range)
                m1 = min(mat(:));
                m2 = max(mat(:));
%                 scale = ceil(abs(log10(m2-m1)))/1e2;
%                 range(1) = floor(m1 / scale) * scale;
%                 range(2) = ceil(m2 / scale) * scale;
                range = [m1, m2];
            end
            ind = gray2ind(mat2gray(mat, range), 256) + 1;
            if nargin < 2 || isempty(cmap)
                cmap = Colormaps.Gray(256);
            end
            ind(~valid) = 0;
            im = ind2rgb(ind, [0 0 0; cmap]);
            for y=1:size(mat, 1)
                for x=1:size(mat, 2)
                    Shapes.RectCenter(x,y,1,1,im(y, x, :))
                    if istext & valid(y, x)  
                        Shapes.Text(x, y, sprintf('%.3g', mat(y, x)), 'w');
                    end
                    Fig.Hon
                end
            end
            if size(mat, 1) == size(mat, 2)
                axis square
            end
            Fig.Fix;
            axis tight
            set(gca, 'XMinorTick', 'off', 'YMinorTick', 'off', 'TickLength', [0 0])
            set(gca, 'Ydir', 'reverse');
            Fig.Hoff
            %%
            if iscolorbar
                %%
                colormap(cmap);
                h = colorbar;
                set(h, 'ytick', [0 1]);
                set(h, 'yticklabel', range(1) + get(h, 'ytick') * (range(2) - range(1)));
            end

        end
        
        function ImageSC(mat, cmap, varargin)
            %%
            idx = find(strcmp('range', varargin), 1);
            if ~isempty(idx)
                range = varargin{idx+1};
                varargin(idx:idx+1) = [];
            else
                range = Q.minmax(mat(:));
            end
            %%
            valid = ~isnan(mat);
            mat(~valid) = min(mat(:));
            ind = gray2ind(mat2gray(mat, range), 256) + 1;
            if nargin == 1
                cmap = Colormaps.BlueRed(256);
            end
            ind(~valid) = 0;
            im = ind2rgb(ind, [0 0 0; cmap]);
            imagesc(im);
            if size(mat, 1) == size(mat, 2)
                axis square
            end
            Fig.Fix;
            axis tight
            %%
            return;
            h = colorbar;
            c = strsplit(sprintf('%.3f\n', range(1) + get(h, 'Ticks') * diff(range)), '\n');
            colormap(cmap);
            %h.Limits = Q.minmax(mat(:))
            set(h, 'TickLabels', c(1:end-1))
        end
        
    end
end