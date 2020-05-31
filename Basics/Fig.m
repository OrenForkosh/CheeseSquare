classdef Fig
    % Fig Set of tools to manipulate figures
    properties (Constant)
        Margin = .025;
        Padding = .1;
    end
        
    methods (Static = true)
        function Default
            % Position figures in a default location of the screen (to get
            % consistency when saving)
            screensize = get( groot, 'Screensize');
            h = screensize(4) * .85;
            set(gcf, 'Position', [20 (screensize(4)-h)/2 h/sqrt(2) h], 'Renderer', 'painters');
        end
        
        function S = MarkerSymbols()
            % List of marker symbols that can be used with the 'plot'
            % command
            s = {'o', '+', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
            S = {};
            while length(S) < 256
                S = [S, s];
            end
        end
        
        function SendBack(h)
            % send handle to the back (beyond other objects)
            c = get(gca, 'Children');
            idx = find(c == h);
            if ~isempty(idx)
                set(gca, 'Children', c([idx, Q.exclude(length(c), idx)]));
            end
        end
        
        function Hon
            % 'hold on' for lazy people
            hold on;
        end

        function Hoff
            % 'hold off' for lazy people
            hold off;
        end
        
        function Equal()
            % make the axis equal (old and stupid. don't use)
            axis square;
            x = xlim;
            y = ylim;
            z = [min(x(1), y(1)) max(x(2), y(2))];
            axis([z, z]);
            Fig.Hon
            plot([z(1) z(end)], [z(1) z(end)], 'k:');
            Fig.Hoff
        end

        function Text(str, varargin)
            % Text add text to a specific location of the figure. 
            %
            %   Text(str, options) Write the text 'str' on the figure.
            %   options include:
            %       location - location to plot text. Can be 'nw', 'sw', 
            %       'ne', 'se', 'n', 'w', 's', 'e' (correspond to
            %       nothe-west, south-west, etc.). The default is 'nw'
            %       margin - distance of the ploted text from the boundries
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('margin', 40);
            p.addOptional('location', 'nw', @(x) ismember(x, {'nw', 'sw', 'ne', 'se', 'n', 'w', 's', 'e'}));
            p.parse(varargin{:});
            opt = p.Results;
            %%
            nf = length(fieldnames(p.Unmatched));
            aux(2:2:nf*2) = struct2cell(p.Unmatched);
            aux(1:2:nf*2) = fieldnames(p.Unmatched);
            
            %%
            x = xlim;
            y = ylim;
            d = Fig.DefaultAxesFont;
            switch opt.location
                case 'n'
                    pos = [(x(2)-x(1))/opt.margin, y(2)-(y(2)-y(1))/opt.margin];
                    text(pos(1), pos(2), str, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', d{:}, aux{:});
                case 's'
                    pos = [(x(2)-x(1))/opt.margin,  y(1)+(y(2)-y(1))/opt.margin];
                    text(pos(1), pos(2), ['  ' str], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', d{:}, aux{:});
                case 'w'
                    pos = [x(1)+(x(2)-x(1))/opt.margin,  -2*(y(2)-y(1))/opt.margin];
                    text(pos(1), pos(2), ['  ' str], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', d{:}, aux{:});
                case 'e'
                    pos = [x(2)-(x(2)-x(1))/opt.margin, -2*(y(2)-y(1))/opt.margin];
                    %pos = [x(2)-(x(2)-x(1))/opt.margin,  -(y(2)-y(1))/opt.margin];
                    text(pos(1), pos(2), ['  ' str], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', d{:}, aux{:});
                case 'nw'
                    pos = [x(1)+(x(2)-x(1))/opt.margin, y(2)-(y(2)-y(1))/opt.margin];
                    text(pos(1), pos(2), str, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', d{:}, aux{:});
                case 'sw'
                    pos = [x(1)+(x(2)-x(1))/opt.margin, y(1)+(y(2)-y(1))/opt.margin];
                    text(pos(1), pos(2), str, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', d{:}, aux{:});
                case 'ne'
                    pos = [x(2)-(x(2)-x(1))/opt.margin, y(2)-(y(2)-y(1))/opt.margin];
                    text(pos(1), pos(2), str, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', d{:}, aux{:});
                case 'se'
                    pos = [x(2)-(x(2)-x(1))/opt.margin, y(1)+(y(2)-y(1))/opt.margin];
                    text(pos(1), pos(2), str, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', d{:}, aux{:});
            end
        end
        
        function TextNW(varargin)
            % plot text to north-west corner of the figure
            Fig.Text(varargin{:}, 'location', 'nw');
        end
        

        function TextSE(varargin)
            % plot text to south-east corner of the figure
            Fig.Text(varargin{:}, 'location', 'se');
        end

        function TextSW(varargin)
            % plot text to south-west corner of the figure
            Fig.Text(varargin{:}, 'location', 'sw');
        end
        
        function AddLayer(alpha)
            % make all plots in the figure brighter by adding a
            % semi-transparent layer
            if nargin < 1
                alpha = .5;
            end
            H = findall(gca, '-property', 'xdata');
            for i=1:length(H)
                h = H(i);
                %%
                meta = get(h);
                names = fieldnames(meta);
                %%
                props = ~cellfun(@isempty, regexpi(names, 'color'));
                for k=find(props(:))'
                    v = meta.(names{k});
                    if isnumeric(v) && length(v) == 3
                        color = v;
                        set(h, names{k}, color + (1-color) * alpha);
                    end
                end
                %%
                if false
                props = ~cellfun(@isempty, regexpi(names, 'FaceVertexCData'));
                for k=find(props(:))'
                    v = meta.(names{k});
                    if isnumeric(v) && size(v, ndims(v)) == 1
                        if size(v, 2) == 1
                            v = squeeze(ind2rgb(v, colormap));
                        else
                            v = ind2rgb(v, colormap);
                        end
                    end
                    if isnumeric(v) && size(v, ndims(v)) == 3
                        v = v + (1 - v) * alpha;
                        set(h, names{k}, v);
                    end
                end
                end
            end
            Fig.Hon
        end
        
        function Fix
            % make figures "prettier" using this fixed set of parameters
            set(gca, ...
                'Box'         , 'off'     , ...
                'TickDir'     , 'out'     , ...
                'TickLength'  , [.005 .005] , ...
                'YMinorTick'  , 'off'      , ...
                'XColor'      , [.3 .3 .3], ...
                'YColor'      , [.3 .3 .3], ...
                'XMinorTick'  , 'on'      , ...
                'LineWidth'   , 1         );
            
            set(gca, 'FontName', 'Calibri Light', 'FontSize', 12);
            set(gcf, 'color', [1 1 1]);            
        end

        function FixFonts
            % make figures "prettier" using this fixed set of font
            % parameters
            set(gca, 'FontName', 'Calibri Light', 'FontSize', 12);
        end

        function res = DefaultAxesFont
            % make figures "prettier" using this fixed set of axes font
            % parameters
            res = {'FontName', 'Calibri Light', 'FontSize', 12,  'color', [.3 .3 .3]};
        end
        
        function res = DefaultLabelFont
            % make figures "prettier" using this fixed set of label font
            % parameters
            res = {'FontName', 'Calibri Light', 'FontSize', 14};
        end
        
        function Legend(t, varargin)
            % show the Legend on the plot. Extends Matlab's default tool by
            % supporting legends which are numerical and by keeping things
            % 'prettier'
            if isnumeric(t)
                entries = cell(1, length(t));
                for i=1:length(t)
                    entries{i} = num2str(t(i));
                end
                h = legend(entries, varargin{:});
                legend boxoff;
            else
                h = legend(t{:}, varargin{:});
                legend boxoff;
            end
            fonts = Fig.DefaultLabelFont;
            set(h, fonts{:});
            
        end
        
        function s = IgnoreInLegend()
            % text needed to add to object if it is not to be included in
            % figure legend
            s = {'HandleVisibility','off'};
        end
        
        function x = Axis(varargin)
            % set axis bounds
            if nargin == 4
                Fig.XAxis(varargin{1}, varargin{2});
                Fig.YAxis(varargin{3}, varargin{4});
                return;
            end
            if nargin == 2 && length(varargin{1}) == 2 && length(varargin{2}) == 2
                Fig.XAxis(varargin{1});
                Fig.YAxis(varargin{2});
            end
            error('wrong number of arguments');
        end
        
        function XAxis(x1, x2)
            % set x-axis bounds. If no parameters are assigned fits it to
            % the data
            if nargin == 0
                limits = objbounds(findall(gca));
                if isempty(limits) % no objects in axes with data limits
                    limits = [get(gca,'XLim') get(gca,'YLim') get(gca,'ZLim')];
                end
                set(gca, 'xlim', limits(1:2));
                return;
            end
            a = axis;
            if nargin == 1
                x2 = x1(2);
                x1 = x1(1);
            end
            if isempty(x2); x2 = a(2); end
            if isempty(x1); x1 = a(1); end
            axis([x1 x2 a(3) a(4)]);
        end
        
        function YAxis(y1, y2)
            % set y-axis bounds. If no parameters are assigned fits it to
            % the data            
            if nargin == 0
                limits = objbounds(findall(gca));
                if isempty(limits) % no objects in axes with data limits
                    limits = [get(gca,'XLim') get(gca,'YLim') get(gca,'ZLim')];
                end
                set(gca, 'ylim', limits(3:4));
                return;
            end
            a = axis;
            if nargin == 1
                y2 = y1(2);
                y1 = y1(1);
            end
            if isempty(y2); y2 = a(4); end
            if isempty(y1); y1 = a(3); end
            axis([a(1) a(2) y1 y2]);
        end
        
        function h = HLine(Y, varargin)
            % HLine(Y, options) draws a horizontal line on the current
            % axes. Options can be:
            %   text - string to put next to the line
            %   (other) - parameteres to be passed to Matlab's line cmd
            a = axis;
            t = strcmpi(varargin, 'text');
            if any(t)
                str = [' ' varargin{find(t)+1}];
                varargin([find(t) find(t)+1]) = [];
            else
                str = '';
            end
            if nargin >= 2
                h = zeros(1, length(Y));
                for i=1:length(Y)
                    y = Y(i);
                    h(i) = line([a(1), a(2)], [y, y], varargin{:});
                    if ~isempty(str)
                        text(a(2), y, str, 'HorizontalAlignment', 'left', 'verticalalignment', 'top', 'color', get(h(i), 'color'));
                    end
                end
            else
                h = zeros(1, length(Y));
                for i=1:length(Y)
                    y = Y(i);
                    h(i) = line([a(1), a(2)], [y, y], 'LineStyle', ':', 'Color', 'k');
                    if ~isempty(str)
                        text(a(2), y, str, 'HorizontalAlignment', 'left', 'verticalalignment', 'top', 'color', get(h(i), 'color'));
                    end
                end

            end
        end
        
        function h = VLine(X, varargin)
            % VLine(Y, options) draws a vertical line on the current
            % axes. Options can be:
            %   text - string to put next to the line
            %   (other) - parameteres to be passed to Matlab's line cmd
            a = axis;
            t = strcmpi(varargin, 'text');
            if any(t)
                str = [' ' varargin{find(t)+1}];
                varargin([find(t) find(t)+1]) = [];
            else
                str = '';
            end
            if nargin >= 2
                h = zeros(1, length(X));
                for i=1:length(X)
                    x = X(i);
                    h(i) = line([x, x], [a(3) a(4)], varargin{:}, 'HandleVisibility','off');
                    if ~isempty(str)
                        text(x, a(4), str, 'HorizontalAlignment', 'left', 'verticalalignment', 'top', 'color', get(h(i), 'color'));
                    end
                end
            else
                h = zeros(1, length(X));
                for i=1:length(X)
                    x = X(i);
                    h(i) = line([x, x], [a(3) a(4)], 'LineStyle', ':', 'Color', 'k');
                    if ~isempty(str)
                        text(x, a(4), str, 'HorizontalAlignment', 'left', 'verticalalignment', 'top');
                    end
                end
            end
        end

        function XTick(x, labels)
            % set the xtick labels
            set(gca, 'xtick', x);
            if nargin > 1
                set(gca, 'xticklabel', labels);
            end
        end

        function YTick(y, labels)
            % set the ytick labels
            set(gca, 'ytick', y);
            if nargin > 1
                set(gca, 'yticklabel', labels);
            end
        end
        
        function YFlip(varargin)
            % flip the direction of the y-axis
            if nargin > 0
                if islogical(varargin{1}) || isnumeric(varargin{1})
                    if state
                        set(gca, 'ydir', 'normal')
                    else
                        set(gca, 'ydir', 'reverse')
                    end
                else
                    set(gca, 'ydir', varargin{1})
                end
            else
                switch get(gca, 'ydir')
                    case 'normal'
                        set(gca, 'ydir', 'reverse')
                    case 'reverse'
                        set(gca, 'ydir', 'normal')
                end
            end
        end

        function Labels(x, y, varargin)
            % set x- and y-labels title
            Fig.XLabel(x, varargin{:});
            Fig.YLabel(y, varargin{:});
        end
        
        function XLabel(varargin)
            % set the x-label title
            s = sprintf(varargin{:});
            xlabel(s, 'FontSize', 14, 'FontName', 'Calibri Light');
            set(gca, 'FontName', 'Calibri Light', 'FontSize', 12);
        end

        function YLabel(varargin)
            % set the y-label title
            s = sprintf(varargin{:});
            ylabel(s, 'FontSize', 14, 'FontName', 'Calibri Light');
            set(gca, 'FontName', 'Calibri Light', 'FontSize', 12);
        end
        
        function Suptitle(varargin)
            % main title for the figure (ignoring subplots)
            suptitle(varargin{:});
        end
        
        function Title(varargin)
            % set title to a subplot
            if ~iscell(varargin{1})
                data = {varargin};
            else
                data = varargin{1};
            end
            s = cell(size(data));
            for i=1:length(data)
                if iscell(data{i})
                    curr = regexprep(data{i}{:}, '\\', '\\\\');
                    s{i} = sprintf(curr);
                else
                    curr = regexprep(data{i}, '\\', '\\\\');
                    s{i} = sprintf(curr);
                end
            end
            title(s, 'FontSize', 16, 'FontName', 'Calibri Light');
            set(gca, 'FontName', 'Calibri Light', 'FontSize', 12);
        end
        
        function Square(count, idx)
            % replaces subplot for a plot that tries to place all sub-plots
            % on a square grid. Square(count, idx) where count is the
            % number of plots, and idx of the current plot
            nrows = ceil(sqrt(count));
            ncols = ceil(count / nrows);
            subplot(nrows, ncols, idx);
        end
        
        function Subplot2(rows, cols, idx)
            % Subplot2(rows, cols, idx) with 'rows' rows and 'cols' cols 
            % with index 'idx' 
            row = floor((idx-1)/cols+1);
            col = mod(idx-1, cols)+1;
            sz = (1 - 2*Fig.Margin);
            axes('Position', [...
                Fig.Margin + (col-1)/cols * sz, ...
                Fig.Margin + (rows - row)/rows * sz, ...
                (sz - Fig.Padding)/rows, (sz-Fig.Padding)/cols]);
        end
        
        function Subplot(rows, cols, row, col)
            % Subplot(rows, cols, idx) generate plot with 'rows' rows and 
            % 'cols' cols and plots the current at row 'row' and col 'col'           
            if length(row) ~= length(col)
                if isscalar(row)
                    row = col * 0 + row;
                else
                    col = row * 0 + col;
                end
            end
            idx = sub2ind([cols, rows], col, row);
            subplot(rows, cols, idx);
        end
        
        function high = Capture()
            % take snapshot of the figure
            screen_DPI = get(0,'ScreenPixelsPerInch');
            tempfile = [tempname '.png'];
            self.source_fig = gcf;
            self.K = [4 4];
            current_paperpositionmode = get(self.source_fig,'PaperPositionMode');
            current_inverthardcopy = get(self.source_fig,'InvertHardcopy');
            set(self.source_fig,'PaperPositionMode','auto');
            set(self.source_fig,'InvertHardcopy','off');
            print(self.source_fig,['-r',num2str(screen_DPI*self.K(1))], '-dpng', tempfile);
            set(self.source_fig,'InvertHardcopy',current_inverthardcopy);
            set(self.source_fig,'PaperPositionMode',current_paperpositionmode);
            self.raw_hires = imread(tempfile);
            delete(tempfile);
            %%
            myconv = @imfilter;

            % Subsample hires figure image with standard anti-aliasing using a
            % butterworth filter
            kk = Fig.lpfilter(self.K(2)*3,self.K(2)*0.9,2);
            mm = myconv(ones(size(self.raw_hires(:,:,1))),kk,'same');
            a1 = max(min(myconv(single(self.raw_hires(:,:,1))/(256),kk,'same'),1),0)./mm;
            a2 = max(min(myconv(single(self.raw_hires(:,:,2))/(256),kk,'same'),1),0)./mm;
            a3 = max(min(myconv(single(self.raw_hires(:,:,3))/(256),kk,'same'),1),0)./mm;
            if abs(1-self.K(2)) > 0.001
                raw_lowres = double(cat(3,a1(2:self.K(2):end,2:self.K(2):end),a2(2:self.K(2):end,2:self.K(2):end),a3(2:self.K(2):end,2:self.K(2):end)));
            else
                raw_lowres = self.raw_hires;
            end
            %%
            high = raw_lowres;
        end
        
        %% A simple lowpass filter kernel (Butterworth).
        % sz is the size of the filter
        % subsmp is the downsampling factor to be used later
        % n is the degree of the butterworth filter
        function kk = lpfilter(sz, subsmp, n)
            sz = 2*floor(sz/2)+1; % make sure the size of the filter is odd
            cut_frequency = 0.5 / subsmp;
            range = (-(sz-1)/2:(sz-1)/2)/(sz-1);
            [ii,jj] = ndgrid(range,range);
            rr = sqrt(ii.^2+jj.^2);
            kk = ifftshift(1./(1+(rr./cut_frequency).^(2*n)));
            kk = fftshift(real(ifft2(kk)));
            kk = kk./sum(kk(:));
        end
        
        function SavePNG(filename)
            fn = MyFilename(filename);
            fn = fn.SetExt('.png');
            saveas(gcf, fn.Full);
        end
        
        function SavePDF(filename)
            fn = MyFilename(filename);
            fn = fn.SetExt('.pdf');
            saveas(gcf, fn.Full);
        end
        
        function Save(filename)
            fn = MyFilename(filename);
            fn = fn.SetExt('.png');
            saveas(gcf, fn.Full);
            %%
            fn = fn.SetExt('.eps');
            set(gcf,'PaperPositionMode','auto')
            set(gcf,'InvertHardcopy','off')
            print(gcf, fn.Full, '-dpsc2', '-r300', '-painters');
            % now read in the file
            fid = fopen(fn.Full);
            ff = char(fread(fid))';
            fclose(fid);
            
            %get the actual font
            actualfont = strrep(get(gca,'FontName'), ' ', '-');
            
            %these are the only allowed fonts in MatLab and so we have to weed them out
            %and replace them:
            mlabfontlist = {'AvantGarde','Helvetica-Narrow','Times-Roman','Bookman',...
                'NewCenturySchlbk','ZapfChancery','Courier','Palatino','ZapfDingbats',...
                'Helvetica'};%,'Symbol'};
            
            for k = 1:length(mlabfontlist)
                ff = strrep(ff,mlabfontlist{k},actualfont);
            end
            
            % open the file up and overwrite it
            fid = fopen(fn.Full,'w');
            fprintf(fid,'%s',ff);
            fclose(fid);
            %%
            fn = fn.SetExt('.fig');
            saveas(gcf, fn.Full);
        end
    end
end