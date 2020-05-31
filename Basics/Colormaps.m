classdef Colormaps
    % Many colormaps and the tools to manipulate them
    
    methods (Static)
        % Basic colors:
        function c = Red;   c = [220 073 089] / 255; end
        function c = Green; c = [053 178 087] / 255; end
        function c = Blue;  c = [000 177 229] / 255; end
    end
    
    methods (Static)
        function z = BaseColors
            % Basic colors
            z = {...
                'Purple', [202   152   227]/255, ...
                'Red',    [196    77    88]/255, ...
                'Orange', [253   189    51]/255, ...
                'Yellow', [253   189    51]/255, ...
                'Green',  [179   224   110]/255, ...
                'Blue',   [000   177   229]/255, ...
                'White',  [236   229   206]/255, ...
                'KBlack', [236   229   206]/255 ...
                };
        end
        function c = Generate(colors)
            % Generate color map from specified colors. For example,
            % Generate('PRBY') creates a colormap with 4 colors: purple,
            % red, blue and yellow
            z = Colormaps.BaseColors;
            c = zeros(length(colors), 3);
            for i=1:length(colors)
                idx = find(cellfun(@(x) ~isempty(x), regexp(z(1:2:end), colors(i))), 1);
                c(i, :) = z{2 * (idx-1) + 2};
            end
        end        
        
        %% Basic colormaps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function cmap = BlueRed(nlevels) % Gradient
            % Continous shift from blue to red
            if nargin == 0
                nlevels = 256;
            end
            cmap = Colormaps.makeColorMap([17 113 185]/255,[170 55 73]/255,nlevels);
        end
        
        function cmap = BlueWhiteRed(nlevels) % Gradient
            % Continous shift from blue to white and to red
            if nargin == 0
                nlevels = 256;
            end
            cmap = Colormaps.makeColorMap([17 113 185]/255,[255 255 255]/255,[170 55 73]/255,nlevels);
        end

        function cmap = BlueGrayRed(nlevels) % Gradient
            % Continous shift from blue to gray and to red
            if nargin == 0
                nlevels = 256;
            end
            cmap = Colormaps.makeColorMap([17 113 185]/255,[255 255 255]/255*.8,[170 55 73]/255,nlevels);
        end
        
        function cmap = UglyBlueWhiteRed(nlevels) % Gradient
            % Continous shift from blue to white and to red. An ugly but a 
            % noticeable colormap
            s = 4;
            if nargin == 0
                nlevels = 256;
            end
            mid = floor(nlevels/2);
            l = exp(linspace(0, 1, mid)' * s) / exp(s);
            fl = flip(l, 1);
            cmap = [[1-fl 1-fl fl*0+1];[fl*0+1 1-l 1-l]];
        end
        
        function cmap = UglyBlueRed(nlevels) % Gradient
            % Continous shift from blue to red. An ugly but a noticeable
            % colormap
            if nargin == 0
                nlevels = 256;
            end
            cmap = Colormaps.makeColorMap([0 0 1],[1 0 0],nlevels);
        end
        
        function cmap = Gray(nlevels)
            % Gray colormap
            if nargin == 0
                nlevels = 256;
            end
            cmap = gray(nlevels);
        end

        function cmap = VentanaAzul(nlevels) % Gradient, Spectrum
            % Orange to white to blue gradient
            if nargin == 0
                nlevels = 256;
            end
            %%
            base = {'F2385A','F5A503','E9F1DF','4AD9D9','36B1BF'};
            cmap = Colormaps.Make(nlevels, base{:});
            if nargout == 0
                Colormaps.Show(cmap);
            end
        end        
        
        function cmap = CheerUpEmoKid(nlevels) % Gradient, Spectrum
            if nargin == 0
                nlevels = 256;
            end
            %%
            base = {'556270','4ECDC4','C7F464','FF6B6B','C44D58'};
            cmap = Colormaps.Make(nlevels, base{:});
            if nargout == 0
                Colormaps.Show(cmap);
            end
        end        
        
        function cmap = OuterRings(nlevels) % Gradient
            if nargin == 0
                nlevels = 256;
            end
            %%
            base = {'EFF3CD', 'B2D5BA', '61ADA0', '248F8D', '605063'};
            cmap = Colormaps.Make(nlevels, base{:});
            if nargout == 0
                Colormaps.Show(cmap);
            end
        end

        function cmap = FreshCutDay(nlevels) % Gradient
            if nargin == 0
                nlevels = 256;
            end
            %%
            base = {'00A8C6', '40C0CB', 'F9F2E7', 'AEE239', '8FBE00'};
            cmap = Colormaps.Make(nlevels, base{:});
            if nargout == 0
                Colormaps.Show(cmap);
            end
        end

        function cmap = GiantGoldfish(nlevels) % Gradient
            if nargin == 0
                nlevels = 256;
            end
            %%
            base = {'4EB3DE', 'A7DBD8', 'E0E4CC', 'F38630', 'FA6900'};
            cmap = Colormaps.Make(nlevels, base{:});
            if nargout == 0
                Colormaps.Show(cmap);
            end
        end

        function cmap = Borel(nlevels) % Gradient
            if nargin == 0
                nlevels = 256;
            end
            %%
            base = {'FFF5DE', 'B8D9C8', '917081', '750E49', '4D002B'};
            cmap = Colormaps.Make(nlevels, base{:});
            if nargout == 0
                Colormaps.Show(cmap);
            end
        end
        
        function cmap = OrigamiLuckyStars(nlevels) % Gradient
            if nargin == 0
                nlevels = 256;
            end
            %%
            base = {'4EB3DE', '8DE0A6', 'FCF09F', 'F27C7C', 'DE528C'};
            cmap = Colormaps.Make(nlevels, base{:});
            if nargout == 0
                Colormaps.Show(cmap);
            end
        end

        function cmap = Retro(ncolors) % Discrete
            cmap(1, :) = Colors.ParseHex('80A2CA');
            cmap(2, :) = Colors.ParseHex('C3D254');
            cmap(3, :) = Colors.ParseHex('E45E57');
            cmap(4, :) = Colors.ParseHex('EF8444');
            cmap(5, :) = Colors.ParseHex('A07EBB');
            cmap(6, :) = Colors.ParseHex('E4D859');
            if ~exist('ncolors', 'var')
                ncolors = size(cmap, 1);
            end
            if ncolors > size(cmap, 1)
                cmap = repmat(cmap, ceil(ncolors/size(cmap, 1)), 1);
            end
            cmap = cmap(1:ncolors, :);
        end
        
        function cmap = Categorical(ncolors) % Discrete
            cmap = Colormaps.SubCategorical;
            cmap = cmap(2:2:end, :);
            if nargin < 1
                ncolors = 256;
            end
            cmap = repmat(cmap, ceil(ncolors / size(cmap, 1)), 1);
            cmap = cmap(1:ncolors, :);
        end
        
        function cmap = SubCategorical % Discrete
            cmap = ...
                [148 216 239;
                0 177 229;
                179 213 157;
                53 178 87;
                234 165 193;
                220 73 89;
                240 191 148;
                240 145 55;
                189 190 200;
                98 111 179;
                246 237 170;
                248 231 83]/255;
        end
        
        function cmap = Random(nlevels)
            % Random colormap
            stream = RandStream('mt19937ar', 'Seed', 5489);
            %cmap = stream.rand(nlevels, 3);
            cmap = zeros(nlevels, 3);
            mix = [1 1 1];
            for i=1:nlevels
                cmap(i, :) = stream.rand(1,3);
                cmap(i, :) = (cmap(i, :) + mix) / 2;
            end
        end
        
        %% Colormap manipulation tools %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function cmap = Make(nlevels, varargin)
            % Auxilary tool to create a colormap
            c = zeros(length(varargin), 3);
            for i=1:length(varargin)
                if ischar(varargin{i})
                    c(i, :) = Colors.ParseHex(varargin{i});
                else
                    c(i, :) = varargin{i};
                end
            end
            %%
            offset = 1;
            jump = nlevels / (size(c,1) - 1);
            cmap = zeros(nlevels, 3);
            idx = 1;
            for i=1:size(c,1)-1
                if i+1 < size(c,1)
                    sz = floor(round(offset + jump) - offset);
                    C = Colormaps.Interpolate(sz + 1, c(i, :), c(i+1, :));
                else
                    sz = nlevels - idx + 1;
                    C = Colormaps.Interpolate(sz, c(i, :), c(i+1, :));
                end
                cmap(idx:idx+sz-1, :) = C(1:sz, :);
                idx = idx + sz;
                offset = offset + jump;
            end
        end
        
        function C = Interpolate(nlevels, c1, c2)
            C = zeros(nlevels, 3);
            for ch=1:3
                C(:, ch) = linspace(c1(ch), c2(ch), nlevels);
            end
        end
    end
        
    methods (Static)
        function Show(cmap)
            % Display a colormap
            im = reshape(cmap, [size(cmap, 1), 1, 3]);
            imagesc(im);
            Fig.XTick([]);
        end
        
        function ShowOnMap(cmap)
            % Display the colormap as points on the colorspace
            im = reshape(cmap, [size(cmap, 1), 1, 3]);
            hsvim = rgb2hsv(im);
            [x,y] = meshgrid(linspace(0, 1, 256), linspace(0, 1, 256));
            huesat = hsv2rgb(cat(3, x, y, ones(size(x))));
            imagesc(huesat)
            hold on;
            axis square;
            axis off
            for i=1:size(cmap, 1)
                c = cmap(i, :) * 255;
                plot(hsvim(i, 1)*255, hsvim(i, 2) * 255, 'o', 'MarkerFaceColor', cmap(i, :), 'MarkerEdgeColor', 'k');
            end
            hold off;
        end
        
        function im = MatToIm(mat, cmap, varargin)
            % Turn a colormap to a displayable image
            valid = ~isnan(mat) & isfinite(mat);
            mat(~valid) = min(mat(:));
            ind = gray2ind(mat2gray(mat, varargin{:}), 256) + 1;
            if nargin == 1
                cmap = Colormaps.BlueRed(256);
            end
            ind(~valid) = 0;
            im = ind2rgb(ind, [0 0 0; cmap]);
        end
        
        function cmap = FromTo(varargin)
            % Create colormap from colors
            if isscalar(varargin{end})
                nlevels = varargin{end};
                varargin(end) = [];
            else
                nlevels = 256;
            end
            cmap = [];
            pos = nlevels / max(length(varargin) - 1, 1);
            for i=2:length(varargin)
                src = varargin{i-1};
                dst = varargin{i};
                c = Colormaps.makeColorMap(src,dst,pos*(i-1) - size(cmap, 1));
                cmap = [cmap; c]; %#ok<AGROW>
            end
            
        end
        
        function color = GetColor(idx)
            % Get a specific color by if
            cmap = Colormaps.Categorical;
            color = cmap(mod(idx-1, size(cmap, 1))+1, :);
        end
    end
    
    methods (Static, Access = public)
        function cMap = makeColorMap(varargin)
            %% MAKECOLORMAP makes smoothly varying colormaps
            % a = makeColorMap(beginColor, middleColor, endColor, numSteps);
            % a = makeColorMap(beginColor, endColor, numSteps);
            % a = makeColorMap(beginColor, middleColor, endColor);
            % a = makeColorMap(beginColor, endColor);
            %
            % all colors are specified as RGB triples
            % numSteps is a scalar saying howmany points are in the colormap
            %
            % Examples:
            %
            % peaks;
            % a = makeColorMap([1 0 0],[1 1 1],[0 0 1],40);
            % colormap(a)
            % colorbar
            %
            % peaks;
            % a = makeColorMap([1 0 0],[0 0 1],40);
            % colormap(a)
            % colorbar
            %
            % peaks;
            % a = makeColorMap([1 0 0],[1 1 1],[0 0 1]);
            % colormap(a)
            % colorbar
            %
            % peaks;
            % a = makeColorMap([1 0 0],[0 0 1]);
            % colormap(a)
            % colorbar
            
            % Reference:
            % A. Light & P.J. Bartlein, "The End of the Rainbow? Color Schemes for
            % Improved Data Graphics," Eos,Vol. 85, No. 40, 5 October 2004.
            % http://geography.uoregon.edu/datagraphics/EOS/Light&Bartlein_EOS2004.pdf
            
            defaultNum = 100;
            errorMessage = 'See help MAKECOLORMAP for correct input arguments';
            
            if nargin == 2 %endPoints of colormap only
                color.start  = varargin{1};
                color.middle = [];
                color.end    = varargin{2};
                color.num    = defaultNum;
            elseif nargin == 4 %endPoints, midPoint, and N defined
                color.start  = varargin{1};
                color.middle = varargin{2};
                color.end    = varargin{3};
                color.num    = varargin{4};
            elseif nargin == 3 %endPoints and num OR endpoints and Mid
                if numel(varargin{3}) == 3 %color
                    color.start  = varargin{1};
                    color.middle = varargin{2};
                    color.end    = varargin{3};
                    color.num    = defaultNum;
                elseif numel(varargin{3}) == 1 %numPoints
                    color.start  = varargin{1};
                    color.middle = [];
                    color.end    = varargin{2};
                    color.num    = varargin{3};
                else
                    error(errorMessage)
                end
            else
                error(errorMessage)
            end
            
            if color.num <= 1
                error(errorMessage)
            end
            
            if isempty(color.middle) %no midPoint
                cMap = Colormaps.interpMap(color.start, color.end, color.num);
            else %midpointDefined
                [topN, botN] = Colormaps.sizePartialMaps(color.num);
                cMapTop = Colormaps.interpMap(color.start, color.middle, topN);
                cMapBot = Colormaps.interpMap(color.middle, color.end, botN);
                cMap = [cMapTop(1:end-1,:); cMapBot];
            end
        end
        
        function cMap = interpMap(colorStart, colorEnd, n)
            for i = 1:3
                cMap(1:n,i) = linspace(colorStart(i), colorEnd(i), n);
            end
        end
        
        function [topN, botN] = sizePartialMaps(n)
            n = n + 1;
            topN =  ceil(n/2);
            botN = floor(n/2);
        end
    end
end