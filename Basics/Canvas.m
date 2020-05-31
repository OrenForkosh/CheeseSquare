classdef Canvas
    % CANVAS creates a canvas that can be drawed on
    properties 
        Surface
        Smoothing = 1; % use smoothing?
    end
    
    properties (Dependent)
        Points;
        Lines;
    end
    
    methods
        function obj = Canvas(varargin)
            % CANVAS initiates the object
            %   Canvas(height, width) creates a canvas of size height x
            %       width
            %
            %   Canvas(source) creates a canvas that contains the matrix
            %   source (which can be an image or a matrix of zeros)
            if isscalar(varargin{1})
                obj.Surface = zeros(varargin{:});
            else
                obj.Surface = varargin{1};
            end
        end
        
        function h = Height(obj)
            % height
            h = size(obj.Surface, 1);
        end

        function w = Width(obj)
            % width
            w = size(obj.Surface, 2);
        end
        
        function sz = Size(obj, varargin)
            % size (i.e. height and width)
            sz = size(obj.Surface, varargin{:});
        end
        
        function im = Image(obj, cmap)
            % create Image from Canvas
            if nargin < 2
                cmap = [1 1 1; Colormaps.Categorical];
            end
            if isfloat(obj.Surface)
                im = ind2rgb(obj.Surface + 1, cmap);
            else
                im = ind2rgb(obj.Surface, cmap);
            end
            if obj.Smoothing > 1
                im = imresize(imresize(im, obj.Smoothing), 1/obj.Smoothing);
            end
        end
        
        function Show(obj)
            % Show canvas
            im = obj.Image();
            imagesc(im);
            set(gca, 'YDir', 'normal');
        end
        
        function obj = set.Points(obj, value)
            % add points to canvas (in the form of a polygon). Parameters
            % to matlab's Polygon function should be assigend as cell
            obj = Polygon(obj, value{:});
        end
        
        function obj = set.Lines(obj, value)
            % add lines to canvas (in the form of a polygon). Parameters
            % to matlab's Polygon function should be assigend as cell
            obj = Polygon(obj, value{:});
        end
        
        function obj = Polygon(obj, x, y, c, options, connected)
            % POLYGON draws a polygon on canvas
            %   obj = Polygon(obj, x, y, c, width)
            %       creates a polygon with vertices at (x, y) coordinates
            %       and color 'c' (optional) and line width 'width'
            %       (optional)
            if nargin < 4
                c = 1;
            end
            width = 1;
            if nargin >=5
                if isscalar(options)
                    width = options;
                else
                    if any(strcmpi(options, 'LineWidth'))
                        idx = find(strcmpi(options, 'LineWidth'), 1, 'last');
                        width = options{idx + 1};
                    end
                end
            end
            surf = false(size(obj.Surface));
            for i=2:length(x)
                len = max(abs(x(i) - x(i-1)), abs(y(i) - y(i-1))) + 1;
                X = round(linspace(x(i-1), x(i), len)); X = max(X, 1); X = min(X, size(obj.Surface, 2));
                Y = round(linspace(y(i-1), y(i), len)); Y = max(Y, 1); Y = min(Y, size(obj.Surface, 1));
                idx = sub2ind(size(obj.Surface), round(Y), round(X));
                surf(idx) = true;
            end
            if width > 1
                surf = imdilate(surf, strel('disk', width - 1));
            end
            obj.Surface(surf) = c;
        end
    end
    
    methods (Static)
        function canvas = FromSize(sz)
            % Creates canvas of given size
            cellsz = num2cell(sz);
            canvas = Canvas(cellsz{:});
        end
        
        function im = MergeWithImage(canvas, im, cmap)
            % Merge canvas with image
            if nargin < 3
                cmap = Colormaps.Categorical;
            end

            map = canvas.Surface ~= 0;
            values = canvas.Surface(map);
            for i=1:size(im, 3)
                ch = im(:, :, i);
                ch(map) = cmap(values, i);
                im(:, :, i) = ch;
            end
        end
        
        function im = Checkers(size, res, count)
            % create canvas with checkers pattern
            if length(size) == 1
                size(2) = size(1);
            end
            if isempty(res)
                res = ceil(size / count);
            end
            if length(res) == 1
                res(2) = res(1);
            end
            [x, y] = meshgrid(0:size(1)-1, 0:size(2)-1);
            X = mod(x / res(1), 2) >= 1;
            Y = mod(y / res(2), 2) >= 1;
            im = xor(X, Y);
        end
        
        function im = AddBorder(im, width, color)
            % add border to canvas 
            if nargin == 1
                width = 1;
            end
            if nargin < 3
                color = 1;
            end
            im(1:width, :) = color;
            im(end-width+1:end, :) = color;
            im(:, 1:width) = color;
            im(:, end-width+1:end) = color;
        end
        
        function im = MatToIm(mat, cmap, varargin)
            % convert matrix to image
            valid = ~isnan(mat) & isfinite(mat);
            mat(~valid) = min(mat(:));
            ind = gray2ind(mat2gray(mat, varargin{:}), 256) + 1;
            if nargin == 1
                cmap = gray(256);
            end
            ind(~valid) = 0;
            im = ind2rgb(ind, [0 0 0; cmap]);
        end
        
    end
end
