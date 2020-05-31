classdef Patches
    % Plot differently shaped object in a figure
    methods (Static)
        function ds = DefaultStyle
            % Default style for patches
            ds = {'EdgeColor', 'none'};
        end
        
        function Text(x,y,varargin)
            % Add text to figure
            text(x, y, varargin{:}, 'HorizontalAlignment', 'center', 'verticalAlignment', 'middle');
        end
        
        function Curve(x,y,w,c,varargin)
            % Plot a curve on figure
            if ~ishold
                cla reset
            end
            %ds = Patches.DefaultStyle;
            %patch([x x(end:-1:1)],[y y(end:-1:1)],c,ds{:}, varargin{:});
%             patch('Faces',1:length(x)+1,'Vertices',[x(:), y(:); [nan nan]], ds{:},'EdgeColor', c, varargin{:})
            for i=1:length(x)
                Patches.Line([x(i-1) x(i)], [y(i-1) y(i)],w,c,varargin{:});
            end
        end
        
        function Polygon(x,y,c,varargin)
            % Plot a polygon on figure
            if ~ishold
                cla reset
            end
            ds = Patches.DefaultStyle;
            patch(x,y,c,ds{:}, varargin{:});
        end
        
        function Arc(x,y,anglefr, angleto, r,c,varargin)
            % Plot an arc on figure
            if ~ishold
                cla reset
            end
            N = 100;
            ds = Patches.DefaultStyle;
            d = linspace(anglefr, angleto, N);
            X = x + r * cos(d);
            Y = y + r * sin(d);
            X = round(X);
            Y = round(Y);
            patch((X), (Y), c, ds{:}, varargin{:}, 'FaceColor', 'none');
            patch([x (X) x], [y (Y) y], c, ds{:}, varargin{:}, 'EdgeColor', 'none');
        end
        
        function Circle(x,y,r,c,varargin)
            % Plot a circle on figure            
            if ~ishold
                cla reset
            end
            if length(r) == 1
                r = [r, r];
            end
            N = 100;
            ds = Patches.DefaultStyle;
            d = linspace(0, 2*pi, N);
            X = x + r(1) * cos(d);
            Y = y + r(2) * sin(d);
            %patch(ceil(X), ceil(Y), c, ds{:}, varargin{:});
            patch((X), (Y), c, ds{:}, varargin{:});
        end
        
        function Rect(x,y,w,h,c,varargin)
            % Plot a rectangle on figure            
            if ~ishold
                cla reset
            end
            ds = Patches.DefaultStyle;
            patch([x x+w x+w x], [y y y+h y+h], c, ds{:}, varargin{:});
        end
        
        function Line(x,y,width,c,varargin)
            % Plot a line on figure 
            if ~ishold
                cla reset
            end
            %%
            a = atan2(y(2)-y(1), x(2)-x(1));
            if ~isempty(varargin) && isnumeric(varargin{1})
                shift = varargin{1};
                varargin = varargin(2:end);
                x = x + shift * sin(a);
                y = y - shift * cos(a);
            end
            %%
            dx = width * sin(a);
            dy = width * cos(a);
            %ds = Patches.DefaultStyle;
            ds = {'EdgeColor', c};
            patch([x(1)+dx  x(2)+dx x(2)-dx x(1)-dx], [y(1)-dy, y(2)-dy, y(2)+dy, y(1)+dy], c, ds{:}, varargin{:});
        end
        
    end
end