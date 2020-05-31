classdef Colors
    % Basic colors and the tools to manipuate them
    methods (Static)
        
        %% Basic colors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function c = Black;       c = Colors.FromHex('000000'); end
        function c = White;       c = Colors.FromHex('FFFFFF'); end
        function c = Red;         c = Colors.FromHex('FF0000'); end
        function c = Yellow;      c = Colors.FromHex('FFFF00'); end
        function c = Blue;        c = Colors.FromHex('0000FF'); end
        function c = PrettyRed;   c = [220 073 089] / 255; end
        function c = PrettyGreen; c = [053 178 087] / 255; end
        function c = PrettyBlue;  c = [000 177 229] / 255; end
        function c = DarkGray;    c = [076 076 076] / 255; end
        function c = LightGray;   c = [178 178 178] / 255; end
        
        function c = DutchTeal;   c = Colors.FromHex('1693A5'); end
        function c = HotPink;     c = Colors.FromHex('FF0066'); end
        function c = Serenity;    c = Colors.FromHex('ADD8C7'); end
        function c = Gold;        c = Colors.FromHex('FBB829'); end
        function c = HauntedMilk; c = Colors.FromHex('CDD7B6'); end
        function c = Slate;       c = Colors.FromHex('556270'); end
        function c = Frogs;       c = Colors.FromHex('C3FF68'); end
        function c = Vanilla;     c = Colors.FromHex('FCFBE3'); end
        function c = Bloons;      c = Colors.FromHex('D31996'); end
        function c = VitaminC;    c = Colors.FromHex('FF9900'); end
        
        function c = RetroBlue;   c = Colors.FromHex('80A2CA'); end
        function c = RetroGreen;  c = Colors.FromHex('C3D254'); end
        function c = RetroRed;    c = Colors.FromHex('E45E57'); end
        function c = RetroOrange; c = Colors.FromHex('EF8444'); end
        function c = RetroPurple; c = Colors.FromHex('A07EBB'); end
        function c = RetroYellow; c = Colors.FromHex('E4D859'); end

        %% Color tools %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function c = RGB2CMYK(m)
            % Convert color from RGB to CMYK
            cform = makecform('srgb2cmyk');
            c = reshape(applycform(reshape(m, [size(m,1) 1 3]),cform), [size(m,1) 4]);
        end
        
        function c = CMYK2RGB(m)
            % Convert color from CMYK to RGB
            cform = makecform('cmyk2srgb');
            c = reshape(applycform(reshape(m, [size(m,1) 1 4]),cform), [size(m,1) 3]);
        end
        
        function c = Mix(varargin)
            % Mix a set of colors. Can speficy a 'weights' options to assign 
            % a weight for each color
            if any(strcmpi(varargin, 'weights'))
                stridx = find(strcmpi(varargin, 'weights'));
                weights = varargin{stridx + 1};
                varargin([stridx stridx+1]) = [];
                m = cat(1, varargin{:});
            else
                m = cat(1, varargin{:});
                weights = ones(size(m, 1), 1);
            end
            cmyk = Colors.RGB2CMYK(m);
            c = Colors.CMYK2RGB(min(weights(:)' * cmyk , 1));
        end
        
        function c = Value(c)
            % Color's value (H from HSV)
            c = Colors.RGB2HSV(c);
            c = c(3);
        end
        
        function c = Saturation(c)
            % Color's saturation (S from HSV)
            c = Colors.RGB2HSV(c);
            c = c(2);
        end
        
        function c = CompLight(c)
            % Complementary color in lightness
            c = Colors.RGB2HSV(c);
            c(3) = 1-c(3);
            %c(2) = 1-c(2);
            c = Colors.HSV2RGB(c);
        end
        
        function c = Lighten(c, beta)
            % Make color lighter
            hsv = rgb2hsv(reshape(c, [1 1 3]));
            hsv(2) = hsv(2) * beta;
            c = reshape(hsv2rgb(hsv), [1 3 1]);
        end
        
        function c = Reset(c, v)
            hsv = rgb2hsv(reshape(c, [1 1 3]));
            for i=1:3
                if ~isnan(v(i))
                    hsv(i) = v(i);
                end
            end
            c = reshape(hsv2rgb(hsv), [1 3 1]);
        end
        
        function c = SetLight(c, l)
            % Set color's lightness
            hsv = rgb2hsv(reshape(c, [1 1 3]));
            hsv(2) = l;
            c = reshape(hsv2rgb(hsv), [1 3 1]);
        end
        
        function c = SetSaturation(c, s)
            % Set color's saturation
            hsv = rgb2hsv(reshape(c, [1 1 3]));
            hsv(3) = s;
            c = reshape(hsv2rgb(hsv), [1 3 1]);
        end
        
        function c = Hueten(c, beta)
            if nargin < 2
                beta = .5;
            end
            hsv = rgb2hsv(c);
            hsv(1) = mod(hsv(1) + beta, 1);
            c = hsv2rgb(hsv);
        end
        
        function c = Brighten(c, beta)
            % Brighten a color
            if nargin < 2
                beta = .5;
            end
            c = brighten(c, beta);
        end
        
        function Show(color)
            % Display the color
            Colormaps.Show(color);
        end
        
        function rgb = TSL2RGB(tsl)
            % Convert from TSL (tint,saturation,level) map to RGB
            x = -cot(2*pi*tsl(:, 1));
            g = -sqrt(5./(9*(x.^2+1))).*tsl(:, 2).*(tsl(:, 1) > .5) + ...
                sqrt(5./(9*(x.^2+1))).*tsl(:, 2).*(tsl(:, 1) < .5);
            r = sqrt(5)/3*tsl(:,2).*(tsl(:, 1) == 0) + (x .* g + 1/3).*(tsl(:, 1) ~= 0);
            k = tsl(:, 3) ./ (.185*r + .473*g + .114);
            rgb(:, 1) = k.*r;
            rgb(:, 2) = k.*g;
            rgb(:, 3) = k.*(1 - r - g);
        end
        
        function tsl = RGB2TSL(c)
            % Convert from RGB to TSL (tint,saturation,level) 
            r = bsxfun(@rdivide, c(:, 1), sum(c, 2)) - 1/3;
            g = bsxfun(@rdivide, c(:, 2), sum(c, 2)) - 1/3;
            t = 1/(2*pi)*atan(r./g + .25) .* (g > 0) + ...
                1/(2*pi)*atan(r./g + .75) .* (g < 0);
            s = sqrt(9/5*(r.^2+g.^2));
            l = .299*c(:, 1)+.587*c(:, 2)+.114*c(:, 3);
            tsl = [t,s,l];
        end
        
        function t = GetTint(c)
            % Returns the tint of a color
            r = bsxfun(@rdivide, c(:, 1), sum(c, 2)) - 1/3;
            g = bsxfun(@rdivide, c(:, 2), sum(c, 2)) - 1/3;
            t = 1/(2*pi)*atan(r./g + (g > 0).*.25 + (g<1).*.75);
        end
        
        function t = SetTint(c, t)
            % Sets the tint of a color
        end
        
        function rgb = RGB2HSV(c)
            % Converts RGB color to HSV
            if size(c, 3) == 3
                rgb = rgb2hsv(c);
            else
                rgb = reshape(rgb2hsv(reshape(c, [size(c,1) 1 3])), [size(c,1) 3]);
            end
        end
        
        function hsv = HSV2RGB(c)
            % Convert HSV color to RGB
            if size(c, 3) == 3
                hsv = hsv2rgb(c);
            else
                hsv = reshape(hsv2rgb(reshape(c, [size(c,1) 1 3])), [size(c,1) 3]);
            end
        end
        
        function color = ParseHex(varargin)
            % Color from hex string
            color = Colors.FromHex(varargin{:});
        end

        function color = FromHex(varargin)
            % Color from hex string
            color = zeros(length(varargin), 3);
            for i=1:length(varargin)
                h = varargin{i};
                color(i, :) = [hex2dec(h(1:2)) hex2dec(h(3:4)) hex2dec(h(5:6))] / 255;
            end
        end
        
    end
end