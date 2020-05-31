classdef Shapes
    % Shapes Deprecated! Various shape creation methods
    methods (Static)
        function Text(x, y, str, color, varargin)
            % Place formated text on plot
            text(x, y, str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'color', color, 'FontSize', 14, 'FontName', 'Calibri Light', varargin{:});
        end

        function Rect(x, y, width, height, color, varargin)
            % Creates a filled rectangle
            fill([x x+width x+width x], [y y y+height y+height], color, varargin{:});
        end

        function RectCenter(x, y, width, height, color, varargin)
            % Creates a filled rectangle where the coordinates refer to its
            % center
            Shapes.Rect(x-width/2,y-height/2,width,height, color, varargin{:});
        end
        
        function ChessCenter(x, y, width, height, color1, color2, count, varargin)
            % Create a chess image
            if nargin < 7
                count = 5;
            end
            if nargin < 6
                color2 = 'k';
            end
            sx = width / count;
            sy = height / count;
            idx = 0;
            for i=0:count-1
                for j=0:count-1
                    if mod(idx, 2) == 1
                        Shapes.RectCenter(x + sx * i, y + sy * j, sx, sy, color1);
                    else
                        Shapes.RectCenter(x + sx * i, y + sy * j, sx, sy, color2);
                    end
                    idx = idx + 1;
                    Fig.Hon
                end
            end
            Fig.Hoff;
        end
        
        function Circle(x, y, radios, color, varargin)
            % Draw a circle
            res = 0.02;
            cord = [];
            center = [x, y];
            for angle=0:res:2*pi
                cord = [cord; center + [radios * sin(angle), radios * cos(angle)]];
            end
            fill(cord(:,1), cord(:,2), color, varargin{:});
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
            
            % % % --------------------------------
            % % % Author: Dr. Murtaza Khan
            % % % Email : drkhanmurtaza@gmail.com
            % % % --------------------------------
            
            
        end
    end
end