classdef DateTime
    % DateTime Tools to work with date and time information
    methods (Static = true)
        function d = Now()
            % DateTime now
            d = now;
        end
        
        function res = AddSecs(curr, dt)
            % Add specific number of seconds to object
            res = curr + DateTime.FromSecs(dt);
        end
        
        function res = FromSecs(t)
            % DateTime from seconds
            res = t / 60 / 60 / 24;
        end

        function res = FromMins(t)
            % DateTime from mins
            res = t / 60 / 24;
        end

        function res = FromHours(t)
            % DateTime from hours            
            res = t / 24;
        end
        
        function str = ToString(t, varargin)
            % DateTime to string
            if numel(t) > 1
                str = cell(size(t));
                for i=1:numel(t)
                    str{i} = datestr(t(i), varargin{:});
                end
            else
                str = datestr(t, varargin{:});
            end
        end

        function res = ToSecs(t)
            % DateTime to total number of seconds
            res = t * 60 * 60 * 24;
        end

        function res = ToFrames(t, sr)
            % DateTime to total number of frames (sr is the frame rate)
            res = (t) * 60 * 60 * 24 * (sr);
        end

        function s = SecToString(sec, sep)
            % Convert secs to string
            if nargin < 2
                sep = ':';
            end
            s = sprintf('%02d%s%02d%s%02d,%02d', floor(sec / 3600), sep, mod(floor(sec / 60), 60), sep, mod(floor(sec), 60), round((sec - floor(sec)) * 100));
        end
        
    end
end
