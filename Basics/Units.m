classdef Units
    % Units Standard physical units
    methods (Static)
        function SoS = SpeedOfSound(temperature)
            % speed of sound in m/sec
            if nargin < 1
                temperature = 20;
            end
            SoS = (331.3 + 0.606 * temperature);
        end
    end
end
        