classdef Debug
    % Auxiliary tools for debugging
    properties
        UseDebug = true;
    end
    
    methods
        function IsValid(obj, x)
            % ensures that no members are nan or infinite
            if obj.UseDebug
                dbstop if warning Debug:IsValid:ValidationFailed
                if any(isnan(x(:))) || any(isinf(x(:)))
                    warning('Debug:Assert:AssertFailed', 'validation failed');
                end
            end
        end
        
        function Assert(obj, v)
            % Assertion
            if obj.UseDebug
                dbstop if warning Debug:Assert:AssertFailed
                if ~v
                    warning('Debug:Assert:AssertFailed', 'assertion failed');
                end
            end
        end
    end
end
