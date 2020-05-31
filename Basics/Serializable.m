classdef Serializable
    % Serializable Allow classes to serialize themself to a struct variable
    % Used by inheriting from this class (see 'MyVideoReader', for example)
    methods
        function out = SerializeOut(obj)
            % Serialize out the class. Should be overridden in descendant
            % classes
            out = Serializable.Serialize(obj);
        end
    end
    
    methods (Static = true)
        function out = Serialize(obj)
            % Go over all variables belonging to the class and output them
            % to a struct
            out = struct();
            if ~isscalar(obj)
                %%
                idx = 1:numel(obj);
                for i=idx
                    %%
                    coord = cell(1, ndims(obj));
                    [coord{:}] = ind2sub(size(obj), i);
                    if i == 1
                        out = Serializable.Serialize(obj(i));
                    else
                        out(coord{:}) = Serializable.Serialize(obj(i));
                    end
                end
            else
                f = fieldnames(obj);
                for i=1:length(f)
                    if isobject(obj.(f{i}))
                        if any(strcmp(superclasses(obj.(f{i})), 'Serializable'))
                            out.(f{i}) = obj.(f{i}).SerializeOut;
                        else
                            out.(f{i}) = obj.(f{i});
                        end
                    elseif isstruct(obj.(f{i}))
                        out.(f{i}) = Serializable.Serialize(obj.(f{i}));
                    else
                        out.(f{i}) = obj.(f{i});
                    end
                end
            end
            
        end
    end
end
