function arg = ParseStructArguments(defaults, arg)
% Parse function arguments
f = fieldnames(defaults);
for i=1:length(f)
    if ~isfield(arg, f{i})
        arg.(f{i}) = defaults.(f{i});
    end 
end

