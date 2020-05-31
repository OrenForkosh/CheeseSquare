function obj = TrackSave(obj)
atomic = true;
obj = TrackOptimize(obj);

filename = [obj.OutputPath obj.FilePrefix '.obj.mat'];
fprintf(['# (>) saving tracking object to ''%s''\n'], filename);
%save(filename, 'obj');
if atomic
    save([filename '.temp'], '-struct', 'obj', '-v7');
    movefile([filename '.temp'], filename);
else
    save(filename, '-struct', 'obj', '-v7');
end    