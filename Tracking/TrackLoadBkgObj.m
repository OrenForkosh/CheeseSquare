function obj = TrackLoadBkgObj(bkgobj)
% Loads older versions of CheeseSquare objects that ran in the background
% (using Job)
obj = [];
if isfield(bkgobj, 'Type') && strcmp(bkgobj.Type, 'bkgobj')
    %%
    fprintf('# (<) - loading from temp file...\n');
    count = 1;
    first = true;
    running = true;
    while running
        strud = ReadStrudels([bkgobj.SourceObj.OutputPath bkgobj.SourceObj.FilePrefix '.output']);
        if count > 1
            fprintf('\b\b\b\b\b\b');
        end
        if isfield(strud, 'status')
            switch strud.status
                case 'done'
                    if count > 1
                        fprintf('\n');
                    end
                    obj = TrackLoad([bkgobj.SourceObj.OutputPath bkgobj.SourceObj.FilePrefix '.obj.mat']);
                    obj.FilePrefix = bkgobj.FilePrefix;
                    obj.OutputPath = bkgobj.OutputPath;
                    running = false;
                    fprintf('# (<) done!\n');
                case 'failed'
                    if count > 1
                        fprintf('\n');
                    end
                    throw(MException('Tracking:InvalidBackgroundObject', ['background object generation failed for ' bkgobj.FilePrefix ' (' bkgobj.SourceObj.FilePrefix '), due to ' strud.message]));
                    running = false;
                case 'running'
                    running = true;
                otherwise
                    throw(MException('Tracking:InvalidBackgroundObject', ['background object generation failed for ' bkgobj.FilePrefix ' (' bkgobj.SourceObj.FilePrefix '), due to ' strud.message]));
            end
        end
        if running
            if first
                fprintf(['# waiting for background object: ' bkgobj.FilePrefix ' (' bkgobj.SourceObj.FilePrefix ')   ']);
                logid = fopen(bkgobj.Process.logfile);
            end
            tline = fgetl(logid);
            while ischar(tline)
                fprintf('\n# (<)   . %s', tline);
                tline = fgetl(logid);
            end
            switch mod(count-1,4)
                case 0
                    fprintf('[   ] ');
                case 1
                    fprintf('[#  ] ');
                case 2
                    fprintf('[## ] ');
                case 3
                    fprintf('[###] ');
            end
            count = count + 1;
            pause(.5);
        end
        first = false;
    end
end