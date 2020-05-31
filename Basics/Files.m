classdef Files
    % Files A set of helper tools to working with files
    methods (Static = true)
        function Write(filename, str)
            % write a single line of text to a file
            fileID = fopen(filename,'w');
            fwrite(fileID,str);
            fclose(fileID);
        end

        function Append(filename, str)
            % append a single line of text to a file
            fileID = fopen(filename,'a');
            fwrite(fileID,str);
            fclose(fileID);
        end
        
        function d = CreationDate(filename)
            % get file's creation date
            l = dir(filename);
            d = l.datenum;
        end
        
        function ExportToH5(filename, data, prefix)
            % ExportToH5(filename, data, prefix) write data into H5 file
            % if 'data' is struct then writes all the fields
            if nargin < 3
                prefix = '/';
            end
            if isa(data, 'struct')
                f = fields(data);
                for i=1:length(f)
                    if isa(data.(f{i}), 'struct')
                        Files.ExportToH5(filename, data, [prefix f{i} '/'])
                    else
                        chunksize = max(ones(1, ndims(data.(f{i}))), min(50 * ones(1, ndims(data.(f{i}))), size(data.(f{i}))));
                        h5create(filename, [prefix f{i}], size(data.(f{i})), 'Deflate', 5, 'ChunkSize', chunksize);
                        %h5create(filename, [prefix f{i}], size(data.(f{i})));
                        h5write(filename, [prefix f{i}], data.(f{i}));
                    end
                    
                end
            end
        end
        
        function PrintStrudel(data, prefix)
            % PrintStrudel(data, prefix) prints data in Strudel format
            % if 'data' is struct then writes all the fields. This is a
            % special format to write matrix data to text files.
            if nargin < 2
                prefix = '/';
            end
            if isa(data, 'struct')
                f = fields(data);
                for i=1:length(f)
                    if isa(data.(f{i}), 'struct')
                        Files.PrintStrudel(data.(f{i}), [prefix f{i} '/'])
                    else
                        curr = data.(f{i});
                        %%
                        sz = size(curr);
                        szstr = [num2str(sz(1)) sprintf(',%d', sz(2:end))];
                        fprintf('@%s[%s(%s)] =', [prefix f{i}], class(data), szstr);
                        switch class(curr)
                            case {'double', 'single'}
                                fprintf(' %g', curr);
                            case {'char'}
                                fprintf(' %s', curr);
                        end
                        fprintf('\n');
                    end
                    
                end
            end
        end
        
        function str = GetStrudel(data, prefix)
            % PrintStrudel(data, prefix) prints data in Strudel format
            % if 'data' is struct then writed all the fields
            if nargin < 2
                prefix = '/';
            end
            if isa(data, 'struct')
                f = fields(data);
                str = '';
                for i=1:length(f)
                    if isa(data.(f{i}), 'struct')
                        str = [str Files.GetStrudel(data.(f{i}), [prefix f{i} '/'])];
                    else
                        curr = data.(f{i});
                        %%
                        sz = size(curr);
                        szstr = [num2str(sz(1)) sprintf(',%d', sz(2:end))];
                        str = [str sprintf('@%s[%s(%s)] =', [prefix f{i}], class(data), szstr)];
                        switch class(curr)
                            case {'double', 'single'}
                                str = [str sprintf(' %g', curr)];
                            case {'char'}
                                str = [str sprintf(' %s', curr)];
                        end
                        str = [str sprintf('\n')];
                    end
                    
                end
            end
        end
        
    end

end