classdef MyFilename
    % MyFilename is a class that holds filenames for easy manipulation
    properties
        Path = '.';
        Name = '';
        Ext = '';
    end
    
    methods
        function obj = MyFilename(fname)
            % MyFilename Constructor
            %   MyFilename(filename) starts the class with the provided
            %   filename
            if isa(fname, 'MyFilename')
                obj = fname;
            else
                [obj.Path, obj.Name, obj.Ext] = fileparts(fname);
            end
        end
        
        function fname = Full(obj)
            % Return full filename
            if isempty(obj.Path)
                obj.Path = '.';
            end
            fname = [obj.Path filesep obj.Name obj.Ext];
        end
        
        function obj = AddPrefix(obj, prefix)
            % adds a prefix to the name of the file (between path and name)
            obj.Name = [prefix, obj.Name];
        end
        
        function obj = SetPath(obj, path)
            % Change the path of the file
            obj.Path = path;
        end

        function obj = SetName(obj, name)
            % Change the name of the file
            obj.Name = name;
        end

        function obj = SetExt(obj, ext)
            % Change the extension of the file
            obj.Ext = ext;
        end
        
        function d = CreationDate(obj)
            % Get the creation date of the file
            data = dir(obj.Full);
            d = data(1).datenum;
        end
        
    end
end