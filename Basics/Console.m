classdef Console
    %% Auxulary tools for displaying data in the console (text data)
    
    methods (Static) 
        
        function t = Timer(curr)
            % Global timer for timing events
            persistent start;
            if nargin > 0
                start = curr;
            end
            t = toc(start);
        end
            
        function count = Counter(count)
            % Global counter of events
            persistent nevents;
            if nargin > 0
                if isnan(count)
                    nevents = 0;
                else
                    nevents = nevents + count;
                end
            end
            count = nevents;
        end
        
        
        function Color(varargin)
            % Create output with specific color (see cprintf for args)
            cprintf(varargin{:});
        end
        
        function Orange(varargin)
            % Display text in oragne
            Console.StrFix('[\b', ']\b', varargin{:});
        end
        
        function Red(varargin)
            % Display text in red
            fprintf(2, varargin{:});
        end
        
        function Bold(varargin)
            % Display text in bold
            Console.StrFix('<strong>', '</strong>', varargin{:});
        end
        
        function Green(varargin)
            % Display text in green
            cprintf('green', varargin{:});
        end
        
        function Link(varargin)
            % Display an hypertext
            Console.StrFix('<a href="">', '</a>', varargin{:});
        end

        function MatlabLink(filename, linenumber, varargin)
            % Display a link to specific file and linenumber
            Console.StrFix('<a href="matlab: opentoline(%s, %d)">', '</a>', filename, linenumber, varargin{:});
        end
        
        function StrFix(pre, suf, varargin)
            % add prefix and sufix to string
            fprintf([pre varargin{1} suf], varargin{2:end});
        end
    end
    
    %% styles
    methods (Static)
        function xml = StructToXML(s)
            % Convert strcture to XML
            xml = '';
            xml = sprintf(['%s' '<struct>\n'], xml);
            f = fields(s);
            for i=1:length(s)
                xml = sprintf(['%s' ' <entry>\n'], xml);
                for j=1:length(f)
                    xml = sprintf(['%s' '  <%s>%s</%s>\n'], xml, f{j}, num2str(s(i).(f{j})), f{j});
                end
                xml = sprintf(['%s' ' </entry>\n'], xml);
            end
            xml = sprintf(['%s' '</struct>\n'], xml);
        end

        function str = Format(v)
            % Format can convert data of (almost) any type to string 
            if ischar(v)
                v = regexprep(v, '\n', '\\n');
                v = regexprep(v, '\r', '\\r');
                v = regexprep(v, '\b', '\\b');
                str = ['''' v ''''];
                return;
            end
            if isscalar(v)
                str = sprintf('%g', v);
            elseif iscell(v)
                str = '';
                for i=1:length(v)
                    if i>1
                        str = sprintf('%s, %s', str, Console.Format(v{i}));
                    else
                        str = Console.Format(v{i});
                    end
                end
                str = sprintf('{ %s }', str);
            elseif isvector(v)
                if isempty(v)
                    str = '[]';
                else
                    str = ['[' sprintf('%g, ', v(1:end-1)) sprintf('%g', v(end)) ']'];
                end
            elseif ismatrix(v)
                sz = size(v); 
                str = '[';
                for i=1:size(v, 1)
                    if i<size(v, 1)
                        suff = ', ';
                    else
                        suff = '';
                    end
                    if ndims(v) > 2 %#ok<ISMAT>
                        str = [str, Console.Format(reshape(v(i, :), sz(2:end))), suff]; %#ok<AGROW>
                    else
                        str = [str, Console.Format(v(i, :)), suff]; %#ok<AGROW>
                    end
                end
                str = [str, ']'];
            end
        end    

        function str = MatlabFormat(v)
            % MatlabFormat can convert data of (almost) any type to string
            if ischar(v)
                v = regexprep(v, '\n', '\\n');
                v = regexprep(v, '\r', '\\r');
                v = regexprep(v, '\b', '\\b');
                str = ['''' v ''''];
                return;
            end
            if iscell(v)
                str = '';
                for i=1:length(v)
                    if i>1
                        str = sprintf('%s, %s', str, Console.MatlabFormat(v{i}));
                    else
                        str = Console.MatlabFormat(v{i});
                    end
                end
                str = sprintf('{ %s }', str);
            elseif isscalar(v)
                str = sprintf('%g', v);
            elseif isvector(v)
                if isempty(v)
                    str = '[]';
                else
                    str = ['[' sprintf('%g, ', v(1:end-1)) sprintf('%g', v(end)) ']'];
                end
            elseif ismatrix(v)
                sz = size(v); 
                str = '[';
                for i=1:size(v, 1)
                    if i<size(v, 1)
                        suff = '; ';
                    else
                        suff = '';
                    end
                    if ndims(v) > 2 %#ok<ISMAT>
                        str = [str, Console.MatlabFormat(reshape(v(i, :), sz(2:end))), suff]; %#ok<AGROW>
                    else
                        str = [str, Console.MatlabFormat(v(i, :)), suff]; %#ok<AGROW>
                    end
                end
                str = [str, ']'];
            end
        end    
        
        function str = StrudelError(err)
            % Display errors in a different format that is easier to parse
            % (Strudel format)
            Console.Red('#! %s\n', err.message);
            str = err.message;
            for i=1:length(err.stack)
                Console.Red('#! - <a href="matlab: opentoline(%s, %d)">%s (line %d)</a>\n', err.stack(i).file, err.stack(i).line, err.stack(i).name, err.stack(i).line);
                str = [str, sprintf('<a href="matlab: opentoline(%s, %d)">%s (line %d)</a>', err.stack(i).file, err.stack(i).line, err.stack(i).name, err.stack(i).line)];
            end
            str = Console.Strudel(struct('Error', str));
        end
        
        function str = Strudel(data, prefix)
            % Create a strudel format output
            if nargin < 2
                prefix = '/';
            end
            if isstruct(data) || isobject(data)
                f = fields(data);
                str = '';
                if length(data) > 1
                    for i=1:length(data)
                        str = sprintf('%s%s', str, Console.Strudel(data(i), [prefix f{i} sprintf('(%d)', i)]));
                    end
                elseif isempty(data)
                else
                    for i=1:length(f)
                        if isstruct(data.(f{i}))
                            str = sprintf('%s%s', str, Console.Strudel(data.(f{i}), [prefix f{i} '/']));
                        else
                            curr = data.(f{i});
                            %%
                            sz = size(curr);
                            szstr = [num2str(sz(1)) sprintf(',%d', sz(2:end))];
                            str = sprintf('%s@%s[%s(%s)] =', str, [prefix f{i}], class(curr), szstr);
                            if iscell(curr)
                                str = sprintf('%s%s', str, Console.Strudel(curr, [prefix f{i}]));
                            else
                                str = sprintf('%s %s', str, Console.Format(curr));
                            end
                            str = sprintf('%s\n', str);
                        end
                    end
                end
            elseif iscell(data)
                str = '';
                coord = cell(1, ndims(data));
                for i=1:numel(data)
                    [coord{:}] = ind2sub(size(data), i);
                    curr = data{i};
                    sz = size(curr);
                    szstr = [num2str(sz(1)) sprintf(',%d', sz(2:end))];
                    for j=1:length(coord)
                        if j>1
                            coordstr = sprintf('%s, %d', coordstr, coord{j});
                        else
                            coordstr = sprintf('%d', coord{j});
                        end
                    end
                    str = sprintf('%s@%s/{%s}[%s(%s)] =', str, prefix, coordstr, class(curr), szstr);
                end
            else
                
            end
        end

        function Warning(varargin)
            % Warning    Show warning message.
            %   Warning(me) displays the warning message in the
            %   'MException' object me (including call stack)
            %
            %   Warning(me, format, ...) displays the warning message with the
            %   header as defined in format (same as sprintf syntax)
            %   
            %   Warning(format, ...) displays the warning message as
            %   defined in format (sprintf syntax)
            %%
            if isa(varargin{1}, 'MException')
                %%
                err = varargin{1};
                if nargin > 1
                    str = sprintf(varargin{2:end});
                    Console.Red('#! %s: %s\n', str, err.message);
                else
                    Console.Red('#! %s\n', err.message);
                end
                for i=1:length(err.stack)
                    Console.Red('#! - <a href="matlab: opentoline(''%s'', %d)">%s (line %d)</a>\n', err.stack(i).file, err.stack(i).line, err.stack(i).name, err.stack(i).line);
                end
            else
                Console.Red('#! %s\n', sprintf(varargin{:}));
            end
        end
        
        function Stamp()
            % Time stamp with indiciation of the current file and line of
            % code. Helpful for debugging
            s = dbstack;
            if length(s) >= 2
                fprintf(2, 'STAMP: %s (%s, line %d)\n', s(2).file, s(2).name, s(2).line);
                drawnow('update')
            end
        end
        
        function Message(varargin)
            % Display a message
            styles = {@(x) Console.Bold(x), @(x) fprintf(x)};
            
            level = max(length(dbstack)-1, 1);
            if isnumeric(varargin{1})
                level = level + varargin{1};
                varargin = varargin(2:end);
            end
            s = sprintf(varargin{:});
            fprintf('# ');
            if level > 1
                fprintf([repmat(' ', 1, (level-2)*2) '- ']);
            end
            styles{min(level, length(styles))}(s);
            Console.NewLine;
        end

        function Write(varargin)
            % Display a message            
            styles = {@(x) Console.Bold(x), @(x) fprintf(x)};
            
            level = max(length(dbstack)-1, 1);
            if isnumeric(varargin{1})
                level = level + varargin{1};
                varargin = varargin(2:end);
            end
            s = sprintf(varargin{:});
            if level > 0
                fprintf('# ');
            end
            if level > 1
                fprintf([repmat(' ', 1, (level-2)*2) '- ']);
            end
            styles{min(level, length(styles))}(s);
        end
        
        function WriteLine(varargin)
            % Display a message ending with a newline
            styles = {@(x) Console.Bold(x), @(x) fprintf(x)};
            
            level = max(length(dbstack)-1, 1);
            if isnumeric(varargin{1})
                level = level + varargin{1};
                varargin = varargin(2:end);
            end
            s = sprintf(varargin{:});
            if level > 0
                fprintf('# ');
            else
                level = max(length(dbstack)-1, 1);
            end
            if level > 1
                fprintf([repmat(' ', 1, (level-2)*2) '- ']);
            end
            styles{min(level, length(styles))}(s);
            Console.NewLine;
        end
        
        function Header(varargin)
            % Display a header title            
            s = sprintf(varargin{:});
            fprintf(['#  . ' s '\n']);
        end
        
        function Title(varargin)
            % Display a title            
            s = sprintf(varargin{:});
            fprintf(['# ' s '\n']);
        end

        function Subtitle(varargin)
            % Display a subtitle            
            s = sprintf(varargin{:});
            fprintf(['# - ' s '\n']);
        end

        function SubSubtitle(varargin)
            % Display a sub-subtitle            
            s = sprintf(varargin{:});
            fprintf(['# - ' s '\n']);
        end
        
    end
    
    %% reports
    methods (Static)
        function Start(varargin)
            % Display message when starting a long process
            s = sprintf(varargin{:});
            fprintf(['#  . ' s ' ... ']);
        end
        
        function Done
            % Process done
            fprintf('[done]\n');
        end

        function Failed
            % Process failed
            fprintf('[fail]\n');
        end

        function Progress(idx, count)
            % Show progress report (step 'idx' of 'count')
            Console.Reprintf('# - %d/%d (%.1f%%)', idx, count, idx/count*100);
            if idx == count
                Console.NewLine;
            end
        end
        
        function rep = ProgRep(title, idx, count)
            % Show progress report with title (step 'idx' of 'count')
            persistent info;
            if idx <= 1 || isempty(info.startTime) || count ~= info.count
                Console.Reprintf(0, '');
                info.startTime = tic;
                info.nDigits = num2str(numel(num2str(fix(abs(count)))));
                info.count = count;
                info.Prev = toc(info.startTime);
            end
            curr = toc(info.startTime);
            if idx ~= count && curr - info.Prev < 0.5
                rep = false;
                return;
            end
            rep = true;
            if ~isempty(title)
                Console.Reprintf(['# - %s: %' info.nDigits 'd/%d (%.1f%%) [%.1f/%.1f] [x%.1fRT]'], title, idx, count, idx/count*100, curr, curr * count / idx, idx/curr);
            else
                Console.Reprintf(['# - %' info.nDigits 'd/%d (%.1f%%) [%s/%s] [x%.1fRT]'], idx, count, idx/count*100, sec2time(curr), sec2time(curr * count / idx), idx/curr);
            end
            info.Prev = curr;
        end
        
        function Report(title, idx, count)
            % Show progress report with title (step 'idx' of 'count')
            persistent startTime;
            if idx <= 1 || isempty(startTime)
                Console.Reprintf(0, '');
                startTime = tic;
            end
            Console.Reprintf('# - %s: %d/%d (%.1f%%) [%.1f/%.1f]', title, idx, count, idx/count*100, toc(startTime), toc(startTime) * count / idx);
        end
    end
    
    %% aux
    methods (Static)
        function Deprintf(varargin)
            % Delete characters from console (supply string to determine
            % number of characters to delete)
            s = sprintf(varargin{:});
            for i=1:length(s)
                s = sprintf('%s%s', s, '\b');
            end
            fprintf(s);
        end
        
        function n = Reprintf(varargin)
            % Display text while overwriting previous text
            persistent curr;
            if isempty(curr)
                curr = 0;
            end
            if isnumeric(varargin{1})
                prev = varargin{1};
                in = varargin(2:end);
            else
                if curr > 0
                    prev = curr;
                else
                    prev = 0;
                end
                in = varargin;
            end
            for i=1:prev
                fprintf('\b');
            end
            s = sprintf(in{:});
            n = length(s);
            curr = n;
            s = regexprep(s, '%', '%%');
            fprintf(s);
        end
        
        function NewLine
            % Newline
            Console.Reprintf(0, '');
            fprintf('\n');
        end
        
        function pos = MarkPosition()
            % Experimental. Don't use
            mde = com.mathworks.mde.desk.MLDesktop.getInstance;
            cw = mde.getClient('Command Window');
            xCmdWndView = cw.getComponent(0).getViewport.getComponent(0);
            pos = xCmdWndView.getCaretPosition;
        end
        
        function GoTo(pos)
            % Experimental. Don't use
            curr = Console.MarkPosition;
            for i=1:curr-pos
                fprintf('\b');
            end
        end
    end
end