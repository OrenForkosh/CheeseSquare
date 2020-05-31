classdef Q
    % A Random collection of useful (and sem-useful) functions
    methods (Static = true)
        function s = signif2sym(p)
            % convert p-value to string representing its significance (like
            % '*' or 'ns)
            if p > 0.1
                s = 'ns';
            elseif p > 0.05
                s = '#';
            elseif p > 0.01
                s = '*';
            elseif p > 0.001
                s = '**';
            else
                s = '***';
            end
                
                
        end
        
        function b = isserver()
            % check if matlab is in server (nodesktop) mode
            b = ~usejava('Desktop');
        end
        
        function r = matfun(m, fun)
            % run function on each row of matrix
            r = [];
            for i=1:size(m,1)
                r(i, :) = fun(m(i, :));
            end
        end
        
        function res = accumrows(subs, val, fun, fillval)
            % ACCUMROWS create matrix by accumulating rows of val specified
            % by subs. Similar to matlab's accumarray but works on
            % matrices. Can also specifiy default 'fillval' for empty
            % entries
            if nargin < 3
                fun = {};
            else
                fun = { fun };
            end
            if nargin < 4
                fillval = 0;
            end
            for i=1:size(val, 2)
                if i==1
                    a = accumarray(subs, val(:, i), [], fun{:}, fillval);
                    res = nan(length(a), size(val, 2));
                    res(:, 1) = a;
                else
                    res(:, i) = accumarray(subs, val(:, i), [], fun{:}, fillval);
                end
            end
        end
        
        function c = vprintf(fmt, x)
            % printf for each entry in vector x
            %% supports generating cells if given vectors
            c = cell(1, length(x));
            for i=1:length(x)
                c{i} = sprintf(fmt, x(i));
            end
        end
        
        function TF = isnempty(varargin)
            % is-not-empty
            TF = ~isempty(varargin(:));
        end
        
        function r = iftheneval(cond, res1, res2)
            % If 'cond' is satisfied than evaluate 'res1' else evaluate 'res2'
            if cond
                if char(res1)
                    r = evalin('caller', res1);
                else
                    r = res1;
                end
            else
                if char(res2)
                    r = evalin('caller', res2);
                else
                    r = res2;
                end
            end
        end
        
        function r = ifthen(cond, res1, res2)
            % If 'cond' is satisfied than return 'res1' else return 'res2'
            r = res1;
            r(~cond) = res2(~cond);
        end
        
        function x = torow(x)
            % convert x to a row
            x = x(:)';
        end

        function x = tocol(x)
            % convert x to a column
            x = x(:);
        end
        
        function v = minval(t)
            % minimal possible value of a specified type 't'
            if isinteger(t)
                v = intmin(class(t));
            elseif isfloat(t)
                v = -inf;
            elseif islogical(t)
                v = false;
            else
                error;
            end
        end
        
        function m = maxabs(x, varargin)
            % the maximum of abs(t)
            x = Q.pushdim(x, varargin{:});
            [~, id] = max(abs(x)); 
            m = x(sub2ind(size(x), id, 1:length(id)));
            m = Q.popdim(m, varargin{:});
        end
        
        function f = toframes(y, sz, overlap, padding)
            % divide y into frames of size sz and overlap
            if nargin < 4
                padding = [];
            end
            if ~isempty(padding)
                y = [padding * ones(overlap, 1); y(:); padding * ones(overlap, 1)];
            else
                y = y(:);
            end
            
            ncol = fix((length(y)-overlap)/(sz-overlap));
            colindex = 1 + (0:(ncol-1))*(sz-overlap);
            rowindex = (1:sz)';
            if islogical(y)
                f = false(sz, ncol); 
            else
                f = zeros(sz, ncol, class(y)); %#ok<*ZEROLIKE>
            end
            f(:) = y(rowindex(:,ones(1,ncol))+colindex(ones(sz,1),:)-1);
        end
        
        function m = quickmedian(y)
            % quick algorithm for finding the median
            hi = max(y);
            lo = min(y);
            m = (lo + hi) / 2;
            h = length(y)/2;
            while true
                count = sum(y > m);
                if abs(count - h) < 1
                    break;
                elseif count < h
                    hi = m;
                else
                    lo = m;
                end
                m = (lo + hi) / 2;
            end
            m = y(Q.argmin(abs(y-m)));
        end
        
        function [data, index] = sort(data, greaterthen, dim)
            % uses Bubble sort because of lazyness
            if nargin < 2
                greaterthen = @(x,y) x>y;
            end
            if nargin < 3
                dim = Q.argmax(size(data));
            end
            order = [dim Q.exclude(ndims(data), dim)];
            data = permute(data, order);
            
            index = 1:size(data, 1);
            found = true;
            if iscell(data)
                cmp = @(x, i, j) greaterthen(x{i, :}, x{j, :});
            else
                cmp = @(x, i, j) greaterthen(x(i,:), x(j,:));
            end
            while found
                found = false;
                for i=2:length(index)
                    if cmp(data, i-1, i)
                        found = true;
                        temp = data(i-1, :);
                        data(i-1, :) = data(i, :);
                        data(i, :) = temp;
                        temp = index(i-1);
                        index(i-1) = index(i);
                        index(i) = temp;
                    end
                end
            end
            data = ipermute(data, order);
        end
        
        function order = hiersort(D)
            % runs hierarchical sorting
            if ~isvector(D)
                D = squareform(D,'tovector');
            end
            z = linkage(D);
            order = optimalleaforder(z, D);
            if nargout == 0
                d = squareform(D);
                d = d(order, order);
                Plot.Hinton(d);
            end
        end
        
        function c = combinations(varargin)
            % all possible combinations of entries from several arrays. Use
            % Matlab's combvec instead
            if nargin == 1 && iscell(varargin{1})
                c = Q.combinations(varargin{1}{:});
                return;
            end
            if nargin == 1
                c = varargin{1}(:);
            else
                other = Q.combinations(varargin{2:end});
                nother = size(other, 1);
                c = zeros(length(varargin{1}) * nother, size(other, 2)+1);
                for i=1:length(varargin{1})
                    c((i-1)*nother+1:i*nother, :) = [ones(nother, 1) * varargin{1}(i), other];
                end
            end
        end
        
        function out = struct2xml(param, filename)
            % serialize struct to xml format
            if ischar(param)
                name = param;
                param = evalin('caller', name);
            else
                name = 'root';
            end
            doc  = com.mathworks.xml.XMLUtils.createDocument(name);
            root = doc.getDocumentElement;
            out = xmlwrite(parse(param, doc, root));
            if nargin > 1
                dlmwrite(filename, out);
            end
            function obj = parse(param, doc, root)
                if isstruct(param)
                    fields = fieldnames(param);
                    for i  = 1 : numel(fields)
                        if isstruct(param.(fields{i}))
                            write(doc, root, fields{i}, '')
                            parse(param.(fields{i}), doc, root.getLastChild)
                        else
                            write(doc, root, fields{i}, param.(fields{i}))
                        end
                    end
                else
                    doc = 'could not convert';
                end
                obj = doc;
            end
            function write(doc, root, pname, param)
                if isempty(param) == 1
                    root.appendChild(doc.createElement(pname));
                else
                    elem = doc.createElement(pname);
                    if ismatrix(param)
                        val = Console.Format(param);
                    else
                        val = Console.Format(param);
                    end
                    elem.appendChild(doc.createTextNode(val));
                    root.appendChild(elem);
                end
            end
        end
        
        function str = capital(str)
            % capitalize a string
            if isempty(str)
                return;
            end
            str = [upper(str(1)) str(2:end)];
        end

        function str = capitalize(str)
            % capitalize a string
            if isempty(str)
                return;
            end
            idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
            str(idx)=upper(str(idx));
        end
        
        function args = defaultargs(varargin)
            % process function arguments. Use Matlab's inputParser instead
            stack = dbstack;
            %%
            offset = 1;
            if islogical(varargin{1})
                output = varargin{1};
                offset = 2;
            else
                output = true;
            end

            %%
            if iscell(varargin{offset})
                if length(varargin{offset}) == 1 && isstruct(varargin{offset}{1})
                    args = varargin{offset}{1};
                else
                    args = struct();
                    for i=1:2:length(varargin{offset})
                        args.(varargin{offset}{i}) = varargin{offset}{i+1};
                    end
                end
            else
                args = varargin{offset};
            end
            if isfield(args, 'verbose')
                output = args.verbose;
            end
            if isfield(args, 'Verbose')
                output = args.Verbose;
            end
            %%
            if length(stack) > 1 && output
                Console.Message('setting arguments for <a href="matlab: opentoline(%s, %d)">%s</a>', stack(2).file, stack(2).line, stack(2).name);
            end
            %%
            namelen = 20;
            vallen = 12;
            for i=offset+1:2:length(varargin)
                if isfield(args, varargin{i})
                    if output
                        val = Console.Format(args.(varargin{i}));
                        if length(val) > 50
                            val = [val(1:50), ' ... '];
                        end
                        Console.Message(1, ['%-' num2str(namelen) 's = %-' num2str(vallen) 's'], ['''' varargin{i} ''''], val);
                    end
                else
                    args.(varargin{i}) = varargin{i+1};
                    if output
                        name = ['''' varargin{i} ''''];
                        if length(name) > namelen
                            currvallen = vallen - (length(name) - namelen);
                        else
                            currvallen = vallen;
                        end
                        val = Console.Format(varargin{i+1});
                        if length(val) > 50
                            val = [val(1:50), ' ... '];
                        end
                        Console.Message(1, ['%-' num2str(namelen) 's = %-' num2str(currvallen) 's [default]'], name, val);
                    end
                end
            end
            
        end
        
        function setifneeded(varname, value)
            % Set value of variable
            %evalin('caller', sprintf('if ~exist(''%s'', ''var''); %s = %s; end', varname, varname, Console.MatlabFormat(value)));
            evalin('caller', sprintf('%s = %s', varname, Console.MatlabFormat(value)));
        end

        function a = getindex(arr, varargin)
            % Get specific index of a matrix. Useful for inline functions
            a = arr(varargin{:});
        end
        
        function a = getargout(index, func, varargin)
            % 
            X = cell(1, index);
            if nargin > 2
                [X{:}] = func(varargin{:});
            else
                [X{:}] = func();
            end
            a = X{index};
        end
        
        function c = shift(c, k, pad)
            % c = shift(c, k, pad=0)
            %   like circshift only uses pad as fill value
            if nargin < 3
                pad = 0;
            end
            c = circshift(c, k);
            for i=1:length(k)
                if k(i) == 0
                    continue;
                end
                c = Q.pushdim(c, i);
                if k(i) > 0
                    c(1:k(i), :) = pad;
                else
                    c(end-k(i):end, :) = pad;
                end
                c = Q.popdim(c, i);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % struct operations
        function val = getfield(var, name)
            % get specific field of struct
            fields = textscan(name, '%s', 'Delimiter', '.');
            val = getfield(var, fields{1}{:});
        end
        
        function var = setfield(var, varargin)
            % set field of struct
            for i=1:2:length(varargin)
                fields = textscan(varargin{i}, '%s', 'Delimiter', '.');
                var = setfield(var, fields{1}{:}, varargin{i+1});
            end
        end

        function to = cpfield(from, to, names, ignore)
            % copy fields between structs
            if nargin == 2
                f = fieldnames(from);
                for i=1:length(f)
                    to = Q.cpfield(from, to, f{i});
                end
                return;
            end
            if nargin == 3
                ignore = false;
            end
            if ischar(names) 
                if Q.isfield(from, names)
                    to = Q.setfield(to, names, Q.getfield(from, names));
                elseif ~ignore
                    throw(MException('Q:nonExistentField', 'Reference to non-existent field ''%s''', names));
                end
            else
                for i=1:length(names)
                    to = Q.cpfield(from, to, names{i}, ignore);
                end
            end
        end
        
        function tf = isfield(var, name)
            %% like matlab's isfield but also approves properties and
            %% can support several levels, like: Q.isfield(s, 'a.b.c.d');
            %%
            c = strsplit(name, '.');
            curr = var;
            tf = true;
            for i=1:length(c)
                if ~isfield(curr, c{i}) && ~isprop(curr, c{i})
                    tf = false;
                    break;
                end
                curr = curr.(c{i});
            end
        end
        
        function c = struct2cell(var)
            % convert a struct to a cell including fieldnames
            a = struct2cell(var); 
            a(:, 2) = a(:, 1); 
            a(:, 1) = fieldnames(var); 
            a = a';
            c = a(:);
        end
    end
    methods (Static = true)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % mathematics
        function [pc, vals] = pca(x)
            % Principal component analysis (simply uses Matlab's)
            [coeff,score,latent,tsquared,explained] = pca(x);
            if nargout == 0
                %%
                subplot(2,2,1);
                plot(cumsum(explained), 'o');
                grid minor
                grid on
                Fig.Fix
            end
        end
        
        function y = gcd(x)
            % greatest common divider for a vector
            J = size(x, 2);
            while J > 1
                if mod(J, 2) ~= 0;
                    x(:, 1:(J+1)/2) = [reshape(gcd(Q.tocol(x(:, 1:2:J-1)), Q.tocol(x(:, 2:2:J))), size(x,1), floor(J/2)), x(:, J)];
                    J = (J+1) / 2;
                else
                    x(:, 1:J/2) = reshape(gcd(Q.tocol(x(:, 1:2:J)), Q.tocol(x(:, 2:2:J))), size(x,1), floor(J/2));
                    %x(1:J/2) = gcd(x(:, 1:2:J), x(:, 2:2:J));
                    J = J / 2;
                end
            end
            y = x(:, 1);
            return;
            %%
            J = length(x);
            while J > 1
                j=1;
                for i=1:2:J
                    if i+1 <= J
                        x(j) = gcd(x(i), x(i+1));
                    else
                        x(j) = x(i);
                    end
                    j = j + 1;
                end
                J = j - 1;
            end
            y = x(1);
            %%
            
        end
        
        
        function [d, Y] = derivative(y, x, n)
            % gradient of y(x) smoothed at n points
            if nargin < 3
                n = 1;
            end
            if nargin < 2
                x = 1:length(y);
            end

            if length(y) == 1
                d = 0;
                Y = y;
                return
            elseif length(y) == 2
                d = ones(size(y)) * (y(2)-y(1)) / (x(2) - x(1));
                Y = y;
                return
            end
            
            
            %%
            f = fit(x, y, 'smoothingspline', 'SmoothingParam', 1-1/n);
            if nargout > 1
                Y = f(x);
            end
            %f = fit(x, y, 'poly2');
            d = gradient(f(x), x);

            
            %%
            return;
            
            %%
%             Y = im2col(Q.torow(y), [1, n+1], 'sliding');
%             X = im2col(Q.torow(x), [1, n+1], 'sliding');
            if length(y) == 1
                d = 0;
                return
            elseif length(y) == 2
                d = [1 1] * (y(2)-y(1));
                return
            end
            N = min(n + 1, length(y));
            d = zeros(size(y));
            for i=1:length(y)
                pre = round((i - 1) / length(y) * N);
                range = i-pre:i+(N-pre-1);
                p = polyfit(x(range), y(range), 1);
                d(i) = p(1);
            end
        end
        
        function [c, p] = nancorr(x, y)
            % correlation which ignores NANs
            if nargin == 1
                y = x;
            end
            c = zeros(size(x, 2), size(y, 2));
            p = zeros(size(x, 2), size(y, 2));
            for i=1:size(x, 2)
                X = x(:, i);
                for j=1:size(y, 2)
                    Y = y(:, j);
                    valid = ~isnan(X) & ~isnan(Y);
                    [c(i,j) p(i,j)] = corr(X(valid), Y(valid));
                end
            end
        end
        
        function x = pushdim(x, dim)
            % pushdim(x, dim) makes dim the first dimension
            if nargin == 1 || dim == 1
                return;
            end
            x = permute(x, [dim Q.exclude(ndims(x), dim)]);
        end

        function x = popdim(x, dim)
            % popdim(x, dim) returns first dimension to dimension dim
            if nargin == 1 || dim == 1
                return;
            end
            seq = cat(2, 2:dim, 1, dim+1:ndims(x));
            x = permute(x, seq);
        end
        
        function s = smooth(x, width)
            % smooth vector 
            s = conv2(ones(width,1)/width, 1, x, 'same');
            return
            if dim == 1
                s = cumsum(x - [zeros(len,size(x, 2)); x(1:end-len, :)]);
            elseif dim == 2
                s = cumsum(x - [zeros(len,size(x, 2)); x(1:end-len, :)]);
            else
                error
            end
               
        end
        
        function s = nansmooth(x, len)
            % smoothes s using a runing average smoothing. Ignores nans.
            x1 = x; 
            x1(isnan(x)) = 0;
            s = convn(x1(:), ones(len, 1), 'same');
            m = convn(~isnan(x(:)), ones(len, 1), 'same');
            s = reshape(s ./ m, size(x));
        end

        function s = nanblock(x, len)
            % smoothes s using a block smoothing. Ignores nans.
            n = ceil(length(x) / len);
            x(end+1:n*len) = nan;
            x = reshape(x, len, n);
            s = nanmean(x);
        end
        
        function s = discspace(from, to, count)
            % like linspace only discrete values
            s = round(linspace(from, to, count));
        end
        
        function i = index(x)
            % return indices for each entries in a vector 'x'
            i = 1:length(x(:));
            i = reshape(i, size(x));
        end
        
        function j = exclude(n, i)
            % exclude(n, i) returns the sequence from 1 to n excluding i,
            % i.e. 1 ,.., i-1, i+1 ,..., n
            j = 1:n;
            if nargin > 1
                j(i) = [];
            end
        end
        
        function v = inrange(x, m1, m2)
            % returns true if x is between m1 and m2
            if nargin == 2 && length(m1) == 2
                m2 = m1(2);
                m1 = m1(1);
            end
            v = x >= m1 & x <= m2;
        end

        function map = overlap(x1, x2, m1, m2)
            % overlap(x1, x2, m1, m2), returns map of segments
            %    [m1(i), m2(i)] that overlap with the segments [x1, x2]
            if length(x1) == 1
                map = (x1 >= m1 & x1 <= m2) | (x2 >= m1 & x2 <= m2) | (x1 <= m1 & x2 >= m2);
            else
                gt11 = bsxfun(@ge, x1(:), m1(:)');
                lt12 = bsxfun(@le, x1(:), m2(:)');
                gt21 = bsxfun(@ge, x2(:), m1(:)');
                lt22 = bsxfun(@le, x2(:), m2(:)');
                lt11 = bsxfun(@le, x1(:), m1(:)');
                gt22 = bsxfun(@ge, x2(:), m2(:)');
                map = (gt11 & lt12) | (gt21 & lt22) |  (lt11 & gt22);
            end
        end
        
        
        function x = padtop(x, padsize, padval)
            % pad top part of matrix
            if nargin < 2; padsize = 1; end
            if nargin < 3; padval = 0; end
            x = padarray(x, padsize, padval, 'pre');
        end

        function x = padbottom(x, padsize, padval)
            % pad bottom part of matrix
            if nargin < 2; padsize = 1; end
            if nargin < 3; padval = 0; end
            x = padarray(x, padsize, padval, 'post');
        end

        function x = padleft(x, padsize, padval)
            % pad left part of matrix
            if nargin < 2; padsize = 1; end
            if nargin < 3; padval = 0; end
            x = padarray(x, [0 padsize], padval, 'pre');
        end

        function x = padright(x, padsize, padval)
            % pad right part of matrix
            if nargin < 2; padsize = 1; end
            if nargin < 3; padval = 0; end
            x = padarray(x, [0 padsize], padval, 'post');
        end
        
        function idx = FindFirst(x, dim)
            % find first true entry in x
            if nargin < 2
                dim = 1;
            end
            c = cumsum(x~=0, dim);
            idx = min(sum(c == 0) + 1, size(x, dim));
        end

        function idx = FindLast(x, dim)
            % find last true entry in x
            if nargin < 2
                dim = 1;
            end
            c = cumsum(flipdim(x, dim)~=0, dim);
            idx = max(size(x, dim) - sum(c == 0), 1);
        end
        
        function y = frac(x)
            % fractional part of a number
            y = abs(x - fix(x));
        end
        
        function x = normalize(x)
            % n = Normalize(x) ensures all values of x are between 0 and 1
            m1 = min(x(:));
            m2 = max(x(:));
            x = (x - m1 - eps) / (m2 - m1);
        end
        
        function [r, N] = rank(x, dim, nranks)
            % for each sample returns it's rank
            % for example, for x = [1 4 2 9 5 6], the function will return 
            % rank(x) = [1 3 2 6 4 5]
            if nargin > 1
                x = permute(x, [dim, Q.exclude(ndims(x), dim)]);
            else
                dim = [];
                if size(x, 1) == 1
                    dim = 2;
                    x = permute(x, [dim, Q.exclude(ndims(x), dim)]);
                end
            end
            [sortx, o] = sort(x, 1);
            N = size(x, 1);
            %%
            seq = 1:N;
            r = zeros(size(x));
            idx = cell(1, ndims(r));
            [idx{:}] = ind2sub(size(r), 1:numel(r));
            idx{1} = o(:)';
            r(sub2ind(size(r), idx{:})) = repmat(seq, 1, numel(o)/N);
            %% map identical values to same rank
            seq = (1:size(x, 1))';
            for i=1:size(x, 2)
                c = cumsum(Q.padtop(diff(sortx(:, i)) ~= 0, 1, true));
                valid = Q.padtop(diff(sortx(:, i)) ~= 0, 1, true);
                s = round(Q.accumrows(cumsum(valid), seq, @mean));
                s = s(c);
                r(o(:, i), i) = s;
            end
            %%
            %%
            if ~isempty(dim)
                b = 1:dim-1;
                a = dim+1:ndims(x);
                r = permute(r, [b+1 1 a]);
            end
            %%
            if nargin >= 3
                r = floor((r - 1) / N * nranks) + 1;
            end
        end
        
        
        function y = pcaR(x)
            % rubost pca
            y = exact_alm_rpca(x);
        end
        
        function y = logit(x)
            % logit function
            y = log(x ./ (1-x));
        end
        
        function [nx, x1, k] = normtrans(x)
            %%
            x = double(data(:, 6))';
            f = @(x0, x) x0(4) * (sinh(x0(2) * (x - x0(1))) - x0(3));
            %f = @(x0, x) (x0(2) * x + x0(1)).^3;
            %f = @(x0, x) normcdf(x, x0(1), x0(2))
            y = Q.nwarp(x);
            x0 = nlinfit(x, y, f, [median(x), Q.stdR(x), mean(y), std(y)]);
            nx = f(x0, x);
            subplot(2,2,1);
            h = histogram([nx(:)]);
            subplot(2,2,2);
            h = histogram([x(:)]);
            subplot(2,2,3);
            h = histogram([y(:)]);
            
            %%
            x = rand(1, 10000);
            f = @(x0, x) boxcox(x, x0(1));
            y = Q.nwarp(x);
            x0 = nlinfit(x, y, f, .5);
            nx = f(x0, x);
            subplot(2,2,1);
            h = histogram([nx(:)]);
            subplot(2,2,2);
            h = histogram([x(:)]);
            subplot(2,2,3);
            h = histogram([y(:)]);
            
            
            %%
            x1 = fminsearch(@(x0) (x - y).^2, [mean(x), std(x)]);
            nx = f(x, x1);
            k = kurtosis(f(x, x1))
            subplot(1,2,1);
            h = histogram([nx(:)]);
            subplot(1,2,2);
            h = histogram([x(:)]);
            
        end
        
        function y = meannorm(x, varargin)
            % normalize to have zero mean
            y = bsxfun(@minus, x, mean(x, varargin{:}));
        end
        
        function y = nwarp(x, varargin)
            %y = Q.znormR(x, varargin{:});
            %return;
            % warp samples to match a normal distribution
            [rank, N] = Q.rank(x, varargin{:});
            %[~, rank, N] = Q.quant(x, varargin{:});
            %n = size(x, varargin{:});
            if false
                alpha = .5;
                minval = 1 - alpha^(1/N);
                maxval = 1 - minval;
                vals = norminv(linspace(minval, maxval, N));
            else
                vals = norminv(linspace(1/N, 1-1/N, N));
                %vals = linspace(0, 1, N);
            end
            y = vals(rank);
            y = reshape(y, size(x));
        end

        function y = nwarpold(x, varargin)
            % warp samples to match a normal distribution
            [rank, N] = Q.rank(x, varargin{:});
            vals = norminv(linspace(1/N, 1-1/N, N));
            y = vals(rank);
            y = reshape(y, size(x));
        end
        
        function y = nwarpR(x, dim, n)
            % warp samples to match a normal distribution in a robust way
            if nargin < 2
                dim = 1;
            end
            if nargin < 3
                n = ceil(size(x, dim)/10);
            end
            % warp samples to match a normal distribution
            [rank, N] = Q.rank(x, dim, n);
            vals = norminv(linspace(1/N, 1-1/N, n));
            y = vals(rank);
            y = reshape(y, size(x));
        end
        
        function x = discretize(x, nbins)
            % discretize values in a vector
            m1 = min(x(:));
            m2 = max(x(:));
            ext = max(abs(m1), abs(m2));
            sp = linspace(-ext, ext, nbins);
            mid = (sp(1:end-1)+sp(2:end))/2;
            sp(end) = inf;
            sp(1) = -inf;
            [~,~,ind] = histcounts(x(:), sp);
            
            x = reshape(mid(ind), size(x));

        end
        
        function y = znorm(x, varargin)
            % z-norm (set mean to zero and std to one)
            cx =  bsxfun(@minus, x, nanmean(x, varargin{:}));
            y = bsxfun(@rdivide, cx, nanstd(x, 0, varargin{:}));
        end
        
        function y = znormR(x, varargin)
            % Robust z-norm (set mean to zero and std to one)
            %%
            y = bsxfun(@rdivide, bsxfun(@minus, x, Q.meanR(x, varargin{:})), Q.stdR(x, varargin{:}));
        end
        
        function [m] = agglomerate(m, ithresh, minn)            
            [m, o] = sort(m);
            thresh = ithresh;
            neigs = -1;
            while neigs < minn
                id = cumsum([1; diff(m) ./ m(2:end) > thresh]);
                neigs = max(id);
                thresh = thresh / 2;
            end          
            newm = zeros(size(m));
            for i=1:neigs
                idx = find(id == i, 1, 'last');
                newm(id == i) = m(idx);
            end
            m(o) = newm;
            %%
            return
            %%
            thresh = ithresh;
            neigs = -1;
            while neigs < minn
                id = cumsum([1; diff(e) ./ e(2:end) > thresh]);
                neigs = max(id);
                thresh = thresh / 2;
            end
            newe = zeros(1, neigs);
            newW = zeros(size(m, 1), neigs);
            for i=1:neigs
                idx = find(id == i, 1, 'last');
                newe(i) = e(idx);
                newW(:, i) = W(:, idx);
                %n(:, i) =
            end
            e = diag(newe);
            %%
            return;
            [V, e] = eig(m);
            
            [~,maxind] = max(abs(V), [], 1);
            [d1, d2] = size(V);
            colsign = sign(V(maxind + (0:d1:(d2-1)*d1)));
            colsign = sign(mean(V));
            V = bsxfun(@times, V, colsign);
            
            
            return;
            %%
            for i=1:size(K, 2)
                M = m - e(1) * eye(size(m));
                K = M;
                K(:, i) = 0;
                det(K) / det(M)
            end
            
            %%
            %return;
            %%
            if true
                e = eig(m);
                n = 2;
                for i=1:length(e)
                    %%
                    K = m - e(i) * eye(size(m));
                    B = K(1:end-n, :);
                    %l = lasso(B, [zeros(size(K, 1), 1); ones(size(v, 2), 1)], 'NumLambda', 2);
                    %V(:, i) = l(:, 1);
                    l = null(B);
                    V(:, i) = l(:, Q.argmin(pdist2(null(K)', l')));
                    V(:, i) = V(:, i) / norm(V(:, i));
                end
            elseif false
            e = eig(m);
            n = 5;
            for i=1:length(e)
                %%
                K = m - e(i) * eye(size(m));
                
                v = randn(size(K, 2), n);
                v = bsxfun(@rdivide, v, sqrt(sum(v.^2)));
                B = [K; v'];
                %l = lasso(B, [zeros(size(K, 1), 1); ones(size(v, 2), 1)], 'NumLambda', 2);
                %V(:, i) = l(:, 1);
                V(:, i) = lsqlin(B, [zeros(size(K, 1), 1); ones(size(v, 2), 1)]);
                V(:, i) = V(:, i) / norm(V(:, i));
            end
            else
                [V, e] = eig(m);
            end
            [~,maxind] = max(abs(V), [], 1);
            [d1, d2] = size(V);
            colsign = sign(V(maxind + (0:d1:(d2-1)*d1)));
            V = bsxfun(@times, V, colsign);
            
            %%
            
            return
            %%
            [v, e] = eig(m);
            v = flipdim(v, 2);
            [~,maxind] = max(abs(v), [], 1);
            [d1, d2] = size(v);
            colsign = sign(v(maxind + (0:d1:(d2-1)*d1)));
            v = bsxfun(@times, v, colsign);
        end
        
        function m = matR(m, beta)
            [u,s,v] = svd(m);
            dm = diag(s);
            q = quantile(dm, beta);
            dm(dm > q) = q;
            m = u * diag(dm) * v';
        end
        
        function [im, m2] = invR(m, beta)
            % Robust matrix invers
            %%
            [u,s,v] = svd(m);
            dm = diag(s);
            q = quantile(dm, beta);
            dm(dm < q) = q;
            im = u * diag(1./dm) * v';
            if nargout > 1
                m2 = u * diag(dm) * v';
            end
        end
        
        function m = covR2(x, beta)
            % Robust covariance
            m = cov(x);
            m = m  + beta * eye(size(m));
        end
        
        function m = covR(x, beta)
            % Robust covariance
            m = zeros(size(x, 2));
            x = bsxfun(@minus, x, median(x));
            for i=1:size(x, 2)
                for j=1:size(x, 2)
                    %m(i, j) = sum(x(:, i) .* x(:, j)) / (size(x, 1) - 1);
                    m(i, j) = median(x(:, i) .* x(:, j));
                end
            end
        end
        
        function m = meanR(x, varargin)
            % Robust mean
            m = nanmedian(x, varargin{:});
        end
        
        function m = stdR(x, varargin)
            % computes the std using mad (median absolute deviation), and scales accordingly
            % so that it equals the std for normaly distributed data
            m = nanmedian(abs(bsxfun(@minus, x, nanmedian(x, varargin{:}))), varargin{:}) / norminv(3/4);
        end

        function m = fitR(x, dim, nstd)
            % Robust fit
            thresh = normcdf(0, nstd);
            q = quantile(x, [thresh, 1-thresh], dim);
            m = 2 * bsxfun(@rdivide, bsxfun(@minus, x, Q.meanR(x, dim)), diff(q, 1, dim)) * nstd;
        end
        
        function [bins, idx, count] = binify(data, nbins)
            % put data into bins
            s = sort(data);
            %%
            nvals = nbins - 2;
            bins = [s(1) zeros(1, nvals) inf];
            offset = 1;
            succ = true;
            for i=0:nvals
                idx = offset;
                bins(i+1) = s(idx);
                offset = round(offset + (length(s) - offset + 1) / (nvals - i + 1));
                offset = max(offset, find(s > bins(i+1), 1));
                if isempty(offset)
                    succ = false;
                    break;
                end
            end
            if ~succ
                us = unique(s);
                if length(us) >= nvals
                    bins(2:1+nvals) = us(round(linspace(1, length(us), nvals)));
                else
                    bins(2:1+nvals) = linspace(0, 2 * max(us), nvals);
                end
            end
            %%
            if nargout > 1
                [count, idx] = histc(data, bins);
            end
        end
        
        
        function m = stderr(x, dim)
            % computes standard error
            if nargin == 1
                dim = 1;
            end
            m = std(x, 0, dim) / sqrt(size(x, dim));
        end

        function m = nanstderr(x, dim)
            % computes standard error while ignoring nans
            if nargin == 1
                dim = 1;
            end
            m = nanstd(x, 0, dim) ./ sqrt(sum(~isnan(x), dim));
        end
        
        function p = peaks(mat, dim, minpeak)
            % find peaks of a matrix quickly
            if nargin < 3
                minpeak = -inf;
            end
            if nargin < 2
                dim = 1;
            end
            if dim ~= 1
                seq = 1:ndims(mat);
                seq(dim) = [];
                mat = permute(mat, [dim seq]);
            end
            d = diff(mat, 1);
            p = [false(1, size(mat, 2)); d(1:end-1, :) > 0 & d(2:end, :) < 0; false(1, size(mat, 2))];
            if minpeak ~= -inf
                p(mat < minpeak) = false;
            end
            if dim ~= 1
                p = permute(p, [2:dim 1 dim+1:ndims(mat)]);
            end
        end
        
        function m = softabs(x, alpha)
            % continous versions of abs
            m = sqrt(x.^2 + alpha^2);
        end
        
        function m = softmax(x, k, varargin)
            % m = Softmax(x, k, dim) smooth and differentiable approximation of maximum
            %    k is the scale
            if nargin < 2
                k = 1;
            end
            m = log(sum(exp(k * x), varargin{:}))/k;
            %m = sum(x .* expx, dim) ./ sum(expx, dim);
        end

        function mm = minmax(x, dim)
            % returns both the minimum and maximum of a matrix
            if nargin < 2
                mm = [min(x); max(x)];
            else
                mm = cat(dim, min(x, [], dim), max(x, [], dim));
            end
            
        end
        
        function idx = argmax(x, dim)
            % returns index of maximum
            if nargin <= 1
                [~, idx] = max(x);
            else
                [~, idx] = max(x, [], dim);
            end
        end

        function idx = argmaxi(x, dim)
            % like argmax, but returns the linear index of the cell in the matrix 
            if nargin <= 1
                [~, idx] = max(x);
            else
                [~, idx] = max(x, [], dim);
            end
            
            if nargin <= 1 || dim == 1
                idx = sub2ind(size(x), idx, 1:length(idx));
            else
                idx = sub2ind(size(x), (1:length(idx))', idx);
            end
        end
        
        function idx = argmin(x, dim)
            % returns index of minimum
            if nargin <= 1
                [~, idx] = min(x);
            else
                [~, idx] = min(x, [], dim);
            end
        end
        
        function [idx, coord] = argmaxnd(x)
            % returns index of maximum of a multidimensional array
            [~, idx] = max(x(:));
            a = cell(1, ndims(x));
            [a{1:ndims(x)}] = ind2sub(size(x), idx);
            coord = cell2mat(a);
        end
        
        function [idx, coord] = argminnd(x)
            % returns index of minimum of a multidimensional array
            [~, idx] = min(x(:));
            a = cell(1, ndims(x));
            [a{1:ndims(x)}] = ind2sub(size(x), idx);
            coord = cell2mat(a);
        end
        
        function s = sumnd(x)
            % sum over all dimensions
            s = sum(x(:));
        end
        
        function [val, coord] = maxnd(x)
            % [val, coord] = MaxIJ(x) returns the maximal value (val) in x and it's coordinated (coord)
            % unlike the built-in max x oppearates on all dimensions
            %%
            [val, i] = max(x(:));
            a = cell(1, ndims(x));
            [a{1:ndims(x)}] = ind2sub(size(x), i);
            coord = cell2mat(a);
        end
        
        function [val, coord] = minnd(x)
            % [val, coord] = MinIJ(x) returns the minimal value (val) in x and it's coordinated (coord)
            % unlike the built-in max x oppearates on all dimensions
            %%
            [val, i] = min(x(:));
            a = cell(1, ndims(x));
            [a{1:ndims(x)}] = ind2sub(size(x), i);
            coord = cell2mat(a);
        end
        
        function R = randi(p, n, m)
            % discrete random (use Matlab's)
            r = rand(n, m);
            P = cumsum(p);
            s = zeros(n, m);
            for i=1:length(p)-1
                s = s + (r > P(i));
            end
            R = s + 1;
        end
        
        function str = randstr(length)
            % str = RandStr(length) generate a random string of specified length
            alphabet = 'abcdefghijklmnopqrstuvwxyz0123456789';
            str = alphabet(randi(size(alphabet, 2), 1, length));
        end
        
        function varargout = findsegs(map, data)
            % Find segments in map. See (and use) the Segs class instead
            c = diff([0 map(:)' 0]);
            begf = find(c > 0);
            endf = find(c < 0) - 1;
            
            if nargin >= 2
                events = repmat(struct('data', [], 'beg', 0, 'end', 0), length(begf), 1);
                for i=1:length(begf)
                    if isvector(data)
                        events(i).data = data(begf(i):endf(i));
                    else
                        events(i).data = data(begf(i):endf(i), :);
                    end
                    events(i).beg = begf(i);
                    events(i).end = endf(i);
                end
                varargout{1} = events;
            else
                varargout{1} = begf;
                varargout{2} = endf;
            end
        end
        
        function plotsegs2(map, y, varargin)
            % Plot segments in map. See (and use) the Segs class instead
            [b,e] = Q.findsegs(map);
            for i=1:length(e)
                if nargin > 2
                    plot([b(i) e(i)], [y y], varargin{:});
                else
                    plot([b(i) e(i)], [y y], 'r-', 'LineWidth', 4);
                end
                Fig.Hon;
            end
            Fig.Hoff;
        end
        
        function plotsegs(map, c, factor)
            % Plot segments in map. See (and use) the Segs class instead
            if nargin < 3
                factor = 1;
            end
            [b, e] = Q.findsegs(map);
            ax = axis;
            y = ax(3:4);
            y(1) = y(1) + factor * (y(2)-y(1)) / 100;
            y(2) = y(2) - factor * (y(2)-y(1)) / 100;
            for i=1:length(e)
                patch([b(i) e(i) e(i) b(i)], [y(1) y(1) y(2) y(2)], 'w', 'FaceColor', c, 'EdgeColor', c, 'FaceAlpha', .2);
            end
        end
        
        function map = removegaps(map, maxgap)
            % Remove gaps in segments map. See (and use) the Segs class instead
            [b, e] = Q.findsegs(map);
            merge = (b(2:end) - e(1:end-1) - 1) <= maxgap;
            for i=find(merge)
                map(e(i)+1:b(i+1)-1) = true;
            end
        end
        
        function map = filtergaps(map, minsize)
            % Remove small segmenst in segments map. See (and use) the Segs class instead
            [b, e] = Q.findsegs(map);
            l = e - b - 1;
            b = b(l>minsize);
            e = e(l>minsize);
            map = Q.segstomap(length(map), b, e);
        end
        
        function [b,e] = mergesegs(map, merge)
            % Merge segmenst in segments map. See (and use) the Segs class instead
            %%
            [b,e] = Q.findsegs(map);
            idx = 1;
            valid = false(size(b));
            for i=1:length(merge)
                if idx == i
                    valid(idx) = true;
                end
                if merge(i)
                    e(idx) = e(i+1);
                else
                    idx = i+1;
                end
            end
            b = b(valid);
            e = e(valid);
        end
        
        function map = segstomap(len, b, e)
            % Create segments map. See (and use) the Segs class instead
            map = false(1, len);
            for i=1:length(b)
                map(b(i):e(i)) = true;
            end
        end
        
        function zipdepend(zipfile, mfile)
            % Put all dependent .m files in zip archive
            flist = matlab.codetools.requiredFilesAndProducts(mfile);
            pref = flist{1};
            for i=2:length(flist)
                other = flist{i};
                pref(end+1:length(other)) = ' ';
                other(end+1:length(pref)) = ' ';
                idx = find(pref ~= other, 1);
                if ~isempty(idx)
                    pref = pref(1:find(pref ~= other, 1)-1);
                end
            end
            for i=1:length(flist)
                flist{i} = flist{i}(length(pref)+1:end);
            end
            zip(zipfile,flist,pref);
        end
    end
end
