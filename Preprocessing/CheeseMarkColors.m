function varargout = CheeseMarkColors(varargin)
% CHEESESQRFINDBKG MATLAB code for CheeseSqrFindBkg.fig
%      CHEESESQRFINDBKG, by itself, creates a new CHEESESQRFINDBKG or raises the existing
%      singleton*.
%
%      H = CHEESESQRFINDBKG returns the handle to a new CHEESESQRFINDBKG or the handle to
%      the existing singleton*.
%
%      CHEESESQRFINDBKG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHEESESQRFINDBKG.M with the given input arguments.
%
%      CHEESESQRFINDBKG('Property','Value',...) creates a new CHEESESQRFINDBKG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CheeseSqrFindBkg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CheeseSqrFindBkg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CheeseSqrFindBkg

% Last Modified by GUIDE v2.5 25-Dec-2014 14:10:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @CheeseMarkColors_OpeningFcn, ...
    'gui_OutputFcn',  @CheeseMarkColors_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before CheeseSqrFindBkg is made visible.
function CheeseMarkColors_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CheeseSqrFindBkg (see VARARGIN)

% Choose default command line output for CheeseSqrFindBkg
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CheeseSqrFindBkg wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global CheeseMarkColorsData;
global CheeseSquareData;

opt = struct();
opt.nColorBins = 20;

CheeseMarkColorsData.Options = opt;
CheeseMarkColorsData.Handles = handles;
CheeseMarkColorsData.Points = containers.Map('KeyType', 'int32', 'ValueType', 'any');
if ~isfield(CheeseSquareData.obj.Colors, 'Marks') || isempty(CheeseSquareData.obj.Colors.Marks)
    CheeseMarkColorsData.IsChanged = false;
    CheeseSquareData.obj.Colors.Marks = struct('x', [], 'y', [], 'id', [], 'frame', [], 'color', []);
else
    CheeseMarkColorsData.IsChanged = true;
end
CheeseMarkColorsData.Preview = zeros(0, 4, 3);

axes(CheeseMarkColorsData.Handles.ProgressBar);
imagesc(reshape(Colors.PrettyBlue, [1 1 3]));
%axis off;
set(CheeseMarkColorsData.Handles.ProgressBar, 'XTick', [], 'YTick', [], 'XColor', Colors.PrettyBlue, 'YColor', Colors.PrettyBlue);

axes(CheeseMarkColorsData.Handles.MainAxes);
CheeseMarkColorsData.Image = label2rgb(Canvas.Checkers(400, [], 50), [1 1 1] * 206/255, [1 1 1]);
CheeseMarkColorsData.Map = false(size(CheeseMarkColorsData.Image, 1), size(CheeseMarkColorsData.Image, 2));
imagesc(CheeseMarkColorsData.Image);
set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');
drawnow;

warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jframe=get(handles.figure1, 'javaframe');
jIcon=javax.swing.ImageIcon('.\mouse-square-small.gif');
jframe.setFigureIcon(jIcon);

CheeseMarkColorsData.CurrID = 0;

CheeseMarkColorsData.CurrentFrame = 1;
set(CheeseMarkColorsData.Handles.VideoPosition, 'Min', 1, 'Max', CheeseSquareData.obj.Video.NumberOfFrames, 'value', 1, 'sliderstep', [1 1]/(CheeseSquareData.obj.Video.NumberOfFrames - 1));

f = fields(CheeseMarkColorsData.Handles);
for i=1:length(f)
    try
        set(CheeseMarkColorsData.Handles.(f{i}), 'KeyPressFcn', @figure1_KeyPressFcn);
    catch
    end
end

NextFrame_Callback();

function handles = AddToForm(handles, fig, form, name, label, value)
style = {'FontName', 'Segoe UI', 'FontSize', 9.0, 'BackgroundColor', 'w'};
idx = length(get(handles.Form, 'children'));
pos = get(handles.Form, 'position');
handles.([name 'label']) = uicontrol(handles.Form, ...
    'Style','text',...
    'String', label,...
    'Position', [10 pos(4) - (45 + idx * 15), pos(3)/2-20, 17], ...
    'HorizontalAlignment', 'left', ...
    style{:},...
    'Callback',@p_Callback);

handles.(name) = uicontrol(handles.Form, ...
    'Tag', name, ...
    'Style','edit',...
    'String', num2str(value),...
    'Position', [pos(3)/2 + 10 pos(4) - (48 + idx * 15), pos(3)/2-25, 21], ...
    'Callback', @Option_Callback, ...
    style{:});

guidata(handles.figure1, handles);


function Message(str)
%%
global CheeseMarkColorsData;
mb = CheeseMarkColorsData.Handles.MessageBox;
prev = get(mb, 'string');
if isempty(prev)
    set(mb, 'string', {[' - ' str]});
else
    set(mb, 'string', {prev{:}, [' - ' str]});
end
%% scroll to bottom
try
    jhEdit = findjobj(mb);
    jEdit = jhEdit.getComponent(0).getComponent(0);
    jEdit.setCaretPosition(jEdit.getDocument.getLength);
catch
end

%%


drawnow

function ProgressBarWarning(str, p)
global CheeseMarkColorsData;
mb = CheeseMarkColorsData.Handles.MessageLine;
set(mb, 'string', [' ' str]);
set(mb, 'BackgroundColor', [0.8627    0.2863    0.3490]);
set(mb, 'ForegroundColor', [0.8627    0.2863    0.3490]);
if nargin < 2
    p = 1;
end
axes(CheeseMarkColorsData.Handles.ProgressBar);
z1 = ones(1, round(50 * p));
z2 = ones(1, round(50 * (1-p)));
c1 = Colors.PrettyBlue;
c2 = Colors.PrettyRed;
imagesc([cat(3, z1*c1(1), z1*c1(2), z1*c1(3)) cat(3, z2*c2(1), z2*c2(2), z2*c2(3))]);
axis off;
axes(CheeseMarkColorsData.Handles.MainAxes);

function ProgressBar(str, p)
global CheeseMarkColorsData;
mb = CheeseMarkColorsData.Handles.MessageLine;
set(mb, 'string', [' ' str]);
set(mb, 'BackgroundColor', 'w');
if nargin < 2
    p = 1;
end
axes(CheeseMarkColorsData.Handles.ProgressBar);
z1 = ones(1, round(50 * p));
z2 = ones(1, round(50 * (1-p)));
c1 = Colors.PrettyBlue;
c2 = Colors.PrettyRed;
imagesc([cat(3, z1*c1(1), z1*c1(2), z1*c1(3)) cat(3, z2*c2(1), z2*c2(2), z2*c2(3))]);
axis off;
axes(CheeseMarkColorsData.Handles.MainAxes);


function Option_Callback(hObject, eventdata, handles)
global CheeseMarkColorsData;
tag = get(hObject, 'tag');
CheeseMarkColorsData.Options.(tag) = str2num(get(hObject, 'string'));

function Start()
global CheeseMarkColorsData;
set(CheeseMarkColorsData.Handles.RunButton, 'String', 'Running...');
mb = CheeseMarkColorsData.Handles.MessageBox;
set(mb, 'string', '');

% CheeseMarkColorsData.GUI.nDots = 0;
% CheeseMarkColorsData.GUI.Timer = timer(...
%     'TimerFcn', @IRun,... 
%     'ExecutionMode', 'fixedSpacing', ...
%     'Period', .5);
% start(CheeseMarkColorsData.GUI.Timer)

function IRun(mTimer, ~)
global CheeseMarkColorsData;
CheeseMarkColorsData.GUI.nDots = mod(CheeseMarkColorsData.GUI.nDots + 1, 3);
n = CheeseMarkColorsData.GUI.nDots+1;
set(CheeseMarkColorsData.Handles.RunButton, 'String', ['Running' repmat('.', 1, n) repmat(' ', 1, 3-n)]);
%dbstack

function Finish()
global CheeseMarkColorsData;
CheeseMarkColorsData.GUI.Timer.stop;
set(CheeseMarkColorsData.Handles.RunButton, 'String', 'Run');

function Success()
global CheeseMarkColorsData;
set(CheeseMarkColorsData.Handles.NextButton, 'Enable', 'on');

% --- Outputs from this function are returned to the command line.
function varargout = CheeseMarkColors_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function MessageBox_Callback(hObject, eventdata, handles)
% hObject    handle to MessageBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MessageBox as text
%        str2double(get(hObject,'String')) returns contents of MessageBox as a double


% --- Executes during object creation, after setting all properties.
function MessageBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MessageBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MessageLine_Callback(hObject, eventdata, handles)
% hObject    handle to MessageLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MessageLine as text
%        str2double(get(hObject,'String')) returns contents of MessageLine as a double


% --- Executes during object creation, after setting all properties.
function MessageLine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MessageLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on slider movement.
function VideoPosition_Callback(hObject, eventdata, handles)
% hObject    handle to VideoPosition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

GetFrame(round(get(hObject,'Value')));

% --- Executes during object creation, after setting all properties.
function VideoPosition_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VideoPosition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
global CheeseMarkColorsData;
if length(eventdata.Character) > 1 || isempty(eventdata.Character) || strcmp(eventdata.Character, '')
    return;
end

num = eventdata.Character - '0';
if num >= 0 && num <= 9
    CheeseMarkColorsData.CurrID = num;
    if CheeseMarkColorsData.Points.isKey(num)
        delete(CheeseMarkColorsData.Points(num));
    end
    h = impoint(CheeseMarkColorsData.Handles.MainAxes);
    CheeseMarkColorsData.Points(num) = h;
    h.setColor('g');
    h.setString(num);
end
if eventdata.Character == 'n' || eventdata.Character == 'N' || eventdata.Character == ' '
    NextFrame_Callback();
end

% --------------------------------------------------------------------
function ShowHist()


% --------------------------------------------------------------------
function GetFrame(num)
global CheeseSquareData;
global CheeseMarkColorsData;

ProgressBar('loading frame...');
%%
CheeseMarkColorsData.CurrID = 0;
if ~isfield(CheeseSquareData.obj.Background, 'im');
    ProgressBarWarning('background is not defined!');
    return;
end

points = CheeseMarkColorsData.Points.keys;
hsv = rgb2hsv(CheeseMarkColorsData.Image);
lab = CV.RGB2NormLAB(CheeseMarkColorsData.Image);
label = bwlabel(CheeseMarkColorsData.Map);
subjmap = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
for i=[points{:}]
    if ~isempty(CheeseMarkColorsData.Points(i))
        pos = CheeseMarkColorsData.Points(i).getPosition;
        pt = CheeseMarkColorsData.Points(i);
        x = round(pos(1));
        y = round(pos(2));
        l = label(y, x);
        if l == 0
            Message(sprintf('ignoring mark for %d; out of region', i));
        else
            subjmap(i) = l;
        end
    end
end
%%
values = subjmap.values;
if length(unique([values{:}])) < subjmap.Count
    keys = subjmap.keys;
    for v=unique([values{:}])
        map = [values{:}] == v;
        if sum(map) > 1
            for i=1:length(subjmap.keys)
                k = keys{i};
                if subjmap(k) == v
                    Message(sprintf('ignoring mark for %d; multiple subjects', k));
                    subjmap.remove(k);
                end
            end
            %Message('ignoring ');
        end
    end
end

%%
for i=[points{:}]
    if ~isempty(CheeseMarkColorsData.Points(i))
        try
            %%
            pos = CheeseMarkColorsData.Points(i).getPosition;
            pt = CheeseMarkColorsData.Points(i);
            pt.delete
            %%
            idx = length(CheeseSquareData.obj.Colors.Marks.id) + 1;
            %%
            CheeseSquareData.obj.Colors.Marks.id(idx) = -1;
            x = round(pos(1));
            y = round(pos(2));

            l = subjmap(i);
            if l == 0
                continue;
            end
            curr = label == l;
            CheeseSquareData.obj.Colors.Marks.x(idx) = x;
            CheeseSquareData.obj.Colors.Marks.y(idx) = y;
            CheeseSquareData.obj.Colors.Marks.frame(idx) = CheeseMarkColorsData.CurrentFrame;
            CheeseSquareData.obj.Colors.Marks.color(idx, :) = squeeze(CheeseMarkColorsData.Image(y, x, :))';
            CheeseSquareData.obj.Colors.Marks.size(idx) = sum(curr(:));
            
            CheeseSquareData.obj.Colors.Bins = linspace(0, 1, CheeseMarkColorsData.Options.nColorBins);
            %
            chnames = {'r', 'g', 'b'};
            for c=1:3
                ch = CheeseMarkColorsData.Image(:, :, c);
                CheeseSquareData.obj.Colors.Marks.hist.(chnames{c})(idx, :) = histc(im2double(ch(curr)), CheeseSquareData.obj.Colors.Bins);
            end
            %
            chnames = {'h', 's', 'v'};
            for c=1:3
                ch = hsv(:, :, c);
                CheeseSquareData.obj.Colors.Marks.hist.(chnames{c})(idx, :) = histc(im2double(ch(curr)), CheeseSquareData.obj.Colors.Bins);
            end
            %
%             chnames = {'ll', 'aa', 'bb'};
%             for c=1:3
%                 ch = lab(:, :, c);
%                 CheeseSquareData.obj.Colors.Marks.hist.(chnames{c})(idx, :) = histc(im2double(ch(curr)), CheeseSquareData.obj.Colors.Bins);
%             end
            %
            CheeseSquareData.obj.Colors.Marks.id(idx) = i;
            CheeseMarkColorsData.IsChanged = true;
        catch
        end
        CheeseMarkColorsData.Points.remove(i);
    end
end
CheeseMarkColorsData.Points.remove(CheeseMarkColorsData.Points.keys);

%%
n = max(max(CheeseSquareData.obj.Colors.Marks.id), 4);
count = max(histc(CheeseSquareData.obj.Colors.Marks.id, 0:n));
CheeseMarkColorsData.Preview = nan(n+1, count, 3);
for i=0:n
    idx = 1;
    for j=find(CheeseSquareData.obj.Colors.Marks.id == i)
        CheeseMarkColorsData.Preview(i+1, idx, :) = CheeseSquareData.obj.Colors.Marks.color(j, :);
        idx = idx + 1;
    end
end

%%
axes(CheeseMarkColorsData.Handles.ColorTable);
im = CheeseMarkColorsData.Preview;
im(isnan(im)) = 255;
imagesc((im)/255);
%set(gca, 'Ydir', 'normal');
xtick = unique([count 5:5:25]);
set(CheeseMarkColorsData.Handles.ColorTable, 'yTick', 1:n+1, 'yTicklabel', 0:n, 'xTick', xtick + .5, 'xticklabel', xtick);
set(gca, 'xGrid', 'on')
Fig.Hon
for i=0:double(n)
    c = sum(CheeseSquareData.obj.Colors.Marks.id == i);
    text(count+.5, i+1, [' ' num2str(c)], 'horizontalalignment', 'left');
end
Fig.Hoff
Fig.Fix;
axes(CheeseMarkColorsData.Handles.MainAxes);
drawnow;

%%
if nargin == 0
    map = [];
    while ~any(map(:))
        num = randi(round(CheeseSquareData.obj.Video.NumberOfFrames), 1, 1);
        [map, img] = ImageSegmentClassic(CheeseSquareData.obj, num, struct('Output', false));
%         [map, img] = ImageSegmentClassicBlackFur(CheeseSquareData.obj, num, struct('Output', false));
    end
else
    [map, img] = ImageSegmentClassic(CheeseSquareData.obj, num, struct('Output', false));
end
CheeseMarkColorsData.CurrentFrame = num;
CheeseMarkColorsData.Image = img;
CheeseMarkColorsData.Map = map;

%%

b = bwboundaries(map);
imagesc(img);
Fig.Hon
for i=1:length(b)
    p = bwtraceboundary(map, b{i}(1, :), 'N');
    patch(p(:, 2), p(:, 1), 'y', 'EdgeColor', [.5 .5 .5], 'FaceColor', 'y', 'FaceAlpha', 0);
end
Fig.Hoff

set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');

ProgressBar('');

% --- Executes on button press in NextFrame.
function NextFrame_Callback(hObject, eventdata, handles)
% hObject    handle to NextFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseMarkColorsData;

GetFrame();
set(CheeseMarkColorsData.Handles.VideoPosition, 'Value', CheeseMarkColorsData.CurrentFrame);

% --- Executes on mouse press over axes background.
function MainAxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to MainAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
p = get(gca, 'CurrentPoint');
x = p(1, 1)
y = p(1, 2)

% --- Executes on button press in Clear.
function Clear_Callback(hObject, eventdata, handles)
% hObject    handle to Clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseMarkColorsData;
global CheeseSquareData;
choice = questdlg('This action will wipe out all previous marks', 'Are you sure?', 'Yes','No');
switch choice
    case 'Yes'
        CheeseSquareData.obj.Colors.Marks = struct('x', [], 'y', [], 'id', [], 'frame', [], 'color', []);
        CheeseMarkColorsData.IsChanged = true;
    case 'no'
        
end
GetFrame(CheeseMarkColorsData.CurrentFrame);

function UpdateObject()
global CheeseSquareData;
global CheeseMarkColorsData;
nBins = CheeseMarkColorsData.Options.nColorBins;
colors = CheeseSquareData.obj.Colors;
marks = colors.Marks;
if isempty(CheeseSquareData.obj.nSubjects)
    CheeseSquareData.obj.nSubjects = max(marks.id);
else
    CheeseSquareData.obj.nSubjects = max(max(marks.id), CheeseSquareData.obj.nSubjects);
end

%%
chnames = {'r', 'g', 'b', 'h', 's', 'v'};
for i=1:length(chnames)
    colors.Histrogram.(upper(chnames{i})) = zeros(CheeseSquareData.obj.nSubjects, nBins);
    colors.Background.(upper(chnames{i})) = zeros(1, nBins);
    colors.Histrogram.Avg.(upper(chnames{i})) = zeros(CheeseSquareData.obj.nSubjects, nBins);
    colors.Background.Avg.(upper(chnames{i})) = zeros(1, nBins);
end
colors.Histrogram.Count = zeros(1, CheeseSquareData.obj.nSubjects);
colors.Background.Count = 0;
%
for i=1:length(marks.id)
    tid = marks.id(i);
    if tid > 0
        type = 'Histrogram';
        id = tid;
    elseif tid == 0
        type = 'Background';
        id = 1;
    else
        continue;
    end
    colors.(type).Count(id) = colors.(type).Count(id) + 1;
    for c=1:length(chnames)
        ch = chnames{c};
        colors.(type).(upper(ch))(id, :) = colors.(type).(upper(ch))(id, :) + marks.hist.(ch)(i, :);
        colors.(type).Avg.(upper(ch))(id, :) = colors.(type).Avg.(upper(ch))(id, :) + marks.hist.(ch)(i, :) / sum(marks.hist.(ch)(i, :));
    end
end
CheeseSquareData.obj.Colors = colors;


% --- Executes on button press in NextButton.
function NextButton_Callback(hObject, eventdata, handles)
% hObject    handle to NextButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UpdateObject
close(gcf);
CheeseMarkRegions

% --- Executes on button press in PrevButton.
function PrevButton_Callback(hObject, eventdata, handles)
% hObject    handle to PrevButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UpdateObject
close(gcf);
CheeseBackground



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FromFile.
function FromFile_Callback(hObject, eventdata, handles)
% hObject    handle to FromFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseSquareData;
[filename, pathname] = uigetfile('*.obj.mat');
if filename ~= 0
    try
        src = load(fullfile(pathname, filename));
        if isfield(src.obj.Colors, {'Marks'})
            CheeseSquareData.obj.Colors.Marks = src.obj.Colors.Marks;
        else
            error(MException('failed'));
        end
    catch
        msgbox('unable to load colors from file');
    end
end


% --- Executes on button press in ReapButton.
function ReapButton_Callback(hObject, eventdata, handles)
% hObject    handle to ReapButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseSquareData;
global CheeseMarkColorsData;

if isempty(CheeseSquareData.obj.Colors.Marks.frame)
    warndlg({'no frame information' 'colors probably originated from a different file'});
else
    choice = questdlg('This action should only follow a change in the algorithm', 'Are you sure you know what you are doing?', 'Yes','No');
    switch choice
        case 'Yes'
        case 'No'
            return;
        otherwise
            return;
    end
    
    m = CheeseSquareData.obj.Colors.Marks;
    CheeseSquareData.obj.Colors.Marks = struct('x', [], 'y', [], 'id', [], 'frame', [], 'color', []);
    CheeseMarkColorsData.IsChanged = true;
    %%
    frames = unique(m.frame);
    idx = 1;
    for f=frames
        GetFrame(f);
        map = m.frame == f;
        for i=find(map)
            h = impoint(CheeseMarkColorsData.Handles.MainAxes, m.x(i), m.y(i));
            CheeseMarkColorsData.Points(m.id(i)) = h;
            h.setColor('g');
            h.setString(m.id(i));
%         for i=1:
%             CheeseMarkColorsData.CurrID = num;
%             if CheeseMarkColorsData.Points.isKey(num)
%                 delete(CheeseMarkColorsData.Points(num));
%             end
%             h = impoint(CheeseMarkColorsData.Handles.MainAxes);
%             CheeseMarkColorsData.Points(num) = h;
%             h.setColor('g');
%             h.setString(num);
%         end
            %fprintf('press any key (%d)...\n', idx);
            %pause
            %idx = idx + 1;
        end
        
    end
    GetFrame(f);
end