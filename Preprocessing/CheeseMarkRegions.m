function varargout = CheeseMarkRegions(varargin)
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

% Last Modified by GUIDE v2.5 14-Jan-2015 14:50:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @CheeseMarkRegions_OpeningFcn, ...
    'gui_OutputFcn',  @CheeseMarkRegions_OutputFcn, ...
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
function CheeseMarkRegions_OpeningFcn(hObject, eventdata, handles, varargin)
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

global CheeseMarkRegionsData;
global CheeseSquareData;

opt = struct();

CheeseMarkRegionsData.Polys = {};
CheeseMarkRegionsData.PolyLabels = {};

CheeseMarkRegionsData.Options = opt;
CheeseMarkRegionsData.Handles = handles;

axes(CheeseMarkRegionsData.Handles.ProgressBar);
imagesc(reshape(Colors.PrettyBlue, [1 1 3]));
set(CheeseMarkRegionsData.Handles.ProgressBar, 'XTick', [], 'YTick', [], 'XColor', Colors.PrettyBlue, 'YColor', Colors.PrettyBlue);

axes(CheeseMarkRegionsData.Handles.MainAxes);
CheeseMarkRegionsData.Image = label2rgb(Canvas.Checkers(400, [], 50), [1 1 1] * 206/255, [1 1 1]);
imagesc(CheeseMarkRegionsData.Image);
set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');
drawnow;

warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jframe=get(handles.figure1, 'javaframe');
jIcon=javax.swing.ImageIcon('.\mouse-square-small.gif');
jframe.setFigureIcon(jIcon);

CheeseMarkRegionsData.CurrID = 0;
CheeseMarkRegionsData.CurrentFrame = 1;
set(CheeseMarkRegionsData.Handles.VideoPosition, 'Min', 1, 'Max', CheeseSquareData.obj.Video.NumberOfFrames, 'value', 1, 'sliderstep', [1 1]/(CheeseSquareData.obj.Video.NumberOfFrames - 1));
f = fields(CheeseMarkRegionsData.Handles);
for i=1:length(f)
    try
        set(CheeseMarkRegionsData.Handles.(f{i}), 'KeyPressFcn', @figure1_KeyPressFcn);
    catch
    end
end

%%
if isempty(CheeseSquareData.obj.ROI)
    LoadFromPrototype('..\Tracking\Prototypes\Default.obj.mat');
    %LoadFromPrototype('Prototypes\withramps.obj.mat');
else
    CheeseSquareData.obj.ROI.TempPoints = CheeseSquareData.obj.ROI.Points;
    if ~isfield(CheeseSquareData.obj.ROI, 'IsAvail') || isempty(CheeseSquareData.obj.ROI.IsAvail)
        CheeseSquareData.obj.ROI.IsAvail = true(1, CheeseSquareData.obj.ROI.nRegions);
    end
end
%%
list = dir('Prototypes\*.obj.mat');
names = {list.name};
set(CheeseMarkRegionsData.Handles.ProtoList, 'String', names);
set(CheeseMarkRegionsData.Handles.ProtoList, 'value', 1:length(names));

%%
try
    GetFrame(CheeseSquareData.obj.Background.im);
catch
    GetFrame();
end

function LoadFromPrototype(src)
global CheeseSquareData;
global CheeseMarkRegionsData;
if isempty(CheeseSquareData.obj.ROI)
    if ischar(src)
        temp = load(src);
    else
        temp = src;
    end
    %%
    roi = temp.obj.ROI;
    if ~isfield(roi, 'Points') || isempty(roi.Points)
        %%
        CheeseSquareData.obj.ROI.Points = cell(1, length(roi.Regions));
        for i=1:length(roi.Regions)
            p = mask2poly(roi.Regions{i}, 4);
            CheeseSquareData.obj.ROI.Points{i} = p(1:end-1, :);
        end
    else
        CheeseSquareData.obj.ROI.Points = roi.Points;
    end
    CheeseSquareData.obj.ROI.RegionNames = roi.RegionNames;
    CheeseSquareData.obj.ROI.IsSheltered = roi.IsSheltered;
    CheeseSquareData.obj.ROI.nHidden = sum(roi.IsSheltered);
    CheeseSquareData.obj.ROI.nZones = length(roi.RegionNames) + CheeseSquareData.obj.ROI.nHidden;
    CheeseSquareData.obj.ROI.nRegions = length(roi.RegionNames);
    CheeseSquareData.obj.ROI.IsAvail = true(1, CheeseSquareData.obj.ROI.nRegions);
end
CheeseMarkRegionsData.Polys = {};
CheeseMarkRegionsData.PolyLabels = {};

if ~isfield(CheeseSquareData.obj.ROI, 'TempPoints') || isempty(CheeseSquareData.obj.ROI.TempPoints)
    CheeseSquareData.obj.ROI.TempPoints = CheeseSquareData.obj.ROI.Points;
end

function handles = AddToRegions(handles, fig, form, label, value)
style = {'FontName', 'Segoe UI', 'FontSize', 9.0, 'BackgroundColor', 'w'};
idx = length(get(handles.GroupArea, 'children'));
pos = get(handles.GroupArea, 'position');
handles.([label 'label']) = uicontrol(handles.GroupArea, ...
    'Style','text',...
    'String', label,...
    'Position', [10 pos(4) - (45 + idx * 15), pos(3)/2-20, 17], ...
    'HorizontalAlignment', 'left', ...
    style{:},...
    'Callback',@p_Callback);

handles.(label) = uicontrol(handles.GroupArea, ...
    'Tag', label, ...
    'Style','checkbox',...
    'Value', value,...
    'Position', [pos(3)/2 + 10 pos(4) - (48 + idx * 15), pos(3)/2-25, 21], ...
    'Callback', @Option_Callback, ...
    style{:});

guidata(handles.figure1, handles);

% --------------------------------------------------------------------
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
global CheeseMarkRegionsData;
mb = CheeseMarkRegionsData.Handles.MessageBox;
prev = get(mb, 'string');
if isempty(prev)
    set(mb, 'string', {[' - ' str]});
else
    set(mb, 'string', {prev{:}, [' - ' str]});
end
drawnow

function ProgressBar(str, p)
global CheeseMarkRegionsData;
mb = CheeseMarkRegionsData.Handles.MessageLine;
set(mb, 'string', [' ' str]);

if nargin < 2
    p = 1;
end
axes(CheeseMarkRegionsData.Handles.ProgressBar);
z1 = ones(1, round(50 * p));
z2 = ones(1, round(50 * (1-p)));
c1 = Colors.PrettyBlue;
c2 = Colors.PrettyRed;
imagesc([cat(3, z1*c1(1), z1*c1(2), z1*c1(3)) cat(3, z2*c2(1), z2*c2(2), z2*c2(3))]);
axis off;
axes(CheeseMarkRegionsData.Handles.MainAxes);


function Option_Callback(hObject, eventdata, handles)
global CheeseMarkRegionsData;
tag = get(hObject, 'tag');
CheeseMarkRegionsData.Options.(tag) = str2num(get(hObject, 'string'));

function Start()
global CheeseMarkRegionsData;
set(CheeseMarkRegionsData.Handles.RunButton, 'String', 'Running...');
mb = CheeseMarkRegionsData.Handles.MessageBox;
set(mb, 'string', '');

% CheeseMarkRegionsData.GUI.nDots = 0;
% CheeseMarkRegionsData.GUI.Timer = timer(...
%     'TimerFcn', @IRun,... 
%     'ExecutionMode', 'fixedSpacing', ...
%     'Period', .5);
% start(CheeseMarkRegionsData.GUI.Timer)

function Finish()
global CheeseMarkRegionsData;
CheeseMarkRegionsData.GUI.Timer.stop;
set(CheeseMarkRegionsData.Handles.RunButton, 'String', 'Run');

function Success()
global CheeseMarkRegionsData;
set(CheeseMarkRegionsData.Handles.NextButton, 'Enable', 'on');

% --- Outputs from this function are returned to the command line.
function varargout = CheeseMarkRegions_OutputFcn(hObject, eventdata, handles)
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
global CheeseMarkRegionsData;
if length(eventdata.Character) > 1 || isempty(eventdata.Character) || strcmp(eventdata.Character, '')
    return;
end

num = eventdata.Character - '0';
if num >= 0 && num <= 9
    CheeseMarkRegionsData.CurrID = num;
    if CheeseMarkRegionsData.Points.isKey(num)
        delete(CheeseMarkRegionsData.Points(num));
    end
    h = impoint(CheeseMarkRegionsData.Handles.MainAxes);
    CheeseMarkRegionsData.Points(num) = h;
    h.setColor('g');
    h.setString(num);
end

% --------------------------------------------------------------------
function GetFrame(frame)
global CheeseSquareData;
global CheeseMarkRegionsData;

CheeseMarkRegionsData.CurrID = 0;

ProgressBar('loading frame...');
if nargin < 1
    frame = randi(CheeseSquareData.obj.Video.NumberOfFrames, 1, 1);
end
if isscalar(frame)
    CheeseSquareData.obj.Video.FrameNumber = frame;
    CheeseMarkRegionsData.CurrentFrame = frame;
    img = CheeseSquareData.obj.Video.CurrentFrame;
    set(CheeseMarkRegionsData.Handles.VideoPosition, 'Value', frame);
else
    img = frame;
end
CheeseMarkRegionsData.Image = img;
set(CheeseMarkRegionsData.Handles.FrameNumTB, 'string', sprintf('%d', CheeseSquareData.obj.Video.FrameNumber));
ProgressBar('');


%%
if get(CheeseMarkRegionsData.Handles.EnhanceImage, 'value')
    imagesc(CV.HistEq(img));
else
    imagesc(img);
end
set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');
cmap = repmat(Colormaps.Categorical, [3, 1]);
for r=find(CheeseSquareData.obj.ROI.IsAvail)
    if length(CheeseMarkRegionsData.Polys) >= r
        %delete(CheeseMarkRegionsData.Polys{r});
    end
    h = impoly(gca, CheeseSquareData.obj.ROI.TempPoints{r});
    CheeseMarkRegionsData.Polys{r} = h;
    
    h.setColor(cmap(r, :));
    h.addNewPositionCallback(@(pos) AreaMove(r, pos));
end

RefreshLabels;

%%
ProgressBar('');

function RefreshLabels(R)
global CheeseMarkRegionsData;
global CheeseSquareData;
if nargin < 1
    %R = 1:length(CheeseMarkRegionsData.Polys);
    R = find(CheeseSquareData.obj.ROI.IsAvail);
end

cmap = repmat(Colormaps.Categorical, [3, 1]);
for r=R
    if length(CheeseMarkRegionsData.PolyLabels) >= r
        try
            delete(CheeseMarkRegionsData.PolyLabels{r});
        catch
        end
        CheeseMarkRegionsData.PolyLabels{r} = [];
    end
end
for r=R
    %%
    pt = CheeseSquareData.obj.ROI.TempPoints{r};
    inside = Q.inrange(pt(:, 1), 0, CheeseSquareData.obj.Video.Width) & Q.inrange(pt(:, 2), 0, CheeseSquareData.obj.Video.Height);
    %%
    [~, idx] = max(sum(pt, 2) .* inside);
    pos = CheeseMarkRegionsData.Polys{r}.getPosition;
    str = sprintf('%s', CheeseSquareData.obj.ROI.RegionNames{r});
    CheeseMarkRegionsData.PolyLabels{r} = text(pos(idx, 1)+10, pos(idx, 2)+10, str, 'BackgroundColor', cmap(r, :), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'color', 'w');
end
data = cell(CheeseSquareData.obj.ROI.nRegions, 2);
for i=1:CheeseSquareData.obj.ROI.nRegions
    data{i, 1} = CheeseSquareData.obj.ROI.IsSheltered(i);
    data{i, 2} = CheeseSquareData.obj.ROI.IsAvail(i);
    %AddToRegions(handles, handles.figure1, handles.Form, name, CheeseSquareData.obj.ROI.IsSheltered(i));
end
set(CheeseMarkRegionsData.Handles.RegionList, 'ColumnName', {'Is nest', 'Available'}, 'ColumnFormat', {'logical', 'logical'}, 'ColumnEditable', [true true true], 'RowName', CheeseSquareData.obj.ROI.RegionNames, 'data', data);




function AreaMove(id, pos)
global CheeseSquareData;
CheeseSquareData.obj.ROI.TempPoints{id} = pos;
RefreshLabels(id)


% --- Executes on mouse press over axes background.
function MainAxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to MainAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
p = get(gca, 'CurrentPoint');
x = p(1, 1)
y = p(1, 2)

function UpdateObject()
global CheeseSquareData;
global CheeseMarkRegionsData;
roi = CheeseSquareData.obj.ROI;
if ~isfield(roi, 'TempPoints')
    return;
end

roi.Points = roi.TempPoints;
roi.ZoneNames{1} = 'Open';
idx = 2;
for r=find(roi.IsAvail)
    roi.Regions{r} = CheeseMarkRegionsData.Polys{r}.createMask;
    roi.ZoneNames{idx} = roi.RegionNames{r};
    idx = idx + 1;
    if roi.IsSheltered(r)
        roi.ZoneNames{idx} = ['(' roi.RegionNames{r} ')'];
        idx = idx + 1;
    end
end
roi.nZones = length(roi.ZoneNames);
roi.Hidden = {roi.Regions{roi.IsSheltered & roi.IsAvail}};
roi.HiddenImage = max(cat(3, roi.Regions{roi.IsSheltered & roi.IsAvail}), [], 3);
roi.nHidden = sum(roi.IsSheltered & roi.IsAvail);
for r=1:roi.nHidden
    rp = regionprops(roi.Hidden{r}, 'Centroid');
    roi.HiddenCenters(r, :) = rp.Centroid;
end
roi = rmfield(roi, 'TempPoints');
CheeseSquareData.obj.ROI = roi;

% --- Executes on button press in NextButton.
function NextButton_Callback(hObject, eventdata, handles)
% hObject    handle to NextButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseMarkRegionsData;
UpdateObject()
close(CheeseMarkRegionsData.Handles.figure1);
CheeseSave;

% --- Executes on button press in PrevButton.
function PrevButton_Callback(hObject, eventdata, handles)
% hObject    handle to PrevButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UpdateObject()
close(gcf);
CheeseMarkColors

% --- Executes on button press in Background.
function Background_Callback(hObject, eventdata, handles)
% hObject    handle to Background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseSquareData;
GetFrame(CheeseSquareData.obj.Background.im);


% --- Executes on button press in RandomButton.
function RandomButton_Callback(hObject, eventdata, handles)
% hObject    handle to RandomButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GetFrame();


% --- Executes on button press in EnhanceImage.
function EnhanceImage_Callback(hObject, eventdata, handles)
% hObject    handle to EnhanceImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EnhanceImage
global CheeseMarkRegionsData;
GetFrame(CheeseMarkRegionsData.Image);


% --- Executes on button press in FromFile.
function FromFile_Callback(hObject, eventdata, handles)
% hObject    handle to FromFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global CheeseSquareData;
[filename, pathname] = uigetfile('*.obj.mat');
if filename ~= 0
    try
        LoadFromPrototype(fullfile(pathname, filename));
        %%
        try
            GetFrame(CheeseSquareData.obj.Background.im);
        catch
            GetFrame();
        end
    catch
        msgbox('unable to load regions-of-interest from file');
    end
end



function ProtoList_Callback(hObject, eventdata, handles)
% hObject    handle to ProtoList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ProtoList as text
%        str2double(get(hObject,'String')) returns contents of ProtoList as a double
global CheeseMarkRegionsData;
global CheeseSquareData;
if strcmp(get(CheeseMarkRegionsData.Handles.figure1,'SelectionType'), 'open')
    CheeseSquareData.obj.ROI = [];
    LoadFromPrototype(['Prototypes\' eventdata.Source.String{eventdata.Source.Value}]);
    %%
    try
        GetFrame(CheeseSquareData.obj.Background.im);
    catch
        GetFrame();
    end
end

% --- Executes during object creation, after setting all properties.
function ProtoList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ProtoList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AddRegion.
function AddRegion_Callback(hObject, eventdata, handles)
% hObject    handle to AddRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseSquareData;
answer = inputdlg({'Enter region name:'}, 'Input', 1);
answer = regexprep(answer, '[ ]*', '');
answer = regexprep(answer, '^\d*', '');
h = CheeseSquareData.obj.Video.Height;
w = CheeseSquareData.obj.Video.Width;
if ~isempty(answer) && ~isempty(answer{1})
    CheeseSquareData.obj.ROI.RegionNames{end+1} = answer{1};
    CheeseSquareData.obj.ROI.nRegions = length(CheeseSquareData.obj.ROI.RegionNames);
    CheeseSquareData.obj.ROI.nZones = CheeseSquareData.obj.ROI.nZones + 1;
    CheeseSquareData.obj.ROI.IsSheltered(end+1) = false;
    CheeseSquareData.obj.ROI.Points{CheeseSquareData.obj.ROI.nRegions} = [h/3, w/3; h/3 w*2/3; h*2/3 w*2/3; h*2/3 w/3];
    CheeseSquareData.obj.ROI.TempPoints{CheeseSquareData.obj.ROI.nRegions} = [h/3, w/3; h/3 w*2/3; h*2/3 w*2/3; h*2/3 w/3];
    CheeseSquareData.obj.ROI.IsAvail(end+1) = true;
end
GetFrame(CheeseSquareData.obj.Video.CurrentFrame);


% --- Executes when entered data in editable cell(s) in RegionList.
function RegionList_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to RegionList (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global CheeseSquareData;
idx = eventdata.Indices(1) ;
if eventdata.Indices(2) == 2
    CheeseSquareData.obj.ROI.IsAvail(idx) = eventdata.NewData ~= 0;
    GetFrame(CheeseSquareData.obj.Video.CurrentFrame);
else
    CheeseSquareData.obj.ROI.IsSheltered(idx) = eventdata.NewData ~= 0;
end



function FrameNumTB_Callback(hObject, eventdata, handles)
% hObject    handle to FrameNumTB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrameNumTB as text
%        str2double(get(hObject,'String')) returns contents of FrameNumTB as a double
global CheeseMarkRegionsData;
GetFrame(str2double(get(CheeseMarkRegionsData.Handles.FrameNumTB, 'string')));

% --- Executes during object creation, after setting all properties.
function FrameNumTB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrameNumTB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
