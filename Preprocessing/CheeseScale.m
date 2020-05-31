function varargout = CheeseScale(varargin)
% CHEESESCALE MATLAB code for CheeseScale.fig
%      CHEESESCALE, by itself, creates a new CHEESESCALE or raises the existing
%      singleton*.
%
%      H = CHEESESCALE returns the handle to a new CHEESESCALE or the handle to
%      the existing singleton*.
%
%      CHEESESCALE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHEESESCALE.M with the given input arguments.
%
%      CHEESESCALE('Property','Value',...) creates a new CHEESESCALE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CheeseScale_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CheeseScale_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CheeseScale

% Last Modified by GUIDE v2.5 14-Jan-2015 14:42:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CheeseScale_OpeningFcn, ...
                   'gui_OutputFcn',  @CheeseScale_OutputFcn, ...
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


% --- Executes just before CheeseScale is made visible.
function CheeseScale_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CheeseScale (see VARARGIN)

% Choose default command line output for CheeseScale
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CheeseScale wait for user response (see UIRESUME)
% uiwait(handles.figure1);

warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jframe=get(handles.figure1, 'javaframe');
jIcon=javax.swing.ImageIcon('.\mouse-square-small.gif');
jframe.setFigureIcon(jIcon);

global CheeseScaleData;
global CheeseSquareData;

CheeseScaleData.Handles = handles;
CheeseScaleData.IsPlaying = false;
set(CheeseScaleData.Handles.PlayStop, 'String', '>');

axes(CheeseScaleData.Handles.MainAxes);
imagesc(label2rgb(Canvas.Checkers(400, [], 50), [1 1 1] * 206/255, [1 1 1]));
set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');
drawnow;

if Q.isfield(CheeseSquareData.obj.Meta, 'Scale')
    CheeseScaleData.Scale = CheeseSquareData.obj.Meta.Scale;
else
    CheeseScaleData.Scale = [];
end

CheeseScaleData.rect = [];

if isfield(CheeseSquareData, 'obj') 
    if ~isa(CheeseSquareData.obj, 'Video')
        CheeseSquareData.obj = CheeseSquare(CheeseSquareData.obj);
    end
    ReadMetaData
else
    LoadFile_Callback
end

try
    %axes(CheeseScaleData.Handles.MainAxes);
    %GetFrame(1);
    EnhanceImage_Callback;
catch
end

function ProgressBar(str, p)
global CheeseScaleData;
mb = CheeseScaleData.Handles.MessageLine;
set(mb, 'string', [' ' str]);

% if nargin < 2
%     p = 1;
% end
% axes(CheeseMarkRegionsData.Handles.ProgressBar);
% z1 = ones(1, round(50 * p));
% z2 = ones(1, round(50 * (1-p)));
% c1 = Colors.PrettyBlue;
% c2 = Colors.PrettyRed;
% imagesc([cat(3, z1*c1(1), z1*c1(2), z1*c1(3)) cat(3, z2*c2(1), z2*c2(2), z2*c2(3))]);
% axis off;
axes(CheeseScaleData.Handles.MainAxes);

% --- Outputs from this function are returned to the command line.
function varargout = CheeseScale_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function FileName_Callback(hObject, eventdata, handles)
% hObject    handle to FileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FileName as text
%        str2double(get(hObject,'String')) returns contents of FileName as a double


% --- Executes during object creation, after setting all properties.
function FileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ChooseFile.
function ChooseFile_Callback(hObject, eventdata, handles)
% hObject    handle to ChooseFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseScaleData;

[filename,pathname] = uigetfile('*.avi;*.obj.mat');
if filename ~= 0
    filename = fullfile(pathname, filename);
    set(CheeseScaleData.Handles.FileName, 'String', filename);
end
LoadFile_Callback;

function ReadMetaData
global CheeseScaleData;
global CheeseSquareData;
set(CheeseScaleData.Handles.VideoPosition, 'Min', 1, 'Max', CheeseSquareData.obj.Video.NumberOfFrames, 'value', 1, 'sliderstep', [1 1]/(CheeseSquareData.obj.Video.NumberOfFrames - 1));



% --- Executes on button press in LoadFile.
function LoadFile_Callback(hObject, eventdata, handles)
% hObject    handle to LoadFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseSquareData;
global CheeseScaleData;

c = get(CheeseScaleData.Handles.InfoPanel, 'children');
for i=1:length(c)
    delete(c(i));
end

filename = get(CheeseScaleData.Handles.FileName, 'String');
if isempty(filename)
    return;
end
try
    
    CheeseSquareData.obj = CheeseSquare(filename);
    
    ReadMetaData;
    drawnow;

    GetFrame(1);
    drawnow;
    
catch me
    warndlg({'unable to load file'; me.message});
end


function GetFrame(num)
global CheeseSquareData;
global CheeseScaleData;
try
    ProgressBar('Loading frame...');
    CheeseSquareData.obj.Video.FrameNumber = num;
    CheeseScaleData.Image = CheeseSquareData.obj.Video.CurrentFrame;
    imagesc(CheeseScaleData.Image);
    set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');
    set(CheeseScaleData.Handles.FrameNumTB, 'string', sprintf('%d', CheeseSquareData.obj.Video.FrameNumber));
catch
end
ProgressBar('');
AddMarks(false);

% --- Executes on slider movement.
function VideoPosition_Callback(hObject, eventdata, handles)
% hObject    handle to VideoPosition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global CheeseSquareData;
framenum = round(get(hObject,'Value'));
GetFrame(framenum);
ProgressBar(sec2time(framenum / CheeseSquareData.obj.Video.FrameRate));
drawnow;

% --- Executes during object creation, after setting all properties.
function VideoPosition_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VideoPosition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in PlayStop.
function PlayStop_Callback(hObject, eventdata, handles)
% hObject    handle to PlayStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseScaleData;
global CheeseSquareData;

set(CheeseScaleData.Handles.PlayStop, 'String', '[]');
CheeseScaleData.IsPlaying = ~CheeseScaleData.IsPlaying;
tic
while true
    if ~CheeseScaleData.IsPlaying
        break;
    end
    t = toc;
    if t >= 1/CheeseSquareData.obj.Video.FrameRate
        framenum = round(get(CheeseScaleData.Handles.VideoPosition,'Value') + 1);
        CheeseSquareData.obj.Video.FrameNumber = framenum;
        imagesc(CheeseSquareData.obj.Video.CurrentFrame);
        set(CheeseScaleData.Handles.VideoPosition, 'Value', framenum);
        set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');
        ProgressBar(sec2time(framenum / CheeseSquareData.obj.Video.FrameRate));
        drawnow;
    end
end
set(CheeseScaleData.Handles.PlayStop, 'String', '>');


% --- Executes on button press in JumpLeft.
function JumpLeft_Callback(hObject, eventdata, handles)
% hObject    handle to JumpLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseScaleData;
global CheeseSquareData;
framenum = max(round(get(CheeseScaleData.Handles.VideoPosition,'Value') - 60 * CheeseSquareData.obj.Video.FrameRate, 0));
set(CheeseScaleData.Handles.VideoPosition, 'Value', framenum);
GetFrame(framenum);
ProgressBar(sec2time(framenum / CheeseSquareData.obj.Video.FrameRate));
drawnow;

% --- Executes on button press in JumpRight.
function JumpRight_Callback(hObject, eventdata, handles)
% hObject    handle to JumpRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseScaleData;
global CheeseSquareData;
framenum = min(round(get(CheeseScaleData.Handles.VideoPosition,'Value') + 60 * CheeseSquareData.obj.Video.FrameRate), CheeseSquareData.obj.Video.NumberOfFrames);
set(CheeseScaleData.Handles.VideoPosition, 'Value', framenum);
GetFrame(framenum);
ProgressBar(sec2time(framenum / CheeseSquareData.obj.Video.FrameRate));
drawnow;


% --- Executes on button press in NextButton.
function NextButton_Callback(hObject, eventdata, handles)
% hObject    handle to NextButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseScaleData;
global CheeseSquareData;
CheeseSquareData.obj.Meta.Scale = CheeseScaleData.Scale;
close(gcf);
CheeseBackground

function AddMarks(isnew)
global CheeseSquareData;
global CheeseScaleData;
try
    if ~isempty(CheeseScaleData.rect)
        delete(CheeseScaleData.rect);
        CheeseScaleData.rect = [];
    end
catch
end
if ~isempty(CheeseScaleData.Scale)
    set(CheeseScaleData.Handles.WidthEdit, 'String', num2str(CheeseScaleData.Scale.ArenaWidth));
    set(CheeseScaleData.Handles.HeightEdit, 'String', num2str(CheeseScaleData.Scale.ArenaHeight));
    CheeseScaleData.rect = imrect(gca, CheeseScaleData.Scale.ArenaCoord);
    CheeseScaleData.rect.addNewPositionCallback(@NewPositionCallback);
else
    if isnew
        width = str2double(get(CheeseScaleData.Handles.WidthEdit, 'String'));
        height = str2double(get(CheeseScaleData.Handles.HeightEdit, 'String'));
        w = CheeseSquareData.obj.Video.Width*6/8;
        h = w / width * height;
        wpx = CheeseSquareData.obj.Video.Width;
        hpx = CheeseSquareData.obj.Video.Height;
        %CheeseScaleData.rect.setPosition();
        CheeseScaleData.Scale.ArenaWidth = width;
        CheeseScaleData.Scale.ArenaHeight = height;
        CheeseScaleData.Scale.ArenaCoord = [(wpx - w)/2, (hpx - h)/2, w, h];
        CheeseScaleData.rect = imrect(gca, CheeseScaleData.Scale.ArenaCoord);
        CheeseScaleData.rect.addNewPositionCallback(@NewPositionCallback);
    end
end

function NewPositionCallback(hObject, eventdata) 
global CheeseSquareData;
global CheeseScaleData;
CheeseScaleData.Scale.ArenaCoord = CheeseScaleData.rect.getPosition;

% --- Executes on button press in MarkButton.
function MarkButton_Callback(hObject, eventdata, handles)
% hObject    handle to MarkButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AddMarks(true);

% --- Executes on button press in SegmentButton.
function SegmentButton_Callback(hObject, eventdata, handles)
% hObject    handle to SegmentButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseScaleData;
global CheeseSquareData;
%%
if Q.isfield(CheeseSquareData.obj, 'Background.im')
[map, img] = ImageSegmentClassic(CheeseSquareData.obj, CheeseSquareData.obj.Video.CurrentFrame, struct('Output', false));
r = rgb2gray(im2double(img));
g = r; b = r;
g(map) = 0;
b(map) = 0;
perim = bwperim(map);
perim = convn(double(perim), ones(3), 'same') > 0;
r(perim) = 1;
g(perim) = 0;
b(perim) = 0;
imagesc(cat(3, r, g, b));
end



function WidthEdit_Callback(hObject, eventdata, handles)
% hObject    handle to WidthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WidthEdit as text
%        str2double(get(hObject,'String')) returns contents of WidthEdit as a double
global CheeseScaleData;
%if ~isempty(CheeseScaleData.Scale)
    CheeseScaleData.Scale.ArenaWidth = str2double(get(CheeseScaleData.Handles.WidthEdit, 'String'));
%end


% --- Executes during object creation, after setting all properties.
function WidthEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WidthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function HeightEdit_Callback(hObject, eventdata, handles)
% hObject    handle to HeightEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HeightEdit as text
%        str2double(get(hObject,'String')) returns contents of HeightEdit as a double
global CheeseScaleData;
%if ~isempty(CheeseScaleData.Scale)
    CheeseScaleData.Scale.ArenaHeight = str2double(get(CheeseScaleData.Handles.HeightEdit, 'String'));
%end


% --- Executes during object creation, after setting all properties.
function HeightEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HeightEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit5_CreateFcn(hObject, eventdata, handles)

function edit6_CreateFcn(hObject, eventdata, handles)



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ClearButton.
function ClearButton_Callback(hObject, eventdata, handles)
% hObject    handle to ClearButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseScaleData;
try
    if ~isempty(CheeseScaleData.rect)
        delete(CheeseScaleData.rect);
        CheeseScaleData.rect = [];
    end
catch
end

CheeseScaleData.Scale = [];


% --- Executes on button press in PrevButton.
function PrevButton_Callback(hObject, eventdata, handles)
% hObject    handle to PrevButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf);
CheeseInit


% --- Executes on button press in EnhanceImage.
function EnhanceImage_Callback(hObject, eventdata, handles)
% hObject    handle to EnhanceImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EnhanceImage
global CheeseScaleData;
global CheeseSquareData;
axes(CheeseScaleData.Handles.MainAxes);
if Q.isfield(CheeseSquareData.obj.Background, 'im') && ~isempty(CheeseSquareData.obj.Background.im)
    img = CheeseSquareData.obj.Background.im;
else
    try
        img = CheeseSquareData.obj.Video.CurrentFrame;
    catch
        img = label2rgb(Canvas.Checkers(400, [], 50), [1 1 1] * 206/255, [1 1 1]);
    end
end
if get(CheeseScaleData.Handles.EnhanceImage, 'value')
    imagesc(CV.HistEq(img));
else
    imagesc(img);
end
set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');
drawnow;
AddMarks(false);

function FrameNumTB_Callback(hObject, eventdata, handles)
% hObject    handle to FrameNumTB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrameNumTB as text
%        str2double(get(hObject,'String')) returns contents of FrameNumTB as a double
global CheeseScaleData;
global CheeseSquareData;
%set(CheeseInitData.Handles.VideoPosition, 'value', str2double(get(CheeseInitData.Handles.FrameNumTB, 'string')));
GetFrame(str2double(get(CheeseScaleData.Handles.FrameNumTB, 'string')));


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
