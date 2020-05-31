function varargout = CheeseInit(varargin)
% CHEESEINIT MATLAB code for CheeseInit.fig
%      CHEESEINIT, by itself, creates a new CHEESEINIT or raises the existing
%      singleton*.
%
%      H = CHEESEINIT returns the handle to a new CHEESEINIT or the handle to
%      the existing singleton*.
%
%      CHEESEINIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHEESEINIT.M with the given input arguments.
%
%      CHEESEINIT('Property','Value',...) creates a new CHEESEINIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CheeseInit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CheeseInit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CheeseInit

% Last Modified by GUIDE v2.5 31-Dec-2014 11:53:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CheeseInit_OpeningFcn, ...
                   'gui_OutputFcn',  @CheeseInit_OutputFcn, ...
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


% --- Executes just before CheeseInit is made visible.
function CheeseInit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CheeseInit (see VARARGIN)

% Choose default command line output for CheeseInit
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CheeseInit wait for user response (see UIRESUME)
% uiwait(handles.figure1);

warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jframe=get(handles.figure1, 'javaframe');
jIcon=javax.swing.ImageIcon('.\mouse-square-small.gif');
jframe.setFigureIcon(jIcon);

global CheeseInitData;
global CheeseSquareData;

CheeseInitData.Handles = handles;
CheeseInitData.IsPlaying = false;
set(CheeseInitData.Handles.PlayStop, 'String', '>');

axes(CheeseInitData.Handles.MainAxes);
imagesc(label2rgb(Canvas.Checkers(400, [], 50), [1 1 1] * 206/255, [1 1 1]));
set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');
drawnow;

if isfield(CheeseSquareData, 'obj') 
    if ~isa(CheeseSquareData.obj, 'Video')
        CheeseSquareData.obj = CheeseSquare(CheeseSquareData.obj);
    end
    set(CheeseInitData.Handles.FileName, 'String', CheeseSquareData.obj.Video.Filename);
    ReadMetaData
else
    LoadFile_Callback
end

try
    axes(CheeseInitData.Handles.MainAxes);
    CheeseSquareData.obj.Video.FrameNumber = 1;
    imagesc(CheeseSquareData.obj.Video.CurrentFrame);
    set(CheeseInitData.Handles.FrameNumTB, 'string', '1');
    set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');
    drawnow;
catch
end

function ProgressBar(str, p)
global CheeseInitData;
mb = CheeseInitData.Handles.MessageLine;
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
axes(CheeseInitData.Handles.MainAxes);

function handles = AddToInfo(handles, fig, form, name, label, value)
style = {'FontName', 'Segoe UI', 'FontSize', 9.0, 'BackgroundColor', 'w'};
idx = length(get(handles.InfoPanel, 'children'));
pos = get(handles.InfoPanel, 'position');
handles.([name 'label']) = uicontrol(handles.InfoPanel, ...
    'Style','text',...
    'String', label,...
    'Position', [10 pos(4) - (45 + idx * 15), pos(3)/2-20, 17], ...
    'HorizontalAlignment', 'left', ...
    style{:});

handles.(name) = uicontrol(handles.InfoPanel, ...
    'Tag', name, ...
    'Style','edit',...
    'String', num2str(value),...
    'Position', [pos(3)/2 + 2 pos(4) - (48 + idx * 15), pos(3)/2-17, 21], ...
    style{:});

guidata(handles.figure1, handles);


% --- Outputs from this function are returned to the command line.
function varargout = CheeseInit_OutputFcn(hObject, eventdata, handles) 
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
global CheeseInitData;

[filename,pathname] = uigetfile('*.avi;*.obj.mat');
if filename ~= 0
    filename = fullfile(pathname, filename);
    set(CheeseInitData.Handles.FileName, 'String', filename);
end
LoadFile_Callback;

function ReadMetaData
global CheeseInitData;
global CheeseSquareData;
CheeseInitData.Handles = ...
    AddToInfo(CheeseInitData.Handles, CheeseInitData.Handles.figure1, CheeseInitData.Handles.InfoPanel, 'Duration', 'Total duration', DateTime.SecToString(CheeseSquareData.obj.Video.NumberOfFrames / CheeseSquareData.obj.Video.FrameRate));
CheeseInitData.Handles = ...
    AddToInfo(CheeseInitData.Handles, CheeseInitData.Handles.figure1, CheeseInitData.Handles.InfoPanel, 'nFrames', 'Number of frames', CheeseSquareData.obj.Video.NumberOfFrames);
try
    CheeseInitData.Handles = AddToInfo(CheeseInitData.Handles, CheeseInitData.Handles.figure1, CheeseInitData.Handles.InfoPanel, ...
        'CreationTime', 'Created in ', datestr(CheeseSquareData.obj.Meta.CreationDate, 'dd/mm/yyyy HH:MM:SS'));
catch
end

try
    CheeseInitData.Handles = AddToInfo(CheeseInitData.Handles, CheeseInitData.Handles.figure1, CheeseInitData.Handles.InfoPanel, ...
        'Remarks', 'Remarks', CheeseSquareData.obj.Meta.Remarks);
catch
end
set(CheeseInitData.Handles.VideoPosition, 'Min', 1, 'Max', CheeseSquareData.obj.Video.NumberOfFrames, 'value', 1, 'sliderstep', [1 1]/(CheeseSquareData.obj.Video.NumberOfFrames - 1));



% --- Executes on button press in LoadFile.
function LoadFile_Callback(hObject, eventdata, handles)
% hObject    handle to LoadFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseSquareData;
global CheeseInitData;

c = get(CheeseInitData.Handles.InfoPanel, 'children');
for i=1:length(c)
    delete(c(i));
end

filename = get(CheeseInitData.Handles.FileName, 'String');
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
global CheeseInitData;
try
    ProgressBar('Loading frame...');
    CheeseSquareData.obj.Video.FrameNumber = num;
    set(CheeseInitData.Handles.FrameNumTB, 'string', num2str(num));
    CheeseInitData.Image = CheeseSquareData.obj.Video.CurrentFrame;
    imagesc(CheeseInitData.Image);
    set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');
catch
end
ProgressBar('');

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
ProgressBar(DateTime.SecToString(framenum / CheeseSquareData.obj.Video.FrameRate));
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
global CheeseInitData;
global CheeseSquareData;

set(CheeseInitData.Handles.PlayStop, 'String', '[]');
CheeseInitData.IsPlaying = ~CheeseInitData.IsPlaying;
tic
while true
    if ~CheeseInitData.IsPlaying
        break;
    end
    t = toc;
    if t >= 1/CheeseSquareData.obj.Video.FrameRate
        framenum = round(get(CheeseInitData.Handles.VideoPosition,'Value') + 1);
        CheeseSquareData.obj.Video.FrameNumber = framenum;
        imagesc(CheeseSquareData.obj.Video.CurrentFrame);
        set(CheeseInitData.Handles.VideoPosition, 'Value', framenum);
        set(CheeseInitData.Handles.FrameNumTB, 'string', num2str(num));
        set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');
        ProgressBar(DateTime.SecToString(framenum / CheeseSquareData.obj.Video.FrameRate));
        drawnow;
    end
end
set(CheeseInitData.Handles.PlayStop, 'String', '>');


% --- Executes on button press in JumpLeft.
function JumpLeft_Callback(hObject, eventdata, handles)
% hObject    handle to JumpLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseInitData;
global CheeseSquareData;
framenum = max(round(get(CheeseInitData.Handles.VideoPosition,'Value') - 60 * CheeseSquareData.obj.Video.FrameRate, 0));
set(CheeseInitData.Handles.VideoPosition, 'Value', framenum);
GetFrame(framenum);
ProgressBar(DateTime.SecToString(framenum / CheeseSquareData.obj.Video.FrameRate));
drawnow;

% --- Executes on button press in JumpRight.
function JumpRight_Callback(hObject, eventdata, handles)
% hObject    handle to JumpRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseInitData;
global CheeseSquareData;
framenum = min(round(get(CheeseInitData.Handles.VideoPosition,'Value') + 60 * CheeseSquareData.obj.Video.FrameRate), CheeseSquareData.obj.Video.NumberOfFrames);
set(CheeseInitData.Handles.VideoPosition, 'Value', framenum);
GetFrame(framenum);
ProgressBar(DateTime.SecToString(framenum / CheeseSquareData.obj.Video.FrameRate));
drawnow;


% --- Executes on button press in NextButton.
function NextButton_Callback(hObject, eventdata, handles)
% hObject    handle to NextButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf);
CheeseScale

% --- Executes on button press in ExternalView.
function ExternalView_Callback(hObject, eventdata, handles)
% hObject    handle to ExternalView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseSquareData;
system(CheeseSquareData.obj.Video.Filename);


% --- Executes on button press in SegmentButton.
function SegmentButton_Callback(hObject, eventdata, handles)
% hObject    handle to SegmentButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseInitData;
global CheeseSquareData;
%%
if Q.isfield(CheeseSquareData.obj, 'Background.im')
    classicseg = true;
    if Q.isfield(CheeseSquareData.obj, 'Colors.Histrogram.H')
        try
            im = CheeseSquareData.obj.Video.CurrentFrame;
            [~, ~, seg] = CheeseColorSegment(CheeseSquareData.obj.ToClassical, im);
            seg = imresize(seg, [size(im, 1), size(im, 2)]);
            for i=1:size(im, 3)
                ch = im(:, :, i);
                sch = seg(:,:,i);
                ch(any(seg, 3)) = sch(any(seg, 3));
                im(:, :, i) = ch;
            end
            axes(CheeseInitData.Handles.MainAxes);
            imagesc(im);
            set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');
            drawnow;
            classicseg = false;
        catch
            classicseg = true;
        end
    end
    if classicseg
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
            set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');
            drawnow;
    end
    
end



function FrameNumTB_Callback(hObject, eventdata, handles)
% hObject    handle to FrameNumTB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrameNumTB as text
%        str2double(get(hObject,'String')) returns contents of FrameNumTB as a double
global CheeseInitData;
global CheeseSquareData;
%set(CheeseInitData.Handles.VideoPosition, 'value', str2double(get(CheeseInitData.Handles.FrameNumTB, 'string')));
GetFrame(str2double(get(CheeseInitData.Handles.FrameNumTB, 'string')));


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
