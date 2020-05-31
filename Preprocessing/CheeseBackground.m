function varargout = CheeseBackground(varargin)
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

% Last Modified by GUIDE v2.5 23-Dec-2014 12:44:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @CheeseBackground_OpeningFcn, ...
    'gui_OutputFcn',  @CheeseBackground_OutputFcn, ...
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
function CheeseBackground_OpeningFcn(hObject, eventdata, handles, varargin)
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

global CheeseBackgroundData;
global CheeseSquareData;

opt.nFrames = 40;
opt.Quantile = .5;
opt.StartTime = [];
%opt.EndTime = CheeseSquareData.obj.Video.NumberOfFrames / CheeseSquareData.obj.Video.FrameRate;
opt.EndTime = [];
opt.DarkThresh = .2;

handles = AddToForm(handles, handles.figure1, handles.Form, 'nFrames', 'Number of frames', opt.nFrames);
handles = AddToForm(handles, handles.figure1, handles.Form, 'Quantile', 'Quantile', opt.Quantile);
handles = AddToForm(handles, handles.figure1, handles.Form, 'StartTime', 'Start time [sec]', opt.StartTime);
handles = AddToForm(handles, handles.figure1, handles.Form, 'EndTime', 'End time [sec]', opt.EndTime);
handles = AddToForm(handles, handles.figure1, handles.Form, 'DarkThresh', 'Darkness threshold', opt.DarkThresh);

CheeseBackgroundData.Options = opt;
CheeseBackgroundData.Handles = handles;

axes(CheeseBackgroundData.Handles.ProgressBar);
imagesc(reshape(Colors.PrettyBlue, [1 1 3]));
%axis off;
set(CheeseBackgroundData.Handles.ProgressBar, 'XTick', [], 'YTick', [], 'XColor', Colors.PrettyBlue, 'YColor', Colors.PrettyBlue);
set(CheeseBackgroundData.Handles.VideoPosition, 'Min', 1, 'Max', CheeseSquareData.obj.Video.NumberOfFrames, 'value', 1, 'sliderstep', [1 1]/(CheeseSquareData.obj.Video.NumberOfFrames - 1));

axes(CheeseBackgroundData.Handles.MainAxes);
%imagesc(label2rgb(1+Canvas.AddBorder(1*Canvas.Checkers(400, [], 50), 1, 2), [1 1 1; .8 .8 .8; 0 0 0]));
imagesc(label2rgb(Canvas.Checkers(400, [], 50), [1 1 1] * 206/255, [1 1 1]));
set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');
drawnow

warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jframe=get(handles.figure1, 'javaframe');
jIcon=javax.swing.ImageIcon('.\mouse-square-small.gif');
jframe.setFigureIcon(jIcon);

global CheeseSquareData;
if all(isfield(CheeseSquareData.obj.Background, {'im', 'value', 'mean', 'std'}))
    imagesc(CheeseSquareData.obj.Background.im);
    set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');
end
ProgressBar('Background model');

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
global CheeseBackgroundData;
mb = CheeseBackgroundData.Handles.MessageBox;
prev = get(mb, 'string');
if isempty(prev)
    set(mb, 'string', {[' - ' str]});
else
    set(mb, 'string', {prev{:}, [' - ' str]});
end
drawnow

function ProgressBar(str, p)
global CheeseBackgroundData;
mb = CheeseBackgroundData.Handles.MessageLine;
set(mb, 'string', [' ' str]);

if nargin < 2
    p = 1;
end
axes(CheeseBackgroundData.Handles.ProgressBar);
z1 = ones(1, round(50 * p));
z2 = ones(1, round(50 * (1-p)));
c1 = Colors.PrettyBlue;
c2 = Colors.PrettyRed;
imagesc([cat(3, z1*c1(1), z1*c1(2), z1*c1(3)) cat(3, z2*c2(1), z2*c2(2), z2*c2(3))]);
axis off;
axes(CheeseBackgroundData.Handles.MainAxes);


function Option_Callback(hObject, eventdata, handles)
global CheeseBackgroundData;
tag = get(hObject, 'tag');
CheeseBackgroundData.Options.(tag) = str2num(get(hObject, 'string'));

function Start()
global CheeseBackgroundData;
set(CheeseBackgroundData.Handles.RunButton, 'String', 'Training...');
mb = CheeseBackgroundData.Handles.MessageBox;
set(mb, 'string', '');

% CheeseBackgroundData.GUI.nDots = 0;
% CheeseBackgroundData.GUI.Timer = timer(...
%     'TimerFcn', @IRun,... 
%     'ExecutionMode', 'fixedSpacing', ...
%     'Period', .5);
% start(CheeseBackgroundData.GUI.Timer)

function IRun(mTimer, ~)
global CheeseBackgroundData;
set(CheeseBackgroundData.Handles.RunButton, 'String', ['Training...']);
%dbstack

function Finish()
global CheeseBackgroundData;
set(CheeseBackgroundData.Handles.RunButton, 'String', 'Train');

function Success()
global CheeseBackgroundData;
set(CheeseBackgroundData.Handles.NextButton, 'Enable', 'on');

% --- Outputs from this function are returned to the command line.
function varargout = CheeseBackground_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in NextButton.
function NextButton_Callback(hObject, eventdata, handles)
% hObject    handle to NextButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseSquareData;
if ~all(isfield(CheeseSquareData.obj.Background, {'im', 'value', 'mean', 'std'}))
    res = questdlg('Background model was not trained', 'Do you want to continue?', 'Train', 'Skip', 'Cancel', 'Train');
    if strcmpi(res, 'Train')
        RunButton_Callback();
        return;
    end
    if strcmpi(res, 'Cancel')
        return;
    end
end
close(gcf);
CheeseMarkColors

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

% --- Executes on button press in RunButton.
function RunButton_Callback(hObject, eventdata, handles)
% hObject    handle to RunButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseSquareData;
global CheeseBackgroundData;
Start();

succ = false;
try
    obj = CheeseSquareData.obj;
    vid = obj.Video;
    opt = CheeseBackgroundData.Options;
    opt.Output = true;
    buffer = uint8(zeros(vid.Height, vid.Width, 3, opt.nFrames));
    %%
    Message('finding dark period');
    ProgressBar('finding dark period', 0);
    obj = CheeseFindDarkPeriod(obj, opt.DarkThresh);

    %%
    nchars = 0;
    i = 1;
    Message('reading frames');
    while i <= opt.nFrames
        if isempty(opt.StartTime)
            framenum = round(obj.Meta.ExperimentTimes(1)*vid.FrameRate) + randi(round((obj.Meta.ExperimentTimes(2) - obj.Meta.ExperimentTimes(1))*vid.FrameRate));
        else
            framenum = round(opt.StartTime*vid.FrameRate) + randi(round((opt.EndTime - opt.StartTime)*vid.FrameRate));
        end
        ProgressBar(sprintf('reading frame no. %d', framenum), .5*i/opt.nFrames);
        vid.FrameNumber = framenum;
        img = im2uint8(vid.CurrentFrame);
        buffer(:, :, :, i) = img;
        %%
        if opt.Output
            imagesc(img);
            axis off;
            drawnow;
        end
        %%
        i = i + 1;
    end
    
    %%
    Message('computing median image');
    ProgressBar('computing median image', .5+0/10);
    obj.Background.im = median(buffer, 4);
    if opt.Output
        imagesc(obj.Background.im);
        axis off;
    end
    
    %%
    vbuffer = squeeze(max(buffer, [], 3));
    
    %%
    Message('computing median brightness (value)');
    ProgressBar('computing median brightness (value)', .5+2/10);
    obj.Background.value = median(vbuffer, 3);
    if opt.Output
        imagesc(obj.Background.value);
        axis off;
    end
    
    %%
    Message('matching histograms');
    ProgressBar('matching histograms', .5+3/10);
    for i=1:opt.nFrames
        vbuffer(:, :, i) = imhistmatch(vbuffer(:, :, i), obj.Background.value);
    end
    
    %%
    Message('computing background model');
    ProgressBar('computing background model', .5+4/10);
    obj.Background.mean = im2double(median(vbuffer, 3));
    if opt.Output
        imagesc(obj.Background.mean);
        axis off;
    end
    
    %%
    Message('computing background model variablitity');
    ProgressBar('computing background model variablitity', .5+5/10);
    D = [];
    for i=1:opt.nFrames
        d = im2double(vbuffer(:, :, i)) - obj.Background.mean;
        D = [D; d(:)]; %#ok<AGROW>
    end
    obj.Background.std = quantile(abs(D), opt.Quantile) / norminv(opt.Quantile / 2 + .5);
    %%
    
    %%
    if opt.Output
        imagesc(obj.Background.im);
        axis off;
    end
    %%    
    Message('done');
    ProgressBar('', 1);
    succ = true;
catch err
    %%
     msgbox({'Unable to train a background model:',  ['   ' err.message]}, 'Failed','error');
end
CheeseSquareData.obj = obj;
Finish();
if succ
    Success();
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
        if all(isfield(src.obj.Background, {'im', 'value', 'mean', 'std'}))
            CheeseSquareData.obj.Background = src.obj.Background;
        else
            error(MException('failed'));
        end
    catch
        msgbox('unable to load backgroud model from file');
    end
end


% --- Executes on slider movement.
function VideoPosition_Callback(hObject, eventdata, handles)
% hObject    handle to VideoPosition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global CheeseSquareData;
global CheeseBackgroundData;

framenum = round(get(hObject,'Value'));
ProgressBar(sprintf('reading frame - time %.2f sec', framenum / CheeseSquareData.obj.Video.FrameRate));
CheeseSquareData.obj.Video.FrameNumber = framenum;
img = im2uint8(CheeseSquareData.obj.Video.CurrentFrame);
imagesc(img);
set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');

set(CheeseBackgroundData.Handles.ProgressBar, 'XTick', [], 'YTick', [], 'XColor', Colors.PrettyBlue, 'YColor', Colors.PrettyBlue);


% --- Executes on button press in PreviousButton.
function PreviousButton_Callback(hObject, eventdata, handles)
% hObject    handle to PreviousButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf);
CheeseScale



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
disp 'a'

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
