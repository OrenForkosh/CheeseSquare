function varargout = CheeseSave(varargin)
% CHEESESAVE MATLAB code for CheeseSave.fig
%      CHEESESAVE, by itself, creates a new CHEESESAVE or raises the existing
%      singleton*.
%
%      H = CHEESESAVE returns the handle to a new CHEESESAVE or the handle to
%      the existing singleton*.
%
%      CHEESESAVE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHEESESAVE.M with the given input arguments.
%
%      CHEESESAVE('Property','Value',...) creates a new CHEESESAVE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CheeseSave_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CheeseSave_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CheeseSave

% Last Modified by GUIDE v2.5 23-Sep-2014 16:21:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CheeseSave_OpeningFcn, ...
                   'gui_OutputFcn',  @CheeseSave_OutputFcn, ...
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


% --- Executes just before CheeseSave is made visible.
function CheeseSave_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CheeseSave (see VARARGIN)

% Choose default command line output for CheeseSave
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CheeseSave wait for user response (see UIRESUME)
% uiwait(handles.Finalize);

warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jframe=get(handles.Finalize, 'javaframe');
jIcon=javax.swing.ImageIcon('.\mouse-square-small.gif');
jframe.setFigureIcon(jIcon);

global CheeseSaveData;
global CheeseSquareData;

CheeseSaveData.Handles = handles;

axes(CheeseSaveData.Handles.BkgAxes);
imagesc(label2rgb(Canvas.Checkers(400, [], 50), [1 1 1] * 206/255, [1 1 1]));
set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');
title('Background model', 'FontSize', 12, 'FontName', 'Calibri Light');
Fig.Fix
box on

axes(CheeseSaveData.Handles.RegionsAxes);
imagesc(label2rgb(Canvas.Checkers(400, [], 50), [1 1 1] * 206/255, [1 1 1]));
set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');
title('Regions of interest', 'FontSize', 12, 'FontName', 'Calibri Light');
Fig.Fix
box on
drawnow;

axes(CheeseSaveData.Handles.BkgAxes);
try
    if all(isfield(CheeseSquareData.obj.Background, {'im', 'value', 'mean', 'std'}))
        imagesc(CheeseSquareData.obj.Background.im);
        set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');
        title('Background model', 'FontSize', 12, 'FontName', 'Calibri Light');
        Fig.Fix
        box on
    end
catch
end

axes(CheeseSaveData.Handles.RegionsAxes);
try
    imagesc(CV.HistEq(CheeseSquareData.obj.Background.im));
    set(gca, 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');
    title('Regions of interest', 'FontSize', 12, 'FontName', 'Calibri Light');
    Fig.Fix
    box on
    cmap = repmat(Colormaps.Categorical, [3, 1]);
    for r=find(CheeseSquareData.obj.ROI.IsAvail)
        %%
        b = bwboundaries(CheeseSquareData.obj.ROI.Regions{r});
        p = bwtraceboundary(CheeseSquareData.obj.ROI.Regions{r}, b{1}(1, :), 'N');
        patch(p(:, 2), p(:, 1), 'y', 'EdgeColor', cmap(r, :), 'FaceColor', cmap(r, :), 'FaceAlpha', .2)
        %%
        pt = CheeseSquareData.obj.ROI.Points{r};
        inside = Q.inrange(pt(:, 1), 0, CheeseSquareData.obj.Video.Width) & Q.inrange(pt(:, 2), 0, CheeseSquareData.obj.Video.Height);
        %%
        [~, idx] = max(sum(pt, 2) .* inside);
        pos = CheeseSquareData.obj.ROI.Points{r};
        str = sprintf('%s', CheeseSquareData.obj.ROI.RegionNames{r});
        text(pos(idx, 1)+10, pos(idx, 2)+10, str, 'BackgroundColor', cmap(r, :), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'color', 'w', 'FontSize', 6);
    end
catch
end

%%
try
    [pathstr, name, ext] = fileparts(CheeseSquareData.obj.Video.Filename);
    list = dir([pathstr filesep '*' ext]);
    match = ~cellfun(@isempty, regexp({list.name}, ['^' regexprep(name, 'day\d+', 'day\\d+')]));
    exact = ~cellfun(@isempty, regexp({list.name}, [name ext]));
    names = {list.name};
    matchingnames = names(match & ~exact);
    set(CheeseSaveData.Handles.RepeatList, 'String', matchingnames);
    set(CheeseSaveData.Handles.RepeatList, 'value', 1:length(matchingnames));
catch
end


% --- Outputs from this function are returned to the command line.
function varargout = CheeseSave_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in RepeatList.
function RepeatList_Callback(hObject, eventdata, handles)
% hObject    handle to RepeatList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns RepeatList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from RepeatList


% --- Executes during object creation, after setting all properties.
function RepeatList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RepeatList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NoneButton.
function NoneButton_Callback(hObject, eventdata, handles)
% hObject    handle to NoneButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseSaveData;
set(CheeseSaveData.Handles.RepeatList, 'value', []);

% --- Executes on button press in AllButton.
function AllButton_Callback(hObject, eventdata, handles)
% hObject    handle to AllButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CheeseSaveData;
matchingnames = get(CheeseSaveData.Handles.RepeatList, 'String');
set(CheeseSaveData.Handles.RepeatList, 'value', 1:length(matchingnames));

function Message(str)
%%
global CheeseSaveData;
mb = CheeseSaveData.Handles.MessageBox;
prev = get(mb, 'string');
if isempty(prev)
    set(mb, 'string', {[' - ' str]});
else
    set(mb, 'string', {prev{:}, [' - ' str]});
end
drawnow



% --- Executes on button press in TrackNow.
function TrackNow_Callback(hObject, eventdata, handles)
% hObject    handle to TrackNow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in FinalizeButton.
function FinalizeButton_Callback(hObject, eventdata, handles)
% hObject    handle to FinalizeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%
global CheeseSaveData;
global CheeseSquareData;
CheeseSquareData.obj.Save;

matchingnames = get(CheeseSaveData.Handles.RepeatList, 'String');
idx = get(CheeseSaveData.Handles.RepeatList, 'value');
pathstr = fileparts(CheeseSquareData.obj.Video.Filename);
orig = CheeseSquareData.obj;
for i=idx
    curr = fullfile(pathstr, matchingnames{i});
    Message(['- Processing ' curr]);
    obj = CheeseSquare(curr);
    obj.ROI = orig.ROI;
    obj.nSubjects = orig.nSubjects;
    obj.Colors = orig.Colors;
    %obj.Colors.Marks = struct('x', [], 'y', [], 'id', [], 'frame', [], 'color', []);
    obj.Background = orig.Background;
    obj.Meta = Q.cpfield(orig.Meta, obj.Meta, 'Scale', true);
    obj.Colors.Marks = Q.setfield(obj.Colors.Marks, 'x', []);
    obj.Colors.Marks = Q.setfield(obj.Colors.Marks, 'y', []);
    obj.Colors.Marks = Q.setfield(obj.Colors.Marks, 'frame', []);
    obj.Save;
end
close(CheeseSaveData.Handles.Finalize);

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


% --- Executes on button press in PrevButton.
function PrevButton_Callback(hObject, eventdata, handles)
% hObject    handle to PrevButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf);
CheeseMarkRegions
