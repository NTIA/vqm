function varargout = bvqm_pc_uncertss(varargin)
% BVQM_PC_BVQM_PC_UNCERTSS M-file for bvqm_pc_bvqm_pc_uncertss.fig
%      BVQM_PC_BVQM_PC_UNCERTSS, by itself, creates a new BVQM_PC_BVQM_PC_UNCERTSS or raises the existing
%      singleton*.
%
%      H = BVQM_PC_BVQM_PC_UNCERTSS returns the handle to a new BVQM_PC_BVQM_PC_UNCERTSS or the handle to
%      the existing singleton*.
%
%      BVQM_PC_BVQM_PC_UNCERTSS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BVQM_PC_BVQM_PC_UNCERTSS.M with the given input arguments.
%
%      BVQM_PC_BVQM_PC_UNCERTSS('Property','Value',...) creates a new BVQM_PC_BVQM_PC_UNCERTSS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bvqm_pc_bvqm_pc_uncertss_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bvqm_pc_bvqm_pc_uncertss_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bvqm_pc_bvqm_pc_uncertss

% Last Modified by GUIDE v2.5 03-May-2006 14:33:46

%
% WRITTEN BY : Mark A. McFarland, P.E.    303/497-4132
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bvqm_pc_uncertss_OpeningFcn, ...
                   'gui_OutputFcn',  @bvqm_pc_uncertss_OutputFcn, ...
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


% --- Executes just before bvqm_pc_bvqm_pc_uncertss is made visible.
function bvqm_pc_uncertss_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bvqm_pc_bvqm_pc_uncertss (see VARARGIN)

% Choose default command line output for bvqm_pc_bvqm_pc_uncertss
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Should the spatial scaling checkbox be displayed? - only if RR cal is to
% be performed:
qss=(varargin{1});
if qss==0
    set(handles.checkbox1, 'Enable', 'off');
    set(handles.checkbox1, 'Value' , 0);
else
    set(handles.checkbox1, 'Enable', 'on');
     set(handles.checkbox1, 'Value' , 0);
end

% UIWAIT makes bvqm_pc_bvqm_pc_uncertss wait for user response (see UIRESUME)
 uiwait(handles.bvqm_pc_uncertss);

% % Update handles structure
% guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = bvqm_pc_uncertss_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ~isempty(handles)
    varargout{1} = handles.output;
    delete(handles.bvqm_pc_uncertss);
else
    % assume 'x' was pressed on window
    handles.output = [];
    varargout{1}=handles.output;
%    delete(bvqm_pc_uncertss);
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
% perform_ss=get(handles.checkbox1, 'Value')
% handles.output.ss=perform_ss;
% guidata(hObject, handles);
% uiresume(handles.bvqm_pc_uncertss);

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
% handles.output.uncert=get(handles.edit1, 'Value')
% guidata(hObject, handles);
% uiresume(handles.bvqm_pc_uncertss);


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1. - OK button
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output=[];
perform_ss=get(handles.checkbox1, 'Value');
uncert=str2num(get(handles.edit1, 'String'));
if isempty(uncert) ; uncert=-1; end  % catch type errors

if uncert < 0
    h= errordlg ({'Value for temporal registration uncertainty must be at least 0.';
        'Please reenter.'}, 'Invalid value.', 'on');
    waitfor(h);
    pause(0.5);
else
    handles.output.ss=perform_ss;
    handles.output.uncert=uncert;
    guidata(hObject, handles);
    uiresume(handles.bvqm_pc_uncertss);
end






