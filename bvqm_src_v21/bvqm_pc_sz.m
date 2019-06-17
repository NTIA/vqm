function varargout = bvqm_pc_sz(varargin)
% BVQM_PC_SZ M-file for bvqm_pc_sz.fig
%      BVQM_PC_SZ, by itself, creates a new BVQM_PC_SZ or raises the existing
%      singleton*.
%
%      H = BVQM_PC_SZ returns the handle to a new BVQM_PC_SZ or the handle to
%      the existing singleton*.
%
%      BVQM_PC_SZ('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BVQM_PC_SZ.M with the given input arguments.
%
%      BVQM_PC_SZ('Property','Value',...) creates a new BVQM_PC_SZ or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bvqm_pc_sz_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bvqm_pc_sz_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bvqm_pc_sz

% Last Modified by GUIDE v2.5 20-Jun-2006 13:35:50

%
% WRITTEN BY : Mark A. McFarland, P.E.    303/497-4132
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bvqm_pc_sz_OpeningFcn, ...
                   'gui_OutputFcn',  @bvqm_pc_sz_OutputFcn, ...
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


% --- Executes just before bvqm_pc_sz is made visible.
function bvqm_pc_sz_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bvqm_pc_sz (see VARARGIN)

% Choose default command line output for bvqm_pc_sz
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes bvqm_pc_sz wait for user response (see UIRESUME)
uiwait(handles.bvqm_pc_imgsz);


% --- Executes during object deletion, before destroying properties.
function bvqm_pc_sz_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to bvqm_pc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Remove the directory added to the path in the bvqm_pc_CreateFcn.
%delll=1

% ho=handles.output
% if ~ischar(ho)
%     handles.output='NONE'
% end

% --- Outputs from this function are returned to the command line.
function varargout = bvqm_pc_sz_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure


%ho_dum=handles.output

% if isempty(handles)
%     handles.output='NONE'
% %    guidata(hObject, handles);
%     %uiresume(handles.bvqm_pc_imgsz);
% end

varargout{1} = handles.output;

%delete(handles.bvqm_pc_imgsz);
%close(bvqm_pc_sz);
% button=handles.output;
delete(hObject);

% --- Executes on button press in NTSC.
function NTSC_Callback(hObject, eventdata, handles)
% hObject    handle to NTSC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output='NTSC';
guidata(hObject, handles);
uiresume(handles.bvqm_pc_imgsz);


% --- Executes on button press in PAL.
function PAL_Callback(hObject, eventdata, handles)
% hObject    handle to PAL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output='PAL';
guidata(hObject, handles);
uiresume(handles.bvqm_pc_imgsz);

% --- Executes on button press in HDTV720i.
function HDTV720i_Callback(hObject, eventdata, handles)
% hObject    handle to HDTV720i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output='HDTV720p';
guidata(hObject, handles);
uiresume(handles.bvqm_pc_imgsz);


% --- Executes on button press in HDTV1080i.
function HDTV1080i_Callback(hObject, eventdata, handles)
% hObject    handle to HDTV1080i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output='HDTV1080i';
guidata(hObject, handles);
uiresume(handles.bvqm_pc_imgsz);

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output='manual';
guidata(hObject, handles);
uiresume(handles.bvqm_pc_imgsz);



% --- Executes when user attempts to close bvqm_pc_imgsz.
function bvqm_pc_imgsz_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to bvqm_pc_imgsz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
ho=handles.output;
if ~ischar(ho)
    handles.output='NONE';
end

uiresume(hObject)


