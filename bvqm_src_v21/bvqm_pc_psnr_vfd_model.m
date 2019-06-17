function varargout = bvqm_pc_psnr_vfd_model(varargin)
% BVQM_PC_PSNR_VFD_MODEL MATLAB code for bvqm_pc_psnr_vfd_model.fig
%      BVQM_PC_PSNR_VFD_MODEL, by itself, creates a new BVQM_PC_PSNR_VFD_MODEL or raises the existing
%      singleton*.
%
%      H = BVQM_PC_PSNR_VFD_MODEL returns the handle to a new BVQM_PC_PSNR_VFD_MODEL or the handle to
%      the existing singleton*.
%
%      BVQM_PC_PSNR_VFD_MODEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BVQM_PC_PSNR_VFD_MODEL.M with the given input arguments.
%
%      BVQM_PC_PSNR_VFD_MODEL('Property','Value',...) creates a new BVQM_PC_PSNR_VFD_MODEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bvqm_pc_psnr_vfd_model_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bvqm_pc_psnr_vfd_model_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bvqm_pc_psnr_vfd_model

% Last Modified by GUIDE v2.5 02-Aug-2011 14:04:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bvqm_pc_psnr_vfd_model_OpeningFcn, ...
                   'gui_OutputFcn',  @bvqm_pc_psnr_vfd_model_OutputFcn, ...
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


% --- Executes just before bvqm_pc_psnr_vfd_model is made visible.
function bvqm_pc_psnr_vfd_model_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bvqm_pc_psnr_vfd_model (see VARARGIN)

% Choose default command line output for bvqm_pc_psnr_vfd_model
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.radiobutton1, 'Value' , 1);
set(handles.radiobutton2, 'Value' , 0);

% UIWAIT makes bvqm_pc_psnr_vfd_model wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = bvqm_pc_psnr_vfd_model_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ~isempty(handles)
    varargout{1} = handles.output;
    delete(handles.figure1);
else
    % assume 'x' was pressed on window
    handles.output = [];
    varargout{1}=handles.output;
%    delete(bvqm_pc_uncertss);
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
set(handles.radiobutton2, 'Value' , 0);


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
set(handles.radiobutton1, 'Value' , 0);



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output=[];
which_psnr(1)=get(handles.radiobutton1, 'Value');
which_psnr(2)=get(handles.radiobutton2, 'Value');
uncert=str2num(get(handles.edit1, 'String'));
if isempty(uncert) ; uncert=-1; end  % catch type errors

if(which_psnr(1) == 1)
    %PSNR of 255
    this_psnr = '255';
else
    %PSNR of 235
    this_psnr = '235';
end

if uncert < 0
    h= errordlg ({'Value for temporal uncertainty must be at least 0.';
        'Please reenter.'}, 'Invalid value.', 'on');
    waitfor(h);
    pause(0.5);
else
    handles.output.which_psnr=this_psnr;
    handles.output.uncert=uncert;
    guidata(hObject, handles);
    uiresume(handles.figure1);
end
