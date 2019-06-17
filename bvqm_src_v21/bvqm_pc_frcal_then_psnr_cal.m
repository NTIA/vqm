function varargout = bvqm_pc_frcal_then_psnr_cal(varargin)
% BVQM_PC_FRCAL_THEN_PSNR_CAL MATLAB code for bvqm_pc_frcal_then_psnr_cal.fig
%      BVQM_PC_FRCAL_THEN_PSNR_CAL, by itself, creates a new BVQM_PC_FRCAL_THEN_PSNR_CAL or raises the existing
%      singleton*.
%
%      H = BVQM_PC_FRCAL_THEN_PSNR_CAL returns the handle to a new BVQM_PC_FRCAL_THEN_PSNR_CAL or the handle to
%      the existing singleton*.
%
%      BVQM_PC_FRCAL_THEN_PSNR_CAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BVQM_PC_FRCAL_THEN_PSNR_CAL.M with the given input arguments.
%
%      BVQM_PC_FRCAL_THEN_PSNR_CAL('Property','Value',...) creates a new BVQM_PC_FRCAL_THEN_PSNR_CAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bvqm_pc_frcal_then_psnr_cal_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bvqm_pc_frcal_then_psnr_cal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bvqm_pc_frcal_then_psnr_cal

% Last Modified by GUIDE v2.5 03-Aug-2011 10:24:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bvqm_pc_frcal_then_psnr_cal_OpeningFcn, ...
                   'gui_OutputFcn',  @bvqm_pc_frcal_then_psnr_cal_OutputFcn, ...
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


% --- Executes just before bvqm_pc_frcal_then_psnr_cal is made visible.
function bvqm_pc_frcal_then_psnr_cal_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bvqm_pc_frcal_then_psnr_cal (see VARARGIN)

% Choose default command line output for bvqm_pc_frcal_then_psnr_cal
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes bvqm_pc_frcal_then_psnr_cal wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = bvqm_pc_frcal_then_psnr_cal_OutputFcn(hObject, eventdata, handles) 
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



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


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

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output=[];
frac_samp=str2num(get(handles.edit1, 'String'));
spat_x=str2num(get(handles.edit2, 'String'));
spat_y=str2num(get(handles.edit3, 'String'));
uncert=str2num(get(handles.edit4, 'String'));
frcal_uncert=str2num(get(handles.edit5, 'String'));
if isempty(uncert) ; uncert=-1; end  % catch type errors
if isempty(spat_y) ; spat_y=-1; end  % catch type errors
if isempty(spat_x) ; spat_x=-1; end  % catch type errors
if isempty(frac_samp) ; frac_samp=-1; end  % catch type errors
if isempty(frcal_uncert) ; frcal_uncert=-1; end  % catch type errors

if uncert < 0
    h= errordlg ({'Value for temporal uncertainty must be at least 0.';
        'Please reenter.'}, 'Invalid value.', 'on');
    waitfor(h);
    pause(0.5);
elseif(frcal_uncert < 0)
    h= errordlg ({'Value for temporal registration uncertainty must be at least 0.';
        'Please reenter.'}, 'Invalid value.', 'on');
    waitfor(h);
    pause(0.5);
elseif(frac_samp < .001 || frac_samp > 1)
    h= errordlg ({'Value for fraction sampled must be between .001 and 1.';
        'Please reenter.'}, 'Invalid value.', 'on');
    waitfor(h);
    pause(0.5);
elseif(spat_y < 0)
    h= errordlg ({'Value for spatial uncertainty in the y direction must be at least 0.';
        'Please reenter.'}, 'Invalid value.', 'on');
    waitfor(h);
    pause(0.5);
elseif(spat_x < 0)
    h= errordlg ({'Value for spatial uncertainty in the x direction must be at least 0.';
        'Please reenter.'}, 'Invalid value.', 'on');
    waitfor(h);
    pause(0.5);
else
    handles.output.frac_samp=frac_samp;
    handles.output.uncert=uncert;
    handles.output.spatial_y=spat_y;
    handles.output.spatial_x=spat_x;
    handles.output.frcal_uncert=frcal_uncert;
    guidata(hObject, handles);
    uiresume(handles.figure1);
end
