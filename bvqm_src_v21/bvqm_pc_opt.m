function [varargout] = bvqm_pc_opt(varargin)
% BVQM_PC_OPT M-file for bvqm_pc_opt.fig
%      BVQM_PC_OPT, by itself, creates a new BVQM_PC_OPT or raises the existing
%      singleton*.
%
%      H = BVQM_PC_OPT returns the handle to a new BVQM_PC_OPT or the handle to
%      the existing singleton*.
%
%      BVQM_PC_OPT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BVQM_PC_OPT.M with the given input arguments.
%
%      BVQM_PC_OPT('Property','Value',...) creates a new BVQM_PC_OPT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bvqm_pc_opt_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bvqm_pc_opt_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bvqm_pc_opt

%
% WRITTEN BY : Mark A. McFarland, P.E.    303/497-4132
%


%Begin initialization code - DO NOT EDITl
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bvqm_pc_opt_OpeningFcn, ...
                   'gui_OutputFcn',  @bvqm_pc_opt_OutputFcn, ...
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


% --- Executes just before bvqm_pc_opt is made visible.
function bvqm_pc_opt_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bvqm_pc_opt (see VARARGIN)

% Choose default command line output for bvqm_pc_opt
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes bvqm_pc_opt wait for user response (see UIRESUME)
% uiwait(handles.bvqm_pc_opt);

%Populate boxes with known info::
% NOTE: varargin(1) = lb2_index
%       varargin(2) = lb2_filename
%       varargin(3) = num_frames
%       varargin(4) = loc_start
%       varargin(5) = loc_stop
%       varargin(6) = rows
%       varargin(7) = cols
%       varargin(8) = parsed_file flag (1=parsed)
%       varargin(9) = working_dir

% NOTE: some variables not visible - just used to store and access data
set(handles.listbox2_index_text22, 'String', varargin(1)); % lb2 index
set(handles.listbox2_filename_text23, 'String', varargin(2)); % lb2 fn
set(handles.parsed_file_text24, 'String', varargin(8)); % set parsed file flag
set(handles.working_dir_text25, 'String', varargin(9)); % working_dir
working_dir=varargin(9);

% If this is an 'original' clip, don't allow user to change all settings:
% disallow: gain offset, shift and scaling
clip_name=varargin(2);
[test_name, rem] = strtok(clip_name, '_');
[scene_name, rem] = strtok(rem, '_');
[hrc_name_temp, rem] = strtok(rem, '_');   % rem should be empty at this point.
[hrc_name,ext] = strtok(hrc_name_temp,'.');

if strcmp(hrc_name, 'original')
   set(handles.luninance_gain_edit5, 'Enable', 'off');
   set(handles.luninance_offset_edit6, 'Enable', 'off');
   set(handles.spatial_horiz_edit3, 'Enable', 'off'); 
   set(handles.spatial_vert_edit4, 'Enable', 'off');
   set(handles.horiz_scale_edit9, 'Enable', 'off');
   set(handles.vert_scale_edit10, 'Enable', 'off');
end




% Was optional data already entered? Look for saved bvqm_pc_opt .mat file:
lb2_index    = str2num(mat2str(cell2mat(get(handles.listbox2_index_text22, 'String'))));
lb2_filename = get(handles.listbox2_filename_text23, 'String');
[filename_pre, filename_ext] = strtok(lb2_filename, '.');
[filename_post] = strtok(filename_ext, '.');
if ispc ; path_sep = '\'; else; path_sep='/'; end
%opt_filename = mat2str(cell2mat(strcat(working_dir,path_sep, 'bvqm-opt_', filename_pre, '_', filename_post, '.mat')));
opt_filename = mat2str(cell2mat(strcat(working_dir,path_sep, 'bvqm-opt_', lb2_filename, '.mat')));
opt_exists = exist(opt_filename, 'file');
if opt_exists == 2  % Options were already entered, so reload them.
    % Load bvqm_pc_opt data from mat file:
    load(opt_filename);
    [dum, opt_index] = size(oclips);

    set(handles.start_frame_edit1, 'String', oclips(opt_index).loc_start);
    set(handles.stop_frame_edit2, 'String', oclips(opt_index).loc_stop);

    set(handles.spatial_horiz_edit3, 'String', oclips(opt_index).spatial.horizontal);
    set(handles.spatial_vert_edit4, 'String', oclips(opt_index).spatial.vertical);

    set(handles.luninance_gain_edit5, 'String', oclips(opt_index).luminance_gain);
    set(handles.luninance_offset_edit6, 'String', oclips(opt_index).luminance_offset);
    
    % Set aling_start/stop to NaN (or ''?) - 26May06...- NO- 05July06
     set(handles.align_start_edit7, 'String', oclips(opt_index).align_start);
     set(handles.align_stop_edit8, 'String', oclips(opt_index).align_stop);
%    set(handles.align_start_edit7, 'String', '');
%   set(handles.align_stop_edit8, 'String', '');
    
    set(handles.horiz_scale_edit9, 'String', oclips(opt_index).scale.horizontal);
    set(handles.vert_scale_edit10, 'String', oclips(opt_index).scale.vertical);
    
    set(handles.viewers_edit16, 'String', oclips(opt_index).viewers);
    set(handles.mos_edit15, 'String', oclips(opt_index).mos);
    set(handles.stdev_edit20, 'String', oclips(opt_index).stdev);
    set(handles.inlsa_mos_edit2, 'String', oclips(opt_index).inlsa_mos);
%    set(handles.inlsa_mos_2003_edit22, 'String', oclips(opt_index).inlsa_mos_2003);
    
    set(handles.cvr_top_edit11, 'String', oclips(opt_index).cvr.top);
    set(handles.cvr_bottom_edit12, 'String', oclips(opt_index).cvr.bottom);
    set(handles.cvr_left_edit13, 'String', oclips(opt_index).cvr.left);
    set(handles.cvr_right_edit14, 'String', oclips(opt_index).cvr.right);

    set(handles.subj_method_edit17, 'String', oclips(opt_index).subj_system);
    set(handles.hrc_desc_edit18, 'String', oclips(opt_index).hrc_definition);
    set(handles.scene_descr_edit19, 'String', oclips(opt_index).scene_definition);
    
else  % 1st time entering optional data -

    % Populate boxes with known info:
    set(handles.stop_frame_edit2, 'String', varargin(3));  % num_frames
    set(handles.start_frame_edit1, 'String', varargin(4)); % loc_start
    set(handles.stop_frame_edit2, 'String', varargin(5));  % loc_stop

    % Set align_start/stop to NaN (or ''?)- 26May06
    %set(handles.align_start_edit7, 'String', varargin(4)); % loc_start
    %set(handles.align_stop_edit8, 'String', varargin(5));  % loc_stop
    set(handles.align_start_edit7, 'String', ''); % loc_start
    set(handles.align_stop_edit8, 'String', '');  % loc_stop
   
    set(handles.cvr_bottom_edit12, 'String', varargin(6)); % image_size.rows
    set(handles.cvr_right_edit14, 'String', varargin(7));  % image_size.cols
   
    % set the cvrs to the values returned by the fcn default_sroi
    sz.rows=cell2mat(varargin(6));
    sz.cols=cell2mat(varargin(7));
    [ds] = default_sroi(sz);
    set(handles.cvr_top_edit11, 'String', ds.top);
    set(handles.cvr_bottom_edit12, 'String', ds.bottom);
    set(handles.cvr_left_edit13, 'String', ds.left);
    set(handles.cvr_right_edit14, 'String', ds.right);
end

% --- Outputs from this function are returned to the command line.
function varargout = bvqm_pc_opt_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handlistbox2_indexle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object deletion, before destroying properties.
function bvqm_pc_opt_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to bvqm_pc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%close_pushbutton1_Callback(hObject, eventdata, handles)


% --- Executes on button press in close_pushbutton1.
function close_pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to close_pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% save optional settings...

% oclips(lb2_index).video_standard = popupmenu1_string_list{get(handles.video_std_popupmenu1, 'Value')};
% oclips(lb2_index).fps            = get(handles.fps_edit, 'String');

% [oclips_filename, rem] = strtok(lb2_filename, '@');
% oclips(lb2_index).file_name     = {oclips_filename};  % CELL

%lb2_index    = str2num(mat2str(cell2mat(get(handles.listbox2_index_text22, 'String'))));
lb2_index    = str2num(cell2mat(get(handles.listbox2_index_text22, 'String')));
lb2_filename = get(handles.listbox2_filename_text23, 'String');

% oclips(lb2_index).loc_start = str2num(get(handles.start_frame_edit1, 'String'));
% oclips(lb2_index).loc_stop  = str2num(get(handles.stop_frame_edit2, 'String'));
oclips(lb2_index).loc_start = get(handles.start_frame_edit1, 'String');
oclips(lb2_index).loc_stop  = get(handles.stop_frame_edit2, 'String');


spatial.horizontal = str2num(get(handles.spatial_horiz_edit3, 'String')); 
spatial.vertical   = str2num(get(handles.spatial_vert_edit4, 'String'));
%spatial.horizontal = get(handles.spatial_horiz_edit3, 'String'); 
%spatial.vertical   = get(handles.spatial_vert_edit4, 'String');
oclips(lb2_index).spatial = spatial;

oclips(lb2_index).luminance_gain   = str2num(get(handles.luninance_gain_edit5, 'String'));
oclips(lb2_index).luminance_offset = str2num(get(handles.luninance_offset_edit6, 'String'));

% oclips(lb2_index).align_start = str2num(cell2mat(get(handles.align_start_edit7, 'String')));
% oclips(lb2_index).align_stop  = str2num(cell2mat(get(handles.align_stop_edit8, 'String')));
oclips(lb2_index).align_start = get(handles.align_start_edit7, 'String');
oclips(lb2_index).align_stop  = get(handles.align_stop_edit8, 'String');

scale.horizontal = str2num(get(handles.horiz_scale_edit9, 'String'));
scale.vertical   = str2num(get(handles.vert_scale_edit10, 'String'));
oclips(lb2_index).scale = scale;  % STRUCT

oclips(lb2_index).viewers     = str2num(get(handles.viewers_edit16, 'String')); 
oclips(lb2_index).mos = str2num(get(handles.mos_edit15, 'String'));
oclips(lb2_index).stdev = str2num(get(handles.stdev_edit20, 'String'));
oclips(lb2_index).inlsa_mos = str2num(get(handles.inlsa_mos_edit2, 'String'));

cvr.top    = str2num(get(handles.cvr_top_edit11, 'String'));
cvr.bottom = str2num(get(handles.cvr_bottom_edit12, 'String'));
cvr.left   = str2num(get(handles.cvr_left_edit13, 'String'));
cvr.right  = str2num(get(handles.cvr_right_edit14, 'String'));
oclips(lb2_index).cvr = cvr;

oclips(lb2_index).subj_system = get(handles.subj_method_edit17, 'String');  % CELL
oclips(lb2_index).hrc_definition   = get(handles.hrc_desc_edit18, 'String'); %CELL
oclips(lb2_index).scene_definition = get(handles.scene_descr_edit19, 'String');  % CELL

% Create new filename and save information in .mat file:
[filename_pre, filename_ext] = strtok(lb2_filename, '.');
[filename_post] = strtok(filename_ext, '.');

working_dir=cell2mat(get(handles.working_dir_text25, 'String'));
if ispc ; path_sep = '\'; else; path_sep='/'; end
%savefilename=mat2str(cell2mat(strcat('bvqm-opt_', filename_pre,  filename_ext, '.mat')));
savefilename=cell2mat(strcat('bvqm-opt_', filename_pre,  filename_ext, '.mat'));
%save_file=mat2str(cell2mat(strcat(working_dir, path_sep, savefilename)));
save_file=fullfile(working_dir, savefilename);

save(save_file, 'oclips');
close;
pause(1);

function start_frame_edit1_Callback(hObject, eventdata, handles)
% hObject    handle to start_frame_edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_frame_edit1 as text
%        str2double(get(hObject,'String')) returns contents of start_frame_edit1 as a double


% --- Executes during object creation, after setting all properties.
function start_frame_edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_frame_edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function stop_frame_edit2_Callback(hObject, eventdata, handles)
% hObject    handle to stop_frame_edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stop_frame_edit2 as text
%        str2double(get(hObject,'String')) returns contents of stop_frame_edit2 as a double

%set(handles.stop_frame_edit2, 'String', num_frames)

% --- Executes during object creation, after setting all properties.
function stop_frame_edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stop_frame_edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function spatial_horiz_edit3_Callback(hObject, eventdata, handles)
% hObject    handle to spatial_horiz_edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spatial_horiz_edit3 as text
%        str2double(get(hObject,'String')) returns contents of spatial_horiz_edit3 as a double


% --- Executes during object creation, after setting all properties.
function spatial_horiz_edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spatial_horiz_edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function spatial_vert_edit4_Callback(hObject, eventdata, handles)
% hObject    handle to spatial_vert_edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spatial_vert_edit4 as text
%        str2double(get(hObject,'String')) returns contents of spatial_vert_edit4 as a double


% --- Executes during object creation, after setting all properties.
function spatial_vert_edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spatial_vert_edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function luninance_gain_edit5_Callback(hObject, eventdata, handles)
% hObject    handle to luninance_gain_edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of luninance_gain_edit5 as text
%        str2double(get(hObject,'String')) returns contents of luninance_gain_edit5 as a double


% --- Executes during object creation, after setting all properties.
function luninance_gain_edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to luninance_gain_edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function luninance_offset_edit6_Callback(hObject, eventdata, handles)
% hObject    handle to luninance_offset_edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of luninance_offset_edit6 as text
%        str2double(get(hObject,'String')) returns contents of luninance_offset_edit6 as a double


% --- Executes during object creation, after setting all properties.
function luninance_offset_edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to luninance_offset_edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function align_start_edit7_Callback(hObject, eventdata, handles)
% hObject    handle to align_start_edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of align_start_edit7 as text
%        str2double(get(hObject,'String')) returns contents of align_start_edit7 as a double

% --- Executes during object creation, after setting all properties.
function align_start_edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to align_start_edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function align_stop_edit8_Callback(hObject, eventdata, handles)
% hObject    handle to align_stop_edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of align_stop_edit8 as text
%        str2double(get(hObject,'String')) returns contents of align_stop_edit8 as a double


% --- Executes during object creation, after setting all properties.
function align_stop_edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to align_stop_edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function horiz_scale_edit9_Callback(hObject, eventdata, handles)
% hObject    handle to horiz_scale_edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of horiz_scale_edit9 as text
%        str2double(get(hObject,'String')) returns contents of horiz_scale_edit9 as a double


% --- Executes during object creation, after setting all properties.
function horiz_scale_edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to horiz_scale_edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function vert_scale_edit10_Callback(hObject, eventdata, handles)
% hObject    handle to vert_scale_edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vert_scale_edit10 as text
%        str2double(get(hObject,'String')) returns contents of vert_scale_edit10 as a double


% --- Executes during object creation, after setting all properties.
function vert_scale_edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vert_scale_edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cvr_top_edit11_Callback(hObject, eventdata, handles)
% hObject    handle to cvr_top_edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cvr_top_edit11 as text
%        str2double(get(hObject,'String')) returns contents of cvr_top_edit11 as a double


% --- Executes during object creation, after setting all properties.
function cvr_top_edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cvr_top_edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function cvr_bottom_edit12_Callback(hObject, eventdata, handles)
% hObject    handle to cvr_bottom_edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cvr_bottom_edit12 as text
%        str2double(get(hObject,'String')) returns contents of cvr_bottom_edit12 as a double


% --- Executes during object creation, after setting all properties.
function cvr_bottom_edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cvr_bottom_edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function cvr_left_edit13_Callback(hObject, eventdata, handles)
% hObject    handle to cvr_left_edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cvr_left_edit13 as text
%        str2double(get(hObject,'String')) returns contents of cvr_left_edit13 as a double


% --- Executes during object creation, after setting all properties.
function cvr_left_edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cvr_left_edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function cvr_right_edit14_Callback(hObject, eventdata, handles)
% hObject    handle to cvr_right_edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cvr_right_edit14 as text
%        str2double(get(hObject,'String')) returns contents of cvr_right_edit14 as a double


% --- Executes during object creation, after setting all properties.
function cvr_right_edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cvr_right_edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function mos_edit15_Callback(hObject, eventdata, handles)
% hObject    handle to mos_edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mos_edit15 as text
%        str2double(get(hObject,'String')) returns contents of mos_edit15 as a double


% --- Executes during object creation, after setting all properties.
function mos_edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mos_edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function viewers_edit16_Callback(hObject, eventdata, handles)
% hObject    handle to viewers_edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of viewers_edit16 as text
%        str2double(get(hObject,'String')) returns contents of viewers_edit16 as a double


% --- Executes during object creation, after setting all properties.
function viewers_edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to viewers_edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function subj_method_edit17_Callback(hObject, eventdata, handles)
% hObject    handle to subj_method_edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subj_method_edit17 as text
%        str2double(get(hObject,'String')) returns contents of subj_method_edit17 as a double

%get(handles.bvqm_pc_opt

% --- Executes during object creation, after setting all properties.
function subj_method_edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subj_method_edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function hrc_desc_edit18_Callback(hObject, eventdata, handles)
% hObject    handle to hrc_desc_edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hrc_desc_edit18 as text
%        str2double(get(hObject,'String')) returns contents of hrc_desc_edit18 as a double


% --- Executes during object creation, after setting all properties.
function hrc_desc_edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hrc_desc_edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function scene_descr_edit19_Callback(hObject, eventdata, handles)
% hObject    handle to scene_descr_edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scene_descr_edit19 as text
%        str2double(get(hObject,'String')) returns contents of scene_descr_edit19 as a double


% --- Executes during object creation, after setting all properties.
function scene_descr_edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scene_descr_edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function stdev_edit20_Callback(hObject, eventdata, handles)
% hObject    handle to stdev_edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stdev_edit20 as text
%        str2double(get(hObject,'String')) returns contents of stdev_edit20 as a double


% --- Executes during object creation, after setting all properties.
function stdev_edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stdev_edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function inlsa_mos_edit2_Callback(hObject, eventdata, handles)
% hObject    handle to inlsa_mos_edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inlsa_mos_edit2 as text
%        str2double(get(hObject,'String')) returns contents of inlsa_mos_edit2 as a double


% --- Executes during object creation, after setting all properties.
function inlsa_mos_edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inlsa_mos_edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function inlsa_mos_2003_edit22_Callback(hObject, eventdata, handles)
% hObject    handle to inlsa_mos_2003_edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inlsa_mos_2003_edit22 as text
%        str2double(get(hObject,'String')) returns contents of inlsa_mos_2003_edit22 as a double


% --- Executes during object creation, after setting all properties.
function inlsa_mos_2003_edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inlsa_mos_2003_edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


