function varargout = bvqm_pc_cal(varargin)
% BVQM_PC_CAL M-file for bvqm_pc_cal.fig
%      BVQM_PC_CAL, by itself, creates a new BVQM_PC_CAL or raises the existing
%      singleton*.
%
%      H = BVQM_PwC_CAL returns the handle to a new BVQM_PC_CAL or the handle to
%      the existing singleton*.
%
%      BVQM_PC_CAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BVQM_PC_CAL.M with the given input arguments.
%
%      BVQM_PC_CAL('Property','Value',...) creates a new BVQM_PC_CAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bvqm_pc_cal_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bvqm_pc_cal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bvqm_pc_cal



% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @bvqm_pc_cal_OpeningFcn, ...
    'gui_OutputFcn',  @bvqm_pc_cal_OutputFcn, ...
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


% --- Executes just before bvqm_pc_cal is made visible.
function bvqm_pc_cal_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bvqm_pc_cal (see VARARGIN)

% Choose default command line output for bvqm_pc_cal
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% set red. ref v2 cal as default:
set(handles.radiobutton2, 'Value', 0)  %full ref
set(handles.radiobutton3, 'Value', 0)  %red. ref
set(handles.radiobutton4, 'Value', 0)  % no cal
set(handles.radiobutton1, 'Value', 0)  % manual cal
set(handles.radiobutton12, 'Value', 0) % TR,VR-FR or RR
set(handles.radiobutton13, 'Value', 0) % manual cal then TR
set(handles.radiobutton15, 'Value', 1) % red. ref, version 2
set(handles.radiobutton16, 'Value', 0) % PSNR
set(handles.radiobutton17, 'Value', 0) % rrcal v2 then PSNR
set(handles.radiobutton18, 'Value', 0) % full ref cal then PSNR

% crash_rec=(varargin{2});

% %Set-up working directory acces:
if ispc ; path_sep = '\'; else; path_sep='/'; end
passed=varargin{2};
working_dir=passed{1};
crash_rec=passed{2};

working_dir=fullfile(working_dir, path_sep);

set(handles.text5_WorkingDirectory, 'String', working_dir);

jcpf=strcat(working_dir, path_sep, 'jclips.mat');  %jclips path/file name
load(jcpf);

%NOTE: set the clip_dir to the first entry in jtests.path
set(handles.text7_clip_dir, 'String', jtests(1).path);

% UIWAIT makes bvqm_pc_cal wait for user response (see UIRESUME)
% uiwait(handles.bvqm_pc_cal);

% Crash Recovery - Reset screen to prior selections, then process clips:
if isequal(crash_rec, 1)
    crashed=1;
    % unset default calibration selection:
    set(handles.radiobutton2, 'Value', 0);
    h=msgbox('Press OK to continue with crash recovery.', 'Automatic Crash Recovery', 'help');
    waitfor(h);
    % continue processing from crash:
    pushbutton1_Continue_Callback(hObject, eventdata, handles);
end


function bvqm_pc_cal_passer(pass)
% use command 'findobj' to find object handles....

handles=guidata('bvqm_pc_cal')

%Set-up working directory acces:
if ispc ; path_sep = '\'; else; path_sep='/'; end
passed=varargin{1};
working_dir=passed{1};
crash_rec=passed{2};

working_dir=fullfile(working_dir, path_sep);

set(handles.text5_WorkingDirectory, 'String', working_dir);

jcpf=strcat(working_dir, path_sep, 'jclips.mat');  %jclips path/file name
load(jcpf);

%NOTE: set the clip_dir to the first entry in jtests.path; this will work
%because jtest.path is the same for all clips...
set(handles.text7_clip_dir, 'String', jtests(1).path);

% UIWAIT makes bvqm_pc_cal wait for user response (see UIRESUME)
% uiwait(handles.bvqm_pc_cal);

% Crash Recovery - Reset screen to prior selections, then process clips:
if isequal(crash_rec, 1)
    crashed=1;
    % unset default calibration selection:
    set(handles.radiobutton2, 'Value', 0);
    h=msgbox('Press OK to continue with crash recovery.', 'Automatic Crash Recovery', 'help');
    waitfor(h);
    % continue processing from crash:
    pushbutton1_Continue_Callback(hObject, eventdata, handles);
end

%handles.working_dir=(varargin{2});
    
% --- Outputs from this function are returned to the command line.
function varargout = bvqm_pc_cal_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;

% --- Executes when user attempts to close bvqm_pc.
function bvqm_pc_cal_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to bvqm_pc_imgsz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Hint: delete(hObject) closes the figure
warning off all;

pushbutton4_exit_Callback(hObject, eventdata, handles)

% --- Executes on button press in radiobutton1.
% --- Manual Calibration Radiobattun
function radiobutton1_Callback(hObject, eventdata, handles) % MANUAL CAL
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1

set(handles.radiobutton2, 'Value', 0)
set(handles.radiobutton3, 'Value', 0)
set(handles.radiobutton4, 'Value', 0)
set(handles.radiobutton12, 'Value', 0)
set(handles.radiobutton13, 'Value', 0)
set(handles.radiobutton15, 'Value', 0)
set(handles.radiobutton16, 'Value', 0)
set(handles.radiobutton17, 'Value', 0)
set(handles.radiobutton18, 'Value', 0)


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles) % FULL REF CAL
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
set(handles.radiobutton1, 'Value', 0)
set(handles.radiobutton3, 'Value', 0)
set(handles.radiobutton4, 'Value', 0)
set(handles.radiobutton12, 'Value', 0)
set(handles.radiobutton13, 'Value', 0)
set(handles.radiobutton15, 'Value', 0)
set(handles.radiobutton16, 'Value', 0)
set(handles.radiobutton17, 'Value', 0)
set(handles.radiobutton18, 'Value', 0)


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles) % RED. REF CAL
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
set(handles.radiobutton1, 'Value', 0)
set(handles.radiobutton2, 'Value', 0)
set(handles.radiobutton4, 'Value', 0)
set(handles.radiobutton12, 'Value', 0)
set(handles.radiobutton13, 'Value', 0)
set(handles.radiobutton15, 'Value', 0)
set(handles.radiobutton16, 'Value', 0)
set(handles.radiobutton17, 'Value', 0)
set(handles.radiobutton18, 'Value', 0)


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles) % NO CAL
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
set(handles.radiobutton1, 'Value', 0)
set(handles.radiobutton2, 'Value', 0)
set(handles.radiobutton3, 'Value', 0)
set(handles.radiobutton12, 'Value', 0)
set(handles.radiobutton13, 'Value', 0)
set(handles.radiobutton15, 'Value', 0)
set(handles.radiobutton16, 'Value', 0)
set(handles.radiobutton17, 'Value', 0)
set(handles.radiobutton18, 'Value', 0)


% --- Executes on button press in radiobutton12.
function radiobutton12_Callback(hObject, eventdata, handles) % TR,VR- FR
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
set(handles.radiobutton1, 'Value', 0)
set(handles.radiobutton2, 'Value', 0)
set(handles.radiobutton3, 'Value', 0)
set(handles.radiobutton4, 'Value', 0)
set(handles.radiobutton13, 'Value', 0)
set(handles.radiobutton15, 'Value', 0)
set(handles.radiobutton16, 'Value', 0)
set(handles.radiobutton17, 'Value', 0)
set(handles.radiobutton18, 'Value', 0)


% --- Executes on button press in radiobutton13.
function radiobutton13_Callback(hObject, eventdata, handles) % TR,VR - RR
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
set(handles.radiobutton1, 'Value', 0)
set(handles.radiobutton2, 'Value', 0)
set(handles.radiobutton3, 'Value', 0)
set(handles.radiobutton4, 'Value', 0)
set(handles.radiobutton12, 'Value', 0)
set(handles.radiobutton15, 'Value', 0)
set(handles.radiobutton16, 'Value', 0)
set(handles.radiobutton17, 'Value', 0)
set(handles.radiobutton18, 'Value', 0)

% --- Executes on button press in radiobutton13.
function radiobutton15_Callback(hObject, eventdata, handles) % TR,VR - RR
% hObject    handle to radiobutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
set(handles.radiobutton1, 'Value', 0)
set(handles.radiobutton2, 'Value', 0)
set(handles.radiobutton3, 'Value', 0)
set(handles.radiobutton4, 'Value', 0)
set(handles.radiobutton12, 'Value', 0)
set(handles.radiobutton13, 'Value', 0)
set(handles.radiobutton16, 'Value', 0)
set(handles.radiobutton17, 'Value', 0)
set(handles.radiobutton18, 'Value', 0)

function radiobutton16_Callback(hObject, eventdata, handles) % TR,VR - RR
% hObject    handle to radiobutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
set(handles.radiobutton1, 'Value', 0)
set(handles.radiobutton2, 'Value', 0)
set(handles.radiobutton3, 'Value', 0)
set(handles.radiobutton4, 'Value', 0)
set(handles.radiobutton12, 'Value', 0)
set(handles.radiobutton13, 'Value', 0)
set(handles.radiobutton15, 'Value', 0)
set(handles.radiobutton17, 'Value', 0)
set(handles.radiobutton18, 'Value', 0)

function radiobutton17_Callback(hObject, eventdata, handles) % TR,VR - RR
% hObject    handle to radiobutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
set(handles.radiobutton1, 'Value', 0)
set(handles.radiobutton2, 'Value', 0)
set(handles.radiobutton3, 'Value', 0)
set(handles.radiobutton4, 'Value', 0)
set(handles.radiobutton12, 'Value', 0)
set(handles.radiobutton13, 'Value', 0)
set(handles.radiobutton15, 'Value', 0)
set(handles.radiobutton16, 'Value', 0)
set(handles.radiobutton18, 'Value', 0)

function radiobutton18_Callback(hObject, eventdata, handles) % TR,VR - RR
% hObject    handle to radiobutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
set(handles.radiobutton1, 'Value', 0)
set(handles.radiobutton2, 'Value', 0)
set(handles.radiobutton3, 'Value', 0)
set(handles.radiobutton4, 'Value', 0)
set(handles.radiobutton12, 'Value', 0)
set(handles.radiobutton13, 'Value', 0)
set(handles.radiobutton15, 'Value', 0)
set(handles.radiobutton16, 'Value', 0)
set(handles.radiobutton17, 'Value', 0)


% --- Executes on button press in pushbutton1_Continue.
function pushbutton1_Continue_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1_Continue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% initialize directory variables and file names.
if ispc ; path_sep = '\'; else; path_sep='/'; end
working_dir = get(handles.text5_WorkingDirectory, 'String');
clip_dir = mat2str(cell2mat(get(handles.text7_clip_dir, 'String')));
temp_file = fullfile(working_dir, 'bvqm-temp_file');
status_file=fullfile(working_dir, 'bvqm-status.mat');
lgo_status=[];  % luminance gain offset status initialization

jcpf=strcat(working_dir, path_sep, 'jclips.mat');  %jclips path/file name
status =[];

load(jcpf); % load jclips

if exist(status_file, 'file') % if status file exists, gather info from file...
    load (status_file);
    crash=1;

    if strcmp(cal_type, 'No Calibration')
        no_cal=1;
        manual_cal=0;
        full_ref_cal=0;
        reduced_ref_cal=0;
        reduced_ref2_cal=0;
        mc_then_tr=0;
        tr_vr=0;
        psnr = 0;
        rrcal2_then_psnr = 0;
        full_ref_cal_then_psnr = 0;
%         set(handles.radiobutton4, 'Value', 1)  % no cal
    end

    if strcmp(cal_type, 'Manual Calibration')
        no_cal=0;
        manual_cal=1;
        full_ref_cal=0;
        reduced_ref_cal=0;
        reduced_ref2_cal=0;
        mc_then_tr=0;
        tr_vr=0;
        psnr = 0;
        rrcal2_then_psnr = 0;
        full_ref_cal_then_psnr = 0;
%         set(handles.radiobutton1, 'Value', 0)  % manual cal
    end

    if strncmp(cal_type, 'Full Reference Calibration',26)
        no_cal=0;
        manual_cal=0;
        full_ref_cal=1;
        reduced_ref_cal=0;
        reduced_ref2_cal=0;
        mc_then_tr=0;
        tr_vr=0;
        psnr = 0;
        rrcal2_then_psnr = 0;
        full_ref_cal_then_psnr = 0;
%         set(handles.radiobutton2, 'Value', 1)  %full ref
    end

    if strcmp(cal_type, 'Reduced Reference Calibration Version 1')
        no_cal=0;
        manual_cal=0;
        full_ref_cal=0;
        reduced_ref_cal=1;
        reduced_ref2_cal=0;
        mc_then_tr=0;
        tr_vr=0;
        psnr = 0;
        rrcal2_then_psnr = 0;
        full_ref_cal_then_psnr = 0;
%         set(handles.radiobutton3, 'Value', 0)  %red. ref

    end
 
    if strcmp(cal_type, 'Reduced Reference Calibration Version 2')
        no_cal=0;
        manual_cal=0;
        full_ref_cal=0;
        reduced_ref_cal=0;
        reduced_ref2_cal=1;
        mc_then_tr=0;
        tr_vr=0;
        psnr = 0;
        rrcal2_then_psnr = 0;
        full_ref_cal_then_psnr = 0;
%         set(handles.radiobutton3, 'Value', 0)  %red. ref

    end

    if strcmp(cal_type,'Temporal Registration and Valid Region Only');
        no_cal=0;
        manual_cal=0;
        full_ref_cal=0;
        reduced_ref_cal=0;
        reduced_ref2_cal=0;
        mc_then_tr=0;
        tr_vr=1;
        psnr = 0;
        rrcal2_then_psnr = 0;
        full_ref_cal_then_psnr = 0;
%         set(handles.radiobutton3, 'Value', 0)  %red. ref

    end

    if strcmp(cal_type,'Temporal Registration after Manual Calibration');
        no_cal=0;
        manual_cal=0;
        full_ref_cal=0;
        reduced_ref_cal=0;
        reduced_ref2_cal=0;
        mc_then_tr=1;
        tr_vr=0;
        psnr = 0;
        rrcal2_then_psnr = 0;
        full_ref_cal_then_psnr = 0;
%         set(handles.radiobutton3, 'Value', 0)  %red. ref
    end
    
    if strncmp(cal_type,'Peak Signal to Noise Ratio',26);
        no_cal=0;
        manual_cal=0;
        full_ref_cal=0;
        reduced_ref_cal=0;
        reduced_ref2_cal=0;
        mc_then_tr=0;
        tr_vr=0;
        psnr = 1;
        rrcal2_then_psnr = 0;
        full_ref_cal_then_psnr = 0;
%         set(handles.radiobutton3, 'Value', 0)  %red. ref
    end
    
    if strcmp(cal_type,'ITU-T J.244 followed by ITU-T J.340');
        no_cal=0;
        manual_cal=0;
        full_ref_cal=0;
        reduced_ref_cal=0;
        reduced_ref2_cal=0;
        mc_then_tr=0;
        tr_vr=0;
        psnr = 0;
        rrcal2_then_psnr = 1;
        full_ref_cal_then_psnr = 0;
%         set(handles.radiobutton3, 'Value', 0)  %red. ref
    end
    
    if strcmp(cal_type,'ITU-T J.144 followed by ITU-T J.340');
        no_cal=0;
        manual_cal=0;
        full_ref_cal=0;
        reduced_ref_cal=0;
        reduced_ref2_cal=0;
        mc_then_tr=0;
        tr_vr=0;
        psnr = 0;
        rrcal2_then_psnr = 0;
        full_ref_cal_then_psnr = 1;
%         set(handles.radiobutton3, 'Value', 0)  %red. ref
    end


else % if status file doesn't exist, get info from user/gui

    % Which type of cal was requested?:
    manual_cal = get(handles.radiobutton1, 'Value');
    full_ref_cal = get(handles.radiobutton2, 'Value');
    reduced_ref_cal = get(handles.radiobutton3, 'Value');
    no_cal = get(handles.radiobutton4, 'Value');
    tr_vr = get(handles.radiobutton12, 'Value');
    mc_then_tr = get(handles.radiobutton13, 'Value');
    reduced_ref2_cal = get(handles.radiobutton15, 'Value');
    psnr = get(handles.radiobutton16, 'Value');
    rrcal2_then_psnr = get(handles.radiobutton17, 'Value');
    full_ref_cal_then_psnr = get(handles.radiobutton18, 'Value');

    % Which model will be run?
    model_list = get(handles.listbox1, 'String');
     index = get(handles.listbox1, 'Value');
    model_selected = model_list(index);
    % IMPT: order of its_model_list below must match the order in the listbox!
    its_model_list = {'model_general'; 'model_lowbw'; 'model_developers'; 'model_television'; 'model_videoconferencing'; 'model_fastlowbw'; 'model_psnr'; 'model_psnr_vfd'; 'model_vqm_vfd'};
    model_to_run = cell2mat(its_model_list(index)); % ITS matlab model name

    % % Will spatial scaling be required? 
    ssf=0;
end


% for video clips with no defined aligned segments, assign align_start and
% align_stop to loc_start and loc_stop:
for count = 1:size(jclips,2)
    if isnan(jclips(count).align_start)
        jclips(count).align_start = jclips(count).loc_start;
    end

    if isnan(jclips(count).align_stop)
        jclips(count).align_stop = jclips(count).loc_stop;
    end
end

%MOVE TO A NEW FUNCTION TO CALL AFTER CALIBRATION QUESTIONS ARE DONE!
% % if chose PSNR, determine peak white value.
% if strcmp(model_to_run, 'model_psnr')
%     [peak_value] = questdlg('What would you like to use for the peak signal:  255 (computer white) or 235 (ITU-R Rec. BT-601 white)?', 'PSNR Model', '255', '235', '235');
% 
%     if strcmp(peak_value,'235'),
%         model_to_run = 'psnr_235';
%     else
%         model_to_run = 'psnr_255';
%     end
% end
% 
% model_additional_options = [];
% if(strcmp(model_to_run, 'model_psnr_vfd'))
%     [output]=bvqm_pc_psnr_vfd_model();
%     if isempty(output)
%         return
%     else
%         model_additional_options.which_psnr=output.which_psnr;
%         model_additional_options.t_uncert=output.uncert ;
%     end
% end
% 
% if((strcmpmodel_to_run, 'model_vqm_vfd'))
%     [output]=bvqm_pc_vqm_vfd_model(jclips(1).image_size);
%     if isempty(output)
%         return
%     else
%         viewing_distance=output.vdist;
%         model_additional_options.t_uncert=output.uncert ;
%         jtests.viewing_distance = viewing_distance;
%     end
% end
model_additional_options = [];

% output = model_gui(model_to_run);
% if strcmp(model_to_run, 'model_psnr')
%     model_to_run = output.model_to_run;
% elseif(strcmp(model_to_run, 'model_psnr_vfd'))
%     model_additional_options.which_psnr=output.model_additional_options.which_psnr;
%     model_additional_options.t_uncert=output.model_additional_options.t_uncert;
% elseif(strcmp(model_to_run, 'model_vqm_vfd'))
%     model_additional_options.t_uncert=output.model_additional_options.t_uncert ;
%     jtests.viewing_distance = output.jtests.viewing_distance;
% end


%--------------------------------------
% ------------ NO CAL -----------------
%--------------------------------------
if (no_cal )

    %Close Calibration window, begin cal
    pause(0.25);
%    delete (handles.bvqm_pc_cal);

    cal_type='No Calibration'; scal_type='no_cal';
    
    % NOTE: var model_to_run exists at this point!
    if exist(status_file, 'file')
        load (status_file);
        crash=1;
    else
        output = model_gui(model_to_run,jclips);
        % logic to handle if 'x' was pressed
        if isempty(output)
            return
        end
        if strcmp(model_to_run, 'model_psnr')
            model_to_run = output.model_to_run;
        elseif(strcmp(model_to_run, 'model_psnr_vfd'))
            model_additional_options.which_psnr=output.model_additional_options.which_psnr;
            model_additional_options.t_uncert=output.model_additional_options.t_uncert;
        elseif(strcmp(model_to_run, 'model_vqm_vfd'))
            model_additional_options.t_uncert=output.model_additional_options.t_uncert ;
            jtests.viewing_distance = output.jtests.viewing_distance;
        end
    end

    % set-up waitbar:
    wb = waitbar(0,'BVQM is running. Please wait...');
    pause(0.25); delete (handles.bvqm_pc_cal);
    waitbar(1/4); pause(0.25);
    
    % Replace all optional/advanced values with the default defaults, then run model:

    % Loop through each item in listbox2 to update clips structure:
    % listbox 2 is also stored in jclips.filename::

    lb2_index =1;
    lb2_size=size(jclips,2);

    lb2_filelist=[jclips.file_name];
    %  lb2_filelist=get(handles.lb2, 'String');

    % Assign all clips default values.
    while lb2_index <= lb2_size

        lb2_filename = lb2_filelist{lb2_index};
%         clip_dir = get(handles.text7_clip_dir, 'String');
%        clip_file=mat2str(cell2mat(strcat(clip_dir, lb2_filename)));
        clip_file=fullfile(clip_dir, lb2_filename);
        [fid,message] = fopen(clip_file);

        if fid ~= -1
            fseek(fid,0,'eof');
            file_size = ftell(fid);
            rows = jclips(lb2_index).image_size.rows;
            cols = jclips(lb2_index).image_size.cols;
%             % note: floor needed b/c sometimes division will result in fractional frames...
%             num_frames = floor(file_size/(rows*cols*2));  % XXX True for all files or just yuv files?
             num_frames = file_size/(rows*cols*2);
        else
            file_size = NaN;
            rows = NaN;
            cols = NaN;
            num_frames = NaN;
        end

        % SET DEFAULTS HERE:
        %        [test_name, scene_name, hrc_name, ext] = parse_clip(lb2_filename,handles);

        %         jclips(lb2_index).test  = {test_name};
        %         jclips(lb2_index).scene = {scene_name};
        %         jclips(lb2_index).hrc   = {hrc_name};
        %         jclips(lb2_index).image_size.rows = str2num(get(handles.image_size_rows_edit, 'String'));
        %         jclips(lb2_index).image_size.cols = str2num(get(handles.image_size_cols_edit, 'String'));
        jclips(lb2_index).spatial.horizontal = 0;
        jclips(lb2_index).spatial.vertical   = 0;
        jclips(lb2_index).luminance_gain   = 1.000;
        jclips(lb2_index).luminance_offset = 0.0000;
        %         jclips(lb2_index).video_standard = jclips(lb2_index).video_standard;
        %         jclips(lb2_index).fps            = str2num(get(handles.fps_edit, 'String'));
        %         jclips(lb2_index).subj_system   = {''};  % XXX don't know!
        %         jclips(lb2_index).mos         = NaN;
        %         jclips(lb2_index).stdev       = NaN;
        %         jclips(lb2_index).inlsa_mos = NaN;

        [jclips_filename, rem] = strtok(lb2_filename, '@');
        jclips(lb2_index).file_name     = {jclips_filename};  % CELL

 
        
 %       jclips(lb2_index).loc_start = 1; % XXX true for parsed ? no!
 %       jclips(lb2_index).loc_stop  = num_frames;

        jclips(lb2_index).align_start = jclips(lb2_index).loc_start;
        jclips(lb2_index).align_stop  = jclips(lb2_index).loc_stop;

        %         jclips(lb2_index).cvr.top     = 1;
        %         jclips(lb2_index).cvr.bottom  = str2num(get(handles.image_size_rows_edit, 'String'));
        %         jclips(lb2_index).cvr.left    = 1;
        %         jclips(lb2_index).cvr.right   = str2num(get(handles.image_size_cols_edit, 'String'));
        %         jclips(lb2_index).viewers     = [];
        %         jclips(lb2_index).hrc_definition   = {''};
        %         jclips(lb2_index).scene_definition = {''};
        %         jclips(lb2_index).subj_system = {''};
        %         jclips(lb2_index).inlsa_mos_2003 = NaN;
        jclips(lb2_index).scale.horizontal = 1000;
        jclips(lb2_index).scale.vertical   = 1000;

        [ds] = default_sroi( jclips(lb2_index).image_size);
        jclips(lb2_index).cvr.top = ds.top;
        jclips(lb2_index).cvr.bottom = ds.bottom;
        jclips(lb2_index).cvr.left = ds.left;
        jclips(lb2_index).cvr.right = ds.right;

        if fid ~= -1; fclose(fid); end
        lb2_index = lb2_index + 1;
    end

    waitbar(2/4); pause(0.25);

    if ~exist ('jclips_ft')
        % run function fix_temporal() with option 'FirstFrame' here
        [jclips_ft] = fix_temporal(jtests, jclips, 'FirstFrame');
        save (status_file, 'status', 'model_selected', 'jclips', 'jclips_ft', 'jtests', 'cal_type', 'model_to_run', 'working_dir', 'ssf','model_additional_options');
    end

    waitbar(3/4); pause(0.25);

    % update jclips with the final calibration report.
    jclips = jclips_ft;

    if ~exist('model_op')
            [model_op, jclips, status] = run_model(model_to_run, jtests, jclips_ft, working_dir, temp_file, cal_type, status, model_additional_options);
            jclips_ft = jclips;
            save (status_file, 'status', 'model_selected', 'jclips', 'jclips_ft', 'jtests', 'model_op', 'cal_type', 'model_to_run', 'working_dir', 'ssf','model_additional_options');
    end

    waitbar(4/4); pause(0.25);

    % if variable 'model_op' is empty, an error occurred runing the model.
    % Warn user about this.
    

    if isempty(model_op)
        model_crash(working_dir, path_sep, clip_dir);
    
         model_op_err=1;
    else
        model_op_err=0;
        [s_rpt, d_rpt] = rpt_names(working_dir, path_sep, model_to_run);
        bvqm_pc_shortrpt(model_selected, jclips_ft, jtests, model_op, cal_type,working_dir, model_to_run,s_rpt,ssf,-1, status);
        bvqm_pc_longrpt (model_selected, jclips_ft, jtests, model_op, cal_type,working_dir, model_to_run,d_rpt,ssf,-1, status);
%        close(wb);
    end
end


%-------------------------------------------------
%------------------ MANUAL CAL -------------------
%---------------------------------------------------
if (manual_cal )
    
  %Close Calibration window, begin cal
  pause(0.25);
%  delete (handles.bvqm_pc_cal);
  
  cal_type='Manual Calibration';
    
    if exist(status_file, 'file')
        load (status_file);
        crash=1;
    else
        output = model_gui(model_to_run,jclips);
        if isempty(output)
            return
        end
        if strcmp(model_to_run, 'model_psnr')
            model_to_run = output.model_to_run;
        elseif(strcmp(model_to_run, 'model_psnr_vfd'))
            model_additional_options.which_psnr=output.model_additional_options.which_psnr;
            model_additional_options.t_uncert=output.model_additional_options.t_uncert;
        elseif(strcmp(model_to_run, 'model_vqm_vfd'))
            model_additional_options.t_uncert=output.model_additional_options.t_uncert ;
            jtests.viewing_distance = output.jtests.viewing_distance;
        end
    end

    wb=waitbar(0,'BVQM is running.  Please wait...');
    pause(0.25); delete (handles.bvqm_pc_cal);
    waitbar(1/3); pause(0.25);
  
    if ~exist('jclips_mc'),
        % set a default alignment -- first frames align
        [jclips] = fix_temporal(jtests, jclips, 'FirstFrame');
        
        % save current jclips to file.  If can't, then manual calibration
        % is unavailable.
        try
            cal_file = strcat(working_dir, path_sep, 'LoadManualCalibration');
            warning off;
            [error_status] = bvqm_pc_calexport(jclips, jtests, cal_file);
        catch
            error_status = 0;
        end
        warning on;
        if ~error_status,
            uiwait( msgbox( 'Could not save to file for manual calibration. BVQM program will exit.','Fatal Error','error'));
            excel_crash(working_dir, path_sep, clip_dir);
            return;
        end
        
        uiwait( msgbox(sprintf('Edit file "%s_sheet3.csv", "%s_sheet2.csv", and "%s_sheet1.csv" with manual calibration values, and save that file.  Press "OK" when done.  Do not change location of scenes, hrcs, or columns. ', cal_file, cal_file, cal_file), ...
            'Stop! Edit Excel File', 'warn'));

        [jclips_mc] = bvqm_pc_calimport(jclips,jtests,cal_file, 1);

        save (status_file, 'status', 'model_selected', 'jclips', 'jclips_mc', 'jtests', 'cal_type', 'model_to_run', 'working_dir', 'ssf','model_additional_options');
    end

    % fix temporal registration, reframing only -- don't expect people to
    % get this right!
    if ~exist('jclips_ft')
        [jclips_ft] = fix_temporal(jtests, jclips_mc, 'Spatial');
        % errors_manual_cal.fix_temporal=status_ft;
        save (status_file, 'status', 'model_selected', 'jclips', 'jclips_mc', 'jclips_ft', 'jtests', 'cal_type', 'model_to_run', 'working_dir', 'ssf','model_additional_options');
    end

    waitbar(2/3); pause(0.25);

    % update jclips with the final calibration report.
    jclips = jclips_ft;
    
    if ~exist('model_op')
        [model_op, jclips, status] = run_model(model_to_run, jtests, jclips_ft, working_dir, temp_file, cal_type, status, model_additional_options);
        save (status_file, 'status', 'model_selected', 'jclips','jclips_mc', 'jclips_ft', 'jtests', 'model_op', 'cal_type', 'model_to_run', 'working_dir', 'ssf','model_additional_options');
    end

    waitbar(3/3); pause(0.25);
    
    % print out scaling values on reports.
    ssf = 1;
    
    % if variable 'model_op' is empty, an error occurred runing the model.
    % Warn user about this.
    if isempty(model_op)
        model_crash(working_dir, path_sep, clip_dir)
         model_op_err=1;
%         warnstr=sprintf('A fatal error occurred while running the model. \n You may have run out of memory.');
%         h=warndlg(warnstr, 'Error Occurred');
%         waitfor(h);
    else
        model_op_err=0;
        [s_rpt, d_rpt] = rpt_names(working_dir, path_sep, model_to_run);
        bvqm_pc_shortrpt(model_selected, jclips_ft, jtests, model_op, cal_type,working_dir, model_to_run,s_rpt,ssf, -1,status);
        bvqm_pc_longrpt (model_selected, jclips_ft, jtests, model_op, cal_type,working_dir, model_to_run,d_rpt,ssf, -1, status, lgo_status);
%        close(wb);
    end

end

%---------------------------------------------------
%----------------- FULL REFERENCE CAL --------------
%---------------------------------------------------
if (full_ref_cal)
    cal_type='Full Reference Calibration using ITU-T J.144/ITU-R BT.1683/ANSI T1.801.03.2003 routines';

    if exist(status_file, 'file')
        load (status_file);
        crash=1;
    end

    % query for spatial scaling option? NO - only for RR cal:
    qss=0;
    
    % get the value for uncertainty. must be >= 0
    if ~exist('uncert')
        [output]=get_uncert(qss);
        % logic to handle if 'x' was pressed
        if isempty(output)
            return
        else
        uncert=output.uncert;
        % ssf=output.ss ;  % ssf not set for full ref cal...
        end
        output = model_gui(model_to_run,jclips);
        if isempty(output)
            return
        end
        if strcmp(model_to_run, 'model_psnr')
            model_to_run = output.model_to_run;
        elseif(strcmp(model_to_run, 'model_psnr_vfd'))
            model_additional_options.which_psnr=output.model_additional_options.which_psnr;
            model_additional_options.t_uncert=output.model_additional_options.t_uncert;
        elseif(strcmp(model_to_run, 'model_vqm_vfd'))
            model_additional_options.t_uncert=output.model_additional_options.t_uncert ;
            jtests.viewing_distance = output.jtests.viewing_distance;
        end
    end

    %Close Calibration window, begin cal
    pause(0.25);
%    delete (handles.bvqm_pc_cal);
    
    wb=waitbar(0,'BVQM is running.  Please wait...');
    pause(0.25); delete (handles.bvqm_pc_cal);
    waitbar(1/6); pause(0.25);
    
    if ~exist('jclips_sr')
        % no status output available...
        [jclips_sr, status.spatial.failed, jclips_sr_unfilt] = spatial_registration(jtests, jclips, 'quiet');
        save (status_file, 'status', 'model_selected','jclips','jclips_sr','jtests','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_sr_unfilt','model_additional_options');
    end

    waitbar(2/6); pause(0.25);

    if ~exist('jclips_vr')
        [jclips_vr, status.valid, jclips_vr_unfilt] = valid_region(jtests, jclips_sr, 'Frequency', 0.5, 'quiet');
        errors_frame_based.valid_region=status.valid;
        if status.valid, % error, continue without valid
            jclips_vr = jclips_sr;
        end
        save (status_file, 'status', 'model_selected','jclips','jclips_sr','jclips_vr','jtests','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_sr_unfilt', 'jclips_vr_unfilt','model_additional_options');
    end

    waitbar(3/6); pause(0.25);

    if ~exist('jclips_lgo')
        % no status output avail.
        [jclips_lgo, status.gain_offset, jclips_lgo_unfilt] = luminance_gain_offset(jtests, jclips_vr, 'Frequency', 0.5 , 'Uncertainty', uncert, 'quiet');
        if status.gain_offset,
            jclips_lgo = jclips_vr;
        end
        % Run check on output from luminance_gain_offset: (Feb. 2007)
        [jclips_lgo, lgo_status]=check_gain_offset(jclips_lgo);
        save (status_file, 'status', 'model_selected','jclips','jclips_sr','jclips_vr','jclips_lgo','jtests','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_sr_unfilt', 'jclips_lgo_unfilt', 'jclips_vr_unfilt','model_additional_options');
    end

    waitbar(4/6); pause(0.25);

    if ~exist('jclips_tr')
        % no status output available
       [jclips_tr,status.temporal] = temporal_registration(jtests, jclips_lgo, 'Uncertainty', uncert, 'quiet');
       if status.temporal.error,
           jclips_tr = jclips_lgo;
       end
       % Per 12 Oct 06 e-mail:
       [jclips_tr] = fix_temporal(jtests, jclips_tr, 'StartPoint'); 
        save (status_file, 'status', 'model_selected','jclips','jclips_sr','jclips_vr','jclips_lgo','jclips_tr','jtests','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_sr_unfilt', 'jclips_lgo_unfilt', 'jclips_vr_unfilt','model_additional_options');
    end

    waitbar(5/6); pause(0.25);

    % update jclips with the final calibration report.
    jclips = jclips_tr;
    
    if ~exist('model_op')
        [model_op, jclips, status] = run_model(model_to_run, jtests, jclips_tr, working_dir, temp_file, cal_type, status, model_additional_options);
        save (status_file, 'status', 'model_selected','jclips','jclips_sr','jclips_vr','jclips_lgo','jclips_tr','jtests','model_op','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_sr_unfilt', 'jclips_lgo_unfilt', 'jclips_vr_unfilt','model_additional_options');
    end

    waitbar(6/6); pause(0.25);

    % if variable 'model_op' is empty, an error occurred runing the model.
    % Warn user about this.
    if isempty(model_op)
        model_crash(working_dir, path_sep, clip_dir);
         model_op_err=1;
%         warnstr=sprintf('A fatal error occurred while running the model. \n You may have run out of memory.');
%         h=warndlg(warnstr, 'Error Occurred');
%         waitfor(h);
    else
        model_op_err=0;
        [s_rpt, d_rpt] = rpt_names(working_dir, path_sep, model_to_run);
        bvqm_pc_shortrpt(model_selected, jclips_tr, jtests, model_op, cal_type,working_dir, model_to_run,s_rpt,ssf,uncert,status, lgo_status, jclips_vr_unfilt, jclips_sr_unfilt, jclips_lgo_unfilt);
        bvqm_pc_longrpt (model_selected, jclips_tr, jtests, model_op, cal_type,working_dir, model_to_run,d_rpt,ssf,uncert, status, lgo_status, jclips_vr_unfilt, jclips_sr_unfilt, jclips_lgo_unfilt);
%        close(wb);
    end
end


%-------------------------------------------------
%--------------- REDUCED REFERENCE CAL -----------
%-------------------------------------------------
if (reduced_ref_cal)
    cal_type='Reduced Reference Calibration Version 1';
%    jclips_orig=jclips;

    if exist(status_file, 'file')
        load (status_file);
        crash=1;
    end

    % Query for spatial scaling option? yes!
    qss=1;
    
    % get the value for uncertainty from the user...
    if ~exist('uncert')  % don't query if recovering from a crash - may already be assigned
        [output]=get_uncert(qss);
        % logic to handle if 'x' was pressed
        if isempty(output)
            return
        else
            uncert=output.uncert;
            ssf=output.ss ;
        end
        output = model_gui(model_to_run,jclips);
        if isempty(output)
            return
        end
        if strcmp(model_to_run, 'model_psnr')
            model_to_run = output.model_to_run;
        elseif(strcmp(model_to_run, 'model_psnr_vfd'))
            model_additional_options.which_psnr=output.model_additional_options.which_psnr;
            model_additional_options.t_uncert=output.model_additional_options.t_uncert;
        elseif(strcmp(model_to_run, 'model_vqm_vfd'))
            model_additional_options.t_uncert=output.model_additional_options.t_uncert ;
            jtests.viewing_distance = output.jtests.viewing_distance;
        end
    end

    %Close Calibration window, begin cal
    pause(0.25);
%    delete (handles.bvqm_pc_cal);
    
    wb = waitbar(0,'BVQM is running.  Please wait...');
    pause(0.25); delete (handles.bvqm_pc_cal);
    waitbar(1/8); pause(0.25);

    % run fix_temporal - just in case frames aren't aligned:
    if ~exist('jclips_ft')
        [jclips_ft] = fix_temporal(jtests, jclips, 'FirstFrame');
        save (status_file, 'status', 'model_selected','jclips','jclips_ft','jtests','cal_type','model_to_run','working_dir','ssf','uncert','model_additional_options');
    end
    
    waitbar(2/8); pause(0.25);
    
    if ~exist('jclips_tr_lbw')
        % log files available...
%        [jclips_tr_lbw, status.rr_cal.tr] = temporal_registration_lowbw(jtests, jclips_ft, working_dir, 'Uncertainty', uncert, 'quiet');
        [jclips_tr_lbw, status.temporal] = temporal_registration_lowbw(jtests, jclips_ft, working_dir, 'Uncertainty', uncert, 'field_select', 'quiet');
       if status.temporal.error,
           jclips_tr_lbw = jclips_ft;
       end
       % Per 12 Oct 06 e-mail:
       [jclips_tr_lbw] = fix_temporal(jtests, jclips_tr_lbw, 'StartPoint', 'HalfOkay'); 
        save (status_file, 'status', 'model_selected','jclips','jclips_ft','jclips_tr_lbw','jtests','cal_type','model_to_run','working_dir','ssf','uncert','model_additional_options');
    end

    waitbar(3/8); pause(0.25);

    if ~exist('jclips_srl')
        if ssf==0
            [jclips_srl, status.spatial, jclips_sr_unfilt] = spatial_registration_lowbw(jtests, jclips_tr_lbw, 'HRCFilter', 'NoScaling','quiet');
        else  %ssf=1
            [jclips_srl, status.spatial, jclips_sr_unfilt] = spatial_registration_lowbw(jtests, jclips_tr_lbw, 'HRCFilter', 'quiet');
        end
        if status.spatial.error,
            jclips_srl = jclips_tr_lbw;
        end
        save (status_file, 'status', 'model_selected','jclips','jclips_ft','jclips_tr_lbw','jclips_srl','jtests','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_sr_unfilt','model_additional_options');
    end

    waitbar(4/8); pause(0.25);

    if ~exist('jclips_vrl')
        [jclips_vrl, status.valid, jclips_vr_unfilt] = valid_region_lowbw(jtests, jclips_srl, 'Frequency', 1.0, 'quiet');
        if status.valid, % error, continue without valid
            jclips_vrl = jclips_srl;
        end
        save (status_file, 'status', 'model_selected','jclips','jclips_ft','jclips_tr_lbw','jclips_srl','jclips_vrl','jtests','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_sr_unfilt' , 'jclips_vr_unfilt','model_additional_options');
    end

    waitbar(5/8) ; pause(0.25);
    
    %New step - 17 May 06:
    if ~exist('jclips_lgo_lbw')
%       [jclips_lgo_lbw, status.rr_cal.lgo_lbw]=luminance_gain_offset_lowbw(jtests, jclips_vrl);
        [jclips_lgo_lbw, status.gain_offset, jclips_lgo_unfilt]=luminance_gain_offset_lowbw(jtests, jclips_vrl, 'HRCFilter');
         if status.gain_offset,
            jclips_lgo_lbw = jclips_vrl;
        end
       % Run check on output from luminance_gain_offset: (Feb. 2007)
        [jclips_lgo_lbw, lgo_status]=check_gain_offset(jclips_lgo_lbw);
        save (status_file, 'status', 'model_selected','jclips','jclips_ft','jclips_tr_lbw','jclips_srl','jclips_vrl','jclips_lgo_lbw','jtests','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_sr_unfilt', 'jclips_lgo_unfilt', 'jclips_vr_unfilt','model_additional_options');
    end
    
    waitbar(6/8) ; pause(0.25);
    
    if ~exist('jclips_tr_lbw2')
        % Delete feature dirs before runing temp_reg_lowbw again..
        try rmdir (fullfile(working_dir,'feature_*'), 's') ; catch; err=1; end
        [jclips_tr_lbw2, status.temporal2] = temporal_registration_lowbw(jtests, jclips_lgo_lbw, working_dir, 'Uncertainty', uncert);
       if status.temporal2.error,
           jclips_tr_lbw2 = jclips_lgo_lbw;
       end

        % Per 12 Oct 06 e-mail:
       [jclips_tr_lbw2] = fix_temporal(jtests, jclips_tr_lbw2, 'StartPoint'); 
        save (status_file, 'status', 'model_selected','jclips','jclips_ft','jclips_tr_lbw','jclips_srl','jclips_vrl','jclips_lgo_lbw','jclips_tr_lbw2','jtests','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_sr_unfilt', 'jclips_lgo_unfilt', 'jclips_vr_unfilt','model_additional_options');
    end

    waitbar(7/8); pause(0.25);

    % update jclips with the final calibration report.
    jclips = jclips_tr_lbw2;

    if ~exist('model_op')
            [model_op, jclips, status] = run_model(model_to_run, jtests, jclips_tr_lbw2, working_dir, temp_file, cal_type, status, model_additional_options);
            save (status_file, 'status', 'model_selected','jclips','jclips_ft', 'jclips_tr_lbw','jclips_srl','jclips_vrl','jclips_lgo_lbw','jclips_tr_lbw2','jtests','model_op','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_sr_unfilt', 'jclips_lgo_unfilt', 'jclips_vr_unfilt','model_additional_options');
    end

    waitbar(8/8); pause(0.25);

    % if variable 'model_op' is empty, an error occurred runing the model.
    % Warn user about this.
    if isempty(model_op)
        model_crash(working_dir, path_sep, clip_dir);
%         
         model_op_err=1;
%         warnstr=sprintf('A fatal error occurred while running the model. \n You may have run out of memory.');
%         h=warndlg(warnstr, 'Error Occurred');
%         waitfor(h);
    else
        model_op_err=0;
        [s_rpt, d_rpt] = rpt_names(working_dir, path_sep, model_to_run);
        bvqm_pc_shortrpt(model_selected, jclips_tr_lbw2,jtests,model_op,cal_type,working_dir,model_to_run,s_rpt,ssf,uncert,status,lgo_status, jclips_vr_unfilt, jclips_sr_unfilt, jclips_lgo_unfilt);
        bvqm_pc_longrpt (model_selected, jclips_tr_lbw2,jtests,model_op,cal_type,working_dir,model_to_run,d_rpt,ssf,uncert,status,lgo_status, jclips_vr_unfilt, jclips_sr_unfilt, jclips_lgo_unfilt);
%        close(wb);
    end  
end

%-------------------------------------------------
%--------------- REDUCED REFERENCE CAL -----------
%-------------------------------------------------
if (reduced_ref2_cal)
    cal_type='Reduced Reference Calibration Version 2';
%    jclips_orig=jclips;

    if exist(status_file, 'file')
        load (status_file);
        crash=1;
    end

    % Query for spatial scaling option? yes!
    qss=1;
    
    % get the value for uncertainty from the user...
    if ~exist('uncert')  % don't query if recovering from a crash - may already be assigned
        [output]=get_uncert(qss);
        % logic to handle if 'x' was pressed
        if isempty(output)
            return
        else
            uncert=output.uncert;
            ssf=output.ss ;
        end
        output = model_gui(model_to_run,jclips);
        if isempty(output)
            return
        end
        if strcmp(model_to_run, 'model_psnr')
            model_to_run = output.model_to_run;
        elseif(strcmp(model_to_run, 'model_psnr_vfd'))
            model_additional_options.which_psnr=output.model_additional_options.which_psnr;
            model_additional_options.t_uncert=output.model_additional_options.t_uncert;
        elseif(strcmp(model_to_run, 'model_vqm_vfd'))
            model_additional_options.t_uncert=output.model_additional_options.t_uncert ;
            jtests.viewing_distance = output.jtests.viewing_distance;
        end
    end

    %Close Calibration window, begin cal
    pause(0.25);
%    delete (handles.bvqm_pc_cal);
    
    wb = waitbar(0,'BVQM is running.  Please wait...');
    pause(0.25); delete (handles.bvqm_pc_cal);
    waitbar(1/8); pause(0.25);

    % run fix_temporal - just in case frames aren't aligned:
    if ~exist('jclips_ft')
        [jclips_ft] = fix_temporal(jtests, jclips, 'FirstFrame');
        save (status_file, 'status', 'model_selected','jclips','jclips_ft','jtests','cal_type','model_to_run','working_dir','ssf','uncert','model_additional_options');
    end
    
    waitbar(2/8); pause(0.25);
    
    if ~exist('jclips_tr_lbw')
        % log files available...
%        [jclips_tr_lbw, status.rr_cal.tr] = temporal_registration_lowbw(jtests, jclips_ft, working_dir, 'Uncertainty', uncert, 'quiet');
        [jclips_tr_lbw, status.temporal] = temporal_registration_lowbw(jtests, jclips_ft, working_dir, 'Uncertainty', uncert, 'field_select', 'quiet');
       if status.temporal.error,
           jclips_tr_lbw = jclips_ft;
       end
       % Per 12 Oct 06 e-mail:
       [jclips_tr_lbw] = fix_temporal(jtests, jclips_tr_lbw, 'StartPoint', 'HalfOkay'); 
        save (status_file, 'status', 'model_selected','jclips','jclips_ft','jclips_tr_lbw','jtests','cal_type','model_to_run','working_dir','ssf','uncert','model_additional_options');
    end

    waitbar(3/8); pause(0.25);

    if ~exist('jclips_srl')
        if ssf==0
            [jclips_srl, status.spatial, jclips_sr_unfilt] = spatial_registration_lowbw(jtests, jclips_tr_lbw, 'HRCFilter', 'NoScaling','quiet');
        else  %ssf=1
            [jclips_srl, status.spatial, jclips_sr_unfilt] = spatial_registration_lowbw(jtests, jclips_tr_lbw, 'HRCFilter', 'quiet');
        end
        if status.spatial.error,
            jclips_srl = jclips_tr_lbw;
        end
        save (status_file, 'status', 'model_selected','jclips','jclips_ft','jclips_tr_lbw','jclips_srl','jtests','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_sr_unfilt','model_additional_options');
    end

    waitbar(4/8); pause(0.25);

    if ~exist('jclips_vrl')
        [jclips_vrl, status.valid, jclips_vr_unfilt] = valid_region_lowbw(jtests, jclips_srl, 'Frequency', 1.0, 'quiet');
        if status.valid, % error, continue without valid
            jclips_vrl = jclips_srl;
        end
        save (status_file, 'status', 'model_selected','jclips','jclips_ft','jclips_tr_lbw','jclips_srl','jclips_vrl','jtests','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_sr_unfilt' , 'jclips_vr_unfilt','model_additional_options');
    end

    waitbar(5/8) ; pause(0.25);
    
    %New step - 17 May 06:
    if ~exist('jclips_lgo_lbw')
%       [jclips_lgo_lbw, status.rr_cal.lgo_lbw]=luminance_gain_offset_lowbw(jtests, jclips_vrl);
        [jclips_lgo_lbw, status.gain_offset, jclips_lgo_unfilt]=gain_offset_lowbw(jtests, jclips_vrl, 'HRCFilter');
        if status.gain_offset,
            jclips_lgo_lbw = jclips_vrl;
        end
        % Run check on output from luminance_gain_offset: (Feb. 2007)
        [jclips_lgo_lbw, lgo_status]=check_gain_offset(jclips_lgo_lbw);
        save (status_file, 'status', 'model_selected','jclips','jclips_ft','jclips_tr_lbw','jclips_srl','jclips_vrl','jclips_lgo_lbw','jtests','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_sr_unfilt', 'jclips_lgo_unfilt', 'jclips_vr_unfilt','model_additional_options');
    end
    
    waitbar(6/8) ; pause(0.25);
    
    if ~exist('jclips_tr_lbw2')
        % Delete feature dirs before runing temp_reg_lowbw again..
        try rmdir (fullfile(working_dir,'feature_*'), 's') ; catch; err=1; end
        [jclips_tr_lbw2, status.temporal2] = temporal_registration_lowbw(jtests, jclips_lgo_lbw, working_dir, 'Uncertainty', uncert);
       if status.temporal.error,
           jclips_tr_lbw2 = jclips_lgo_lbw;
       end

        % Per 12 Oct 06 e-mail:
       [jclips_tr_lbw2] = fix_temporal(jtests, jclips_tr_lbw2, 'StartPoint'); 
        save (status_file, 'status', 'model_selected','jclips','jclips_ft','jclips_tr_lbw','jclips_srl','jclips_vrl','jclips_lgo_lbw','jclips_tr_lbw2','jtests','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_sr_unfilt', 'jclips_lgo_unfilt', 'jclips_vr_unfilt','model_additional_options');
    end

    waitbar(7/8); pause(0.25);

    % update jclips with the final calibration report.
    jclips = jclips_tr_lbw2;

    if ~exist('model_op')
            [model_op, jclips, status] = run_model(model_to_run, jtests, jclips_tr_lbw2, working_dir, temp_file, cal_type, status, model_additional_options);
            save (status_file, 'status', 'model_selected','jclips','jclips_ft', 'jclips_tr_lbw','jclips_srl','jclips_vrl','jclips_lgo_lbw','jclips_tr_lbw2','jtests','model_op','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_sr_unfilt', 'jclips_lgo_unfilt', 'jclips_vr_unfilt','model_additional_options');
    end

    waitbar(8/8); pause(0.25);

    % if variable 'model_op' is empty, an error occurred runing the model.
    % Warn user about this.
    if isempty(model_op)
        model_crash(working_dir, path_sep, clip_dir);
%         
         model_op_err=1;
%         warnstr=sprintf('A fatal error occurred while running the model. \n You may have run out of memory.');
%         h=warndlg(warnstr, 'Error Occurred');
%         waitfor(h);
    else
        model_op_err=0;
        [s_rpt, d_rpt] = rpt_names(working_dir, path_sep, model_to_run);
        bvqm_pc_shortrpt(model_selected, jclips_tr_lbw2,jtests,model_op,cal_type,working_dir,model_to_run,s_rpt,ssf,uncert,status,lgo_status, jclips_vr_unfilt, jclips_sr_unfilt, jclips_lgo_unfilt);
        bvqm_pc_longrpt (model_selected, jclips_tr_lbw2,jtests,model_op,cal_type,working_dir,model_to_run,d_rpt,ssf,uncert,status,lgo_status, jclips_vr_unfilt, jclips_sr_unfilt, jclips_lgo_unfilt);
%        close(wb);
    end  
end

if tr_vr,
    if ~isfield(status,'is_rr'),
        button = questdlg('Which algorithm would you like to use for Temporal Registration and Valid Video Region Estimation?', ...
            'Choose Algoritm', 'Full Reference', 'Reduced Reference', 'Reduced Reference');
        if strcmpi(button, 'Reduced Reference'),
            status.is_rr = 1;
            status.algorithm = 'Reduced Reference';
        else
            status.is_rr = 0;
            status.algorithm = 'Full Reference';
        end
        cal_type='Temporal Registration and Valid Region Only';
        output = model_gui(model_to_run,jclips);
        if isempty(output)
            return
        end
        if strcmp(model_to_run, 'model_psnr')
            model_to_run = output.model_to_run;
        elseif(strcmp(model_to_run, 'model_psnr_vfd'))
            model_additional_options.which_psnr=output.model_additional_options.which_psnr;
            model_additional_options.t_uncert=output.model_additional_options.t_uncert;
        elseif(strcmp(model_to_run, 'model_vqm_vfd'))
            model_additional_options.t_uncert=output.model_additional_options.t_uncert ;
            jtests.viewing_distance = output.jtests.viewing_distance;
        end
        save (status_file, 'status', 'model_selected','jclips','jtests','cal_type','model_to_run','working_dir','ssf','model_additional_options');
    end
    if status.is_rr == 0,
        %--------------------------------------------------------------------------
        %--------------- TEMP REGISTRATION, VALID REGION - FULL REF CAL -----------
        %--------------------------------------------------------------------------
        %      tr_vr = get(handles.radiobutton12, 'Value');
        %      mc_then_tr = get(handles.radiobutton13, 'Value');

        if exist(status_file, 'file')
            load (status_file);
            crash=1;
        end

        % query for spatial scaling option? NO 
        qss=0;

        % get the value for uncertainty. must be >= 0
        if ~exist('uncert')
            [output]=get_uncert(qss);
           % logic to handle if 'x' was pressed
            if isempty(output)
                return
            else
                uncert=output.uncert;
                % ssf=output.ss ;
            end
        end

        %Close Calibration window, begin cal
        pause(0.25);
    %    delete (handles.bvqm_pc_cal);

        wb=waitbar(0,'BVQM is running.  Please wait...');
        pause(0.25); delete (handles.bvqm_pc_cal);
        waitbar(1/3); pause(0.25);


        if ~exist('jclips_vr')
            [jclips_vr, status.valid, jclips_vr_unfilt] = valid_region(jtests, jclips, 'Frequency', 0.5, 'quiet');
            errors_frame_based.valid_region=status.valid;
            if status.valid, % error, continue without valid
                jclips_vr = jclips;
            end
            save (status_file, 'status', 'model_selected','jclips','jclips_vr','jtests','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_vr_unfilt','model_additional_options');
        end

        waitbar(2/3); pause(0.25);

        if ~exist('jclips_tr')
           [jclips_tr,status.temporal] = temporal_registration(jtests, jclips_vr, 'Uncertainty', uncert, 'quiet');
           if status.temporal.error,
               jclips_tr = jclips_vr;
           end
           % Per 12 Oct 06 e-mail:
           [jclips_tr] = fix_temporal(jtests, jclips_tr, 'StartPoint'); 
            save (status_file, 'status', 'model_selected','jclips','jclips_vr', 'jclips_tr','jtests','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_vr_unfilt','model_additional_options');
        end

        % update jclips with the final calibration report.
        jclips = jclips_tr;

        if ~exist('model_op')
            [model_op, jclips, status] = run_model(model_to_run, jtests, jclips_tr, working_dir, temp_file, cal_type, status, model_additional_options);
            save (status_file, 'status', 'model_selected','jclips','jclips_vr','jclips_tr','jtests','model_op','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_vr_unfilt','model_additional_options');
        end

        waitbar(3/3); pause(0.25);

        % if variable 'model_op' is empty, an error occurred runing the model.
        % Warn user about this.
        if isempty(model_op)
            model_crash(working_dir, path_sep, clip_dir);
             model_op_err=1;
    %         warnstr=sprintf('A fatal error occurred while running the model. \n You may have run out of memory.');
    %         h=warndlg(warnstr, 'Error Occurred');
    %         waitfor(h);
        else
            model_op_err=0;
            [s_rpt, d_rpt] = rpt_names(working_dir, path_sep, model_to_run);
            bvqm_pc_shortrpt(model_selected, jclips_tr, jtests, model_op, cal_type,working_dir, model_to_run,s_rpt,ssf, uncert,status, [], jclips_vr_unfilt);
            bvqm_pc_longrpt (model_selected, jclips_tr, jtests, model_op, cal_type,working_dir, model_to_run,d_rpt,ssf, uncert,status, [], jclips_vr_unfilt);
        end


    else

        %--------------------------------------------------------------------------
        %--------------- TEMP REGISTRATION, VALID REGION - REDUCED REF CAL -----------
        %--------------------------------------------------------------------------

        if exist(status_file, 'file')
            load (status_file);
            crash=1;
        end

        % query for spatial scaling option? NO 
        qss=0;

        % get the value for uncertainty. must be >= 0
        if ~exist('uncert')
            [output]=get_uncert(qss);
            % logic to handle if 'x' was pressed
            if isempty(output)
                return
            else
                uncert=output.uncert;
                % ssf=output.ss ;  % not for this cal
            end
        end

        %Close Calibration window, begin cal
        pause(0.25);

        wb=waitbar(0,'BVQM is running.  Please wait...');
        pause(0.25); delete (handles.bvqm_pc_cal);
        waitbar(1/3); pause(0.25);

        % fix any temporal problems first: 
        % no need to save results this step b/c it is very fast.
        [jclips]=fix_temporal(jtests, jclips, 'FirstFrame');

        if ~exist('jclips_vr')
            [jclips_vr, status.valid, jclips_vr_unfilt] = valid_region_lowbw(jtests, jclips, 'Frequency', 0.5, 'quiet');
            errors_frame_based.valid_region=status.valid;
            if status.valid, % error, continue without valid
                jclips_vr = jclips;
            end
            save (status_file, 'status', 'model_selected','jclips','jclips_vr','jtests','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_vr_unfilt','model_additional_options');
        end

        waitbar(2/3); pause(0.25);

        if ~exist('jclips_tr_lbw')
            % no status output available
           [jclips_tr_lbw,status.temporal] = temporal_registration_lowbw(jtests, jclips_vr, working_dir, 'frame_select', 'unaligned', 'Uncertainty', uncert, 'quiet');
           if status.temporal.error,
               jclips_tr_lbw = jclips_vr;
           end
           % Per 12 Oct 06 e-mail:
           [jclips_tr_lbw] = fix_temporal(jtests, jclips_tr_lbw, 'StartPoint'); 
            save (status_file, 'status', 'model_selected','jclips','jclips_vr', 'jclips_tr_lbw','jtests','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_vr_unfilt','model_additional_options');
        end

        % update jclips with the final calibration report.
        jclips = jclips_tr_lbw;

        if ~exist('model_op')
            [model_op, jclips, status] = run_model(model_to_run, jtests, jclips_tr_lbw, working_dir, temp_file, cal_type, status, model_additional_options);
            save (status_file, 'status', 'model_selected','jclips','jclips_vr','jclips_tr_lbw','jtests','model_op','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_vr_unfilt','model_additional_options');
        end

        waitbar(3/3); pause(0.25);

        % if variable 'model_op' is empty, an error occurred runing the model.
        % Warn user about this.
        if isempty(model_op)
            model_crash(working_dir, path_sep, clip_dir);
             model_op_err=1;
    %         warnstr=sprintf('A fatal error occurred while running the model. \n You may have run out of memory.');
    %         h=warndlg(warnstr, 'Error Occurred');
    %         waitfor(h);
        else
            model_op_err=0;
            [s_rpt, d_rpt] = rpt_names(working_dir, path_sep, model_to_run);
            bvqm_pc_shortrpt(model_selected, jclips_tr_lbw, jtests, model_op, cal_type,working_dir, model_to_run,s_rpt,ssf,uncert,status, [], jclips_vr_unfilt);
            bvqm_pc_longrpt (model_selected, jclips_tr_lbw, jtests, model_op, cal_type,working_dir, model_to_run,d_rpt,ssf, uncert,status, [], jclips_vr_unfilt);
    %        close(wb);
        end
    end
end

if mc_then_tr,
    if ~isfield(status,'is_rr'),
        button = questdlg('Which algorithm would you like to use for Temporal Registration and Valid Video Region Estimation?', ...
            'Choose Algorithm','Full Reference', 'Reduced Reference', 'Reduced Reference');
        if strcmpi(button, 'Reduced Reference'),
            status.is_rr = 1;
            status.algorithm = 'Reduced Reference';
        else
            status.is_rr = 0;
            status.algorithm = 'Full Reference';
        end
        cal_type='Temporal Registration after Manual Calibration';
        output = model_gui(model_to_run,jclips);
        if isempty(output)
            return
        end
        if strcmp(model_to_run, 'model_psnr')
            model_to_run = output.model_to_run;
        elseif(strcmp(model_to_run, 'model_psnr_vfd'))
            model_additional_options.which_psnr=output.model_additional_options.which_psnr;
            model_additional_options.t_uncert=output.model_additional_options.t_uncert;
        elseif(strcmp(model_to_run, 'model_vqm_vfd'))
            model_additional_options.t_uncert=output.model_additional_options.t_uncert ;
            jtests.viewing_distance = output.jtests.viewing_distance;
        end
        save (status_file, 'status', 'model_selected','jclips','jtests','cal_type','model_to_run','working_dir','ssf','model_additional_options');
    end
    if status.is_rr == 0,
        %--------------------------------------------------------------------------
        %--------------- TEMP REGISTRATION, MANUAL CALIBRATION - FULL REF CAL -----------
        %--------------------------------------------------------------------------

        if exist(status_file, 'file')
            load (status_file);
            crash=1;
        end

        % query for spatial scaling option? NO 
        qss=0;

        % get the value for uncertainty. must be >= 0
        if ~exist('uncert')
            [output]=get_uncert(qss);
           % logic to handle if 'x' was pressed
            if isempty(output)
                return
            else
                uncert=output.uncert;
                % ssf=output.ss ;
            end
        end

        %Close Calibration window, begin cal
        pause(0.25);
    %    delete (handles.bvqm_pc_cal);

        wb=waitbar(0,'BVQM is running.  Please wait...');
        pause(0.25); delete (handles.bvqm_pc_cal);
        waitbar(1/3); pause(0.25);

        if ~exist('jclips_mc'),
            % set a default alignment -- first frames align
            [jclips] = fix_temporal(jtests, jclips, 'FirstFrame');

            % save current jclips to file.  If can't, then manual calibration
            % is unavailable.
            try
                cal_file = strcat(working_dir, path_sep, 'LoadManualCalibration');
                warning off;
                [error_status] = bvqm_pc_calexport(jclips, jtests, cal_file);
            catch
                error_status = 0;
            end
            warning on;
            if ~error_status,
                uiwait( msgbox( 'Could not save to file for manual calibration.  BVQM program will exit.','Fatal Error','error'));
                excel_crash(working_dir, path_sep, clip_dir);
                return;
            end
            uiwait( msgbox(sprintf('Edit file "%s_sheet3.csv", "%s_sheet2.csv", and "%s_sheet1.csv" with manual calibration values, and save that file.  Press "OK" when done.  Do not change location of scenes, hrcs, or columns. ', cal_file, cal_file, cal_file), ...
                'Stop! Edit Excel File', 'warn'));

            [jclips_mc] = bvqm_pc_calimport(jclips,jtests,cal_file, 0);

            save (status_file, 'status', 'model_selected', 'jclips', 'jclips_mc', 'jtests', 'cal_type', 'model_to_run', 'working_dir', 'ssf','model_additional_options');
        end


        waitbar(2/3); pause(0.25);

        if ~exist('jclips_tr')
           [jclips_tr,status.temporal] = temporal_registration(jtests, jclips_mc, 'Uncertainty', uncert, 'quiet');
           if status.temporal.error,
               jclips_tr = jclips_mc;
           end
           % Per 12 Oct 06 e-mail:
           [jclips_tr] = fix_temporal(jtests, jclips_tr, 'StartPoint'); 
            save (status_file, 'status', 'model_selected','jclips','jclips_mc', 'jclips_tr','jtests','cal_type','model_to_run','working_dir','ssf','uncert','model_additional_options');
        end

        % update jclips with the final calibration report.
        jclips = jclips_tr;

        if ~exist('model_op')
            [model_op, jclips, status] = run_model(model_to_run, jtests, jclips_tr, working_dir, temp_file, cal_type, status, model_additional_options);
            save (status_file, 'status', 'model_selected','jclips','jclips_mc','jclips_tr','jtests','model_op','cal_type','model_to_run','working_dir','ssf','uncert','model_additional_options');
        end

        waitbar(3/3); pause(0.25);

        % if variable 'model_op' is empty, an error occurred runing the model.
        % Warn user about this.
        if isempty(model_op)
            model_crash(working_dir, path_sep, clip_dir);
             model_op_err=1;
    %         warnstr=sprintf('A fatal error occurred while running the model. \n You may have run out of memory.');
    %         h=warndlg(warnstr, 'Error Occurred');
    %         waitfor(h);
        else
            model_op_err=0;
            [s_rpt, d_rpt] = rpt_names(working_dir, path_sep, model_to_run);
            bvqm_pc_shortrpt(model_selected, jclips_tr, jtests, model_op, cal_type,working_dir, model_to_run,s_rpt,ssf, uncert,status, []);
            bvqm_pc_longrpt (model_selected, jclips_tr, jtests, model_op, cal_type,working_dir, model_to_run,d_rpt,ssf, uncert,status, []);
        end


    else

        %--------------------------------------------------------------------------
        %--------------- TEMP REGISTRATION, Manual Calibration - REDUCED REF CAL -----------
        %--------------------------------------------------------------------------

        if exist(status_file, 'file')
            load (status_file);
            crash=1;
        end

        % query for spatial scaling option? NO 
        qss=0;

        % get the value for uncertainty. must be >= 0
        if ~exist('uncert')
            [output]=get_uncert(qss);
            % logic to handle if 'x' was pressed
            if isempty(output)
                return
            else
                uncert=output.uncert;
                % ssf=output.ss ;  % not for this cal
            end
        end

        %Close Calibration window, begin cal
        pause(0.25);

        wb=waitbar(0,'BVQM is running.  Please wait...');
        pause(0.25); delete (handles.bvqm_pc_cal);
        waitbar(1/3); pause(0.25);

        if ~exist('jclips_mc'),
            % set a default alignment -- first frames align
            [jclips] = fix_temporal(jtests, jclips, 'FirstFrame');

            % save current jclips to file.  If can't, then manual calibration
            % is unavailable.
            try
                cal_file = strcat(working_dir, path_sep, 'LoadManualCalibration');
                warning off;
                [error_status] = bvqm_pc_calexport(jclips, jtests, cal_file);
            catch
                error_status = 0;
            end
            warning on;
            if ~error_status,
                uiwait( msgbox( 'Could not save to file for manual calibration. VQM program will exit.','Fatal Error','error'));
                excel_crash(working_dir, path_sep, clip_dir);
                return;
            end
            uiwait( msgbox(sprintf('Edit file "%s_sheet3.csv", "%s_sheet2.csv", and "%s_sheet1.csv" with manual calibration values, and save that file.  Press "OK" when done.  Do not change location of scenes, hrcs, or columns. ', cal_file, cal_file, cal_file), ...
                'Stop! Edit Excel File', 'warn'));

            [jclips_mc] = bvqm_pc_calimport(jclips,jtests,cal_file, 0);

            save (status_file, 'status', 'model_selected', 'jclips', 'jclips_mc', 'jtests', 'cal_type', 'model_to_run', 'working_dir', 'ssf','model_additional_options');
        end

        waitbar(2/3); pause(0.25);

        if ~exist('jclips_tr_lbw')
            % no status output available
           [jclips_tr_lbw,status.temporal] = temporal_registration_lowbw(jtests, jclips_mc, working_dir, 'frame_select', 'unaligned', 'Uncertainty', uncert, 'quiet');
           if status.temporal.error,
               jclips_tr_lbw = jclips_mc;
           end
           % Per 12 Oct 06 e-mail:
           [jclips_tr_lbw] = fix_temporal(jtests, jclips_tr_lbw, 'StartPoint'); 
            save (status_file, 'status', 'model_selected','jclips','jclips_mc', 'jclips_tr_lbw','jtests','cal_type','model_to_run','working_dir','ssf','uncert','model_additional_options');
        end

        % update jclips with the final calibration report.
        jclips = jclips_tr_lbw;

        if ~exist('model_op')
            [model_op, jclips, status] = run_model(model_to_run, jtests, jclips_tr_lbw, working_dir, temp_file, cal_type, status, model_additional_options);
            save (status_file, 'status', 'model_selected','jclips','jclips_tr_lbw','jclips_mc','jtests','model_op','cal_type','model_to_run','working_dir','ssf','uncert','model_additional_options');
        end

        waitbar(3/3); pause(0.25);

        % if variable 'model_op' is empty, an error occurred runing the model.
        % Warn user about this.
        if isempty(model_op)
            model_crash(working_dir, path_sep, clip_dir);
             model_op_err=1;
    %         warnstr=sprintf('A fatal error occurred while running the model. \n You may have run out of memory.');
    %         h=warndlg(warnstr, 'Error Occurred');
    %         waitfor(h);
        else
            model_op_err=0;
            [s_rpt, d_rpt] = rpt_names(working_dir, path_sep, model_to_run);
            bvqm_pc_shortrpt(model_selected, jclips_tr_lbw, jtests, model_op, cal_type,working_dir, model_to_run,s_rpt,ssf,uncert,status, []);
            bvqm_pc_longrpt (model_selected, jclips_tr_lbw, jtests, model_op, cal_type,working_dir, model_to_run,d_rpt,ssf, uncert,status, []);
    %        close(wb);
        end
    end
    
end

%-------------------------------------------------
%--------- Peak Signal to Noise Ratio CAL --------
%-------------------------------------------------
if(psnr)
    %Close Calibration window, begin cal
    pause(0.25);
    %  delete (handles.bvqm_pc_cal);
    
    cal_type='Peak Signal to Noise Ratio';
    
    if exist(status_file, 'file')
        load (status_file);
        crash=1;
    end
    if ~exist('spatial_x')  % don't query if recovering from a crash - may already be assigned
        [output]=bvqm_pc_psnr_cal();
        if isempty(output)
            return
        else
            tuncert_psnr=output.uncert;
            frac_samp = output.frac_samp;
            spatial_y = output.spatial_y;
            spatial_x = output.spatial_x;
        end
        output = model_gui(model_to_run,jclips);
        if isempty(output)
            return
        end
        if strcmp(model_to_run, 'model_psnr')
            model_to_run = output.model_to_run;
        elseif(strcmp(model_to_run, 'model_psnr_vfd'))
            model_additional_options.which_psnr=output.model_additional_options.which_psnr;
            model_additional_options.t_uncert=output.model_additional_options.t_uncert;
        elseif(strcmp(model_to_run, 'model_vqm_vfd'))
            model_additional_options.t_uncert=output.model_additional_options.t_uncert ;
            jtests.viewing_distance = output.jtests.viewing_distance;
        end
        save (status_file, 'status', 'model_selected','jclips','jtests','cal_type','model_to_run','working_dir','tuncert_psnr','frac_samp','spatial_x','spatial_y','model_additional_options');
    end
    
    wb=waitbar(0,'BVQM is running.  Please wait...');
    pause(0.25); delete (handles.bvqm_pc_cal);
    waitbar(1/4); pause(0.25);
    
    if ~exist('crash')
        %Either this is the first time this function is running on the
        %clipset or the user wants to "start fresh" meaning that if PSNR
        %crashed then we need to make sure it removes all of its files
        %before running again.
        %Only calibrated version of PSNR will run in BVQM...
        warning off
        try
            delete([working_dir '\' 'temp_gclips.mat']);
        catch
            %File does not exist!
        end
        try
            delete([working_dir '\' 'temp_struct_orig.mat']);
        catch
            %File does not exist!
        end
        try
            delete([working_dir '\' 'tshifts.mat']);
        catch
            %File does not exist!
        end
        try
            delete([working_dir '\' 'results_psnr.mat']);
        catch
            %File does not exist!
        end
        try
            delete([working_dir '\' 'results_psnr.csv']);
        catch
            %File does not exist!
        end
        warning on
    end
    
    %Run PSNR Search and update jclips.
    if ~exist('jclips_temp')
        jclips_temp = bvqm_psnr_search_gclips(jtests,jclips,'results_psnr.csv',working_dir,'results_psnr',...
            'fraction_sampled',frac_samp,'spatial_uncertainty',spatial_x,spatial_y,...
            'temporal_uncertainty',tuncert_psnr,'silent','calibrated');
        delete([working_dir '\' 'results_psnr.mat']);
        delete([working_dir '\' 'results_psnr.csv']);
        save (status_file, 'status', 'model_selected','jclips','jtests','cal_type','model_to_run','working_dir','tuncert_psnr','frac_samp','spatial_x','spatial_y','jclips_temp','model_additional_options');
    end
    
    waitbar(2/4); pause(0.25);
    
    %Valid Region after search is WRONG.  Needs to be recalculated.
    if ~exist('jclips_vr')
        [jclips_vr, status.valid, jclips_vr_unfilt] = valid_region(jtests, jclips_temp,'quiet');
        if status.valid, % error, continue without valid
            jclips_vr = jclips_temp;
        end
        save (status_file, 'status', 'model_selected','jclips','jtests','cal_type','model_to_run','working_dir','tuncert_psnr','frac_samp','spatial_x','spatial_y','jclips_temp','jclips_vr','jclips_vr_unfilt','model_additional_options');
    end
    
    jclips = jclips_vr;
    
    waitbar(3/4); pause(0.25);

    if ~exist('model_op')
            [model_op, jclips, status] = run_model(model_to_run, jtests, jclips, working_dir, temp_file, cal_type, status, model_additional_options);
            save (status_file, 'status', 'model_selected','jclips','jtests','model_op','cal_type','model_to_run','working_dir','tuncert_psnr','frac_samp','spatial_x','spatial_y','jclips_temp','jclips_vr','jclips_vr_unfilt','model_additional_options');
    end

    waitbar(4/4); pause(0.25);

    % if variable 'model_op' is empty, an error occurred runing the model.
    % Warn user about this.
    if isempty(model_op)
        model_crash(working_dir, path_sep, clip_dir);
%         
         model_op_err=1;
%         warnstr=sprintf('A fatal error occurred while running the model. \n You may have run out of memory.');
%         h=warndlg(warnstr, 'Error Occurred');
%         waitfor(h);
    else
        model_op_err=0;
        [s_rpt, d_rpt] = rpt_names(working_dir, path_sep, model_to_run);
        ssf = 0;
        bvqm_pc_shortrpt(model_selected,jclips,jtests,model_op,cal_type,working_dir,model_to_run,s_rpt,ssf,tuncert_psnr,status,[],[],[],[]);
        bvqm_pc_longrpt (model_selected,jclips,jtests,model_op,cal_type,working_dir,model_to_run,d_rpt,ssf,tuncert_psnr,status,[],jclips_vr_unfilt,[],[]);
    end 
end

%---------------------------------------------------
%--------- Fast PSNR Search 1 (rrcal 2)
%---------------------------------------------------
if (rrcal2_then_psnr)
    cal_type = 'ITU-T J.244 followed by ITU-T J.340';
    if exist(status_file, 'file')
        load (status_file);
        crash=1;
    end

    % Query for spatial scaling option? yes!
    qss=1;
    
    % get values from the user...
    if ~exist('uncert')
        [output]=bvqm_pc_rrcal2_then_psnr_cal();
        if isempty(output)
            return
        else
            tuncert_psnr = output.uncert;
            frac_samp = output.frac_samp;
            spatial_y = output.spatial_y;
            spatial_x = output.spatial_x;
            uncert = output.rrcal2_uncert;
            ssf = output.perform_ss;
        end
        output = model_gui(model_to_run,jclips);
        if isempty(output)
            return
        end
        if strcmp(model_to_run, 'model_psnr')
            model_to_run = output.model_to_run;
        elseif(strcmp(model_to_run, 'model_psnr_vfd'))
            model_additional_options.which_psnr=output.model_additional_options.which_psnr;
            model_additional_options.t_uncert=output.model_additional_options.t_uncert;
        elseif(strcmp(model_to_run, 'model_vqm_vfd'))
            model_additional_options.t_uncert=output.model_additional_options.t_uncert ;
            jtests.viewing_distance = output.jtests.viewing_distance;
        end
    end
    

    %Close Calibration window, begin cal
    pause(0.25);
%    delete (handles.bvqm_pc_cal);
    
    wb = waitbar(0,'BVQM is running.  Please wait...');
    pause(0.25); delete (handles.bvqm_pc_cal);
    waitbar(1/9); pause(0.25);

    % run fix_temporal - just in case frames aren't aligned:
    if ~exist('jclips_ft')
        [jclips_ft] = fix_temporal(jtests, jclips, 'FirstFrame');
        save (status_file, 'status', 'model_selected','jclips','jclips_ft','jtests','cal_type','model_to_run','working_dir','ssf','uncert','tuncert_psnr','frac_samp','spatial_y','spatial_x','model_additional_options');
    end
    
    waitbar(2/9); pause(0.25);
    
    if ~exist('jclips_tr_lbw')
        % log files available...
%        [jclips_tr_lbw, status.rr_cal.tr] = temporal_registration_lowbw(jtests, jclips_ft, working_dir, 'Uncertainty', uncert, 'quiet');
        [jclips_tr_lbw, status.temporal] = temporal_registration_lowbw(jtests, jclips_ft, working_dir, 'Uncertainty', uncert, 'field_select', 'quiet');
       if status.temporal.error,
           jclips_tr_lbw = jclips_ft;
       end
       % Per 12 Oct 06 e-mail:
       [jclips_tr_lbw] = fix_temporal(jtests, jclips_tr_lbw, 'StartPoint', 'HalfOkay'); 
        save (status_file, 'status', 'model_selected','jclips','jclips_ft','jclips_tr_lbw','jtests','cal_type','model_to_run','working_dir','ssf','uncert','tuncert_psnr','frac_samp','spatial_y','spatial_x','model_additional_options');
    end

    waitbar(3/9); pause(0.25);

    if ~exist('jclips_srl')
        if ssf==0
            [jclips_srl, status.spatial, jclips_sr_unfilt] = spatial_registration_lowbw(jtests, jclips_tr_lbw, 'HRCFilter', 'NoScaling','quiet');
        else  %ssf=1
            [jclips_srl, status.spatial, jclips_sr_unfilt] = spatial_registration_lowbw(jtests, jclips_tr_lbw, 'HRCFilter', 'quiet');
        end
        if status.spatial.error,
            jclips_srl = jclips_tr_lbw;
        end
        save (status_file, 'status', 'model_selected','jclips','jclips_ft','jclips_tr_lbw','jclips_srl','jtests','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_sr_unfilt','tuncert_psnr','frac_samp','spatial_y','spatial_x','model_additional_options');
    end

    waitbar(4/9); pause(0.25);

    if ~exist('jclips_vrl')
        [jclips_vrl, status.valid, jclips_vr_unfilt] = valid_region_lowbw(jtests, jclips_srl, 'Frequency', 1.0, 'quiet');
        if status.valid, % error, continue without valid
            jclips_vrl = jclips_srl;
        end
        save (status_file, 'status', 'model_selected','jclips','jclips_ft','jclips_tr_lbw','jclips_srl','jclips_vrl','jtests','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_sr_unfilt' , 'jclips_vr_unfilt','tuncert_psnr','frac_samp','spatial_y','spatial_x','model_additional_options');
    end
    
    waitbar(5/9) ; pause(0.25);
    
    if ~exist('jclips_tr_lbw2')
        % Delete feature dirs before runing temp_reg_lowbw again..
        try rmdir (fullfile(working_dir,'feature_*'), 's') ; catch; err=1; end
        [jclips_tr_lbw2, status.temporal2] = temporal_registration_lowbw(jtests, jclips_vrl, working_dir, 'Uncertainty', uncert);
       if status.temporal.error,
           jclips_tr_lbw2 = jclips_vrl;
       end

        % Per 12 Oct 06 e-mail:
       [jclips_tr_lbw2] = fix_temporal(jtests, jclips_tr_lbw2, 'StartPoint'); 
        save (status_file, 'status', 'model_selected','jclips','jclips_ft','jclips_tr_lbw','jclips_srl','jclips_vrl','jclips_tr_lbw2','jtests','cal_type','model_to_run','working_dir','ssf','uncert', 'jclips_sr_unfilt', 'jclips_vr_unfilt','tuncert_psnr','frac_samp','spatial_y','spatial_x','model_additional_options');
    end

    waitbar(6/9); pause(0.25);

    % update jclips with the final calibration report.
    jclips = jclips_tr_lbw2;
    
    if ~exist('crash')
        %Either this is the first time this function is running on the
        %clipset or the user wants to "start fresh" meaning that if PSNR
        %crashed then we need to make sure it removes all of its files
        %before running again.
        %Only calibrated version of PSNR will run in BVQM...
        warning off
        try
            delete([working_dir '\' 'temp_gclips.mat']);
        catch
            %File does not exist!
        end
        try
            delete([working_dir '\' 'temp_struct_orig.mat']);
        catch
            %File does not exist!
        end
        try
            delete([working_dir '\' 'tshifts.mat']);
        catch
            %File does not exist!
        end
        try
            delete([working_dir '\' 'results_psnr.mat']);
        catch
            %File does not exist!
        end
        try
            delete([working_dir '\' 'results_psnr.csv']);
        catch
            %File does not exist!
        end
        warning on
    end
    
    %Run PSNR Search and update jclips.
    if ~exist('jclips_temp')
        jclips_temp = bvqm_psnr_search_gclips(jtests,jclips,'results_psnr.csv',working_dir,'results_psnr',...
            'fraction_sampled',frac_samp,'spatial_uncertainty',spatial_x,spatial_y,...
            'temporal_uncertainty',tuncert_psnr,'silent','calibrated');
        delete([working_dir '\' 'results_psnr.mat']);
        delete([working_dir '\' 'results_psnr.csv']);
        save (status_file, 'status', 'model_selected','jclips','jtests','cal_type','model_to_run','working_dir','tuncert_psnr','frac_samp','spatial_x','spatial_y','jclips_temp','jclips_ft','jclips_tr_lbw','jclips_srl','jclips_vrl','jclips_tr_lbw2','ssf','uncert', 'jclips_sr_unfilt', 'jclips_vr_unfilt','model_additional_options');
    end
    
    waitbar(7/9); pause(0.25);
    
    %Valid Region after search is WRONG.  Needs to be recalculated.
    if ~exist('jclips_vr')
        [jclips_vr, status.valid, jclips_vr_unfilt] = valid_region(jtests, jclips_temp,'quiet');
        if status.valid, % error, continue without valid
            jclips_vr = jclips_temp;
        end
        save (status_file, 'status', 'model_selected','jclips','jtests','cal_type','model_to_run','working_dir','tuncert_psnr','frac_samp','spatial_x','spatial_y','jclips_temp','jclips_vr','jclips_vr_unfilt','jclips_ft','jclips_tr_lbw','jclips_srl','jclips_vrl','jclips_tr_lbw2','ssf','uncert', 'jclips_sr_unfilt','model_additional_options');
    end
    
    jclips = jclips_vr;
    waitbar(8/9); pause(0.25);
    
        if ~exist('model_op')
            [model_op, jclips, status] = run_model(model_to_run, jtests, jclips, working_dir, temp_file, cal_type, status, model_additional_options);
            save (status_file, 'status', 'model_selected','jclips','jtests','model_op','cal_type','model_to_run','working_dir','tuncert_psnr','frac_samp','spatial_x','spatial_y','jclips_temp','jclips_vr','jclips_vr_unfilt','jclips_ft','jclips_tr_lbw','jclips_srl','jclips_vrl','jclips_tr_lbw2','ssf','uncert', 'jclips_sr_unfilt','model_additional_options');
    end

    waitbar(9/9); pause(0.25);

    % if variable 'model_op' is empty, an error occurred runing the model.
    % Warn user about this.
    if isempty(model_op)
        model_crash(working_dir, path_sep, clip_dir);
%         
         model_op_err=1;
%         warnstr=sprintf('A fatal error occurred while running the model. \n You may have run out of memory.');
%         h=warndlg(warnstr, 'Error Occurred');
%         waitfor(h);
    else
        model_op_err=0;
        [s_rpt, d_rpt] = rpt_names(working_dir, path_sep, model_to_run);
        ssf = 0;
        bvqm_pc_shortrpt(model_selected,jclips,jtests,model_op,cal_type,working_dir,model_to_run,s_rpt,ssf,tuncert_psnr,status,[],[],[],[]);
        bvqm_pc_longrpt (model_selected,jclips,jtests,model_op,cal_type,working_dir,model_to_run,d_rpt,ssf,tuncert_psnr,status,[],jclips_vr_unfilt,jclips_sr_unfilt,[]);
    end 
end


%---------------------------------------------------
%--------- Fast PSNR Search 2 (frcal)
%---------------------------------------------------

if(full_ref_cal_then_psnr)
    cal_type = 'ITU-T J.144 followed by ITU-T J.340';
    if exist(status_file, 'file')
        load (status_file);
        crash=1;
    end
    
    if ~exist('tuncert_frcal')
        [output]=bvqm_pc_frcal_then_psnr_cal();
        if isempty(output)
            return
        else
            tuncert_psnr=output.uncert;
            frac_samp = output.frac_samp;
            spatial_y = output.spatial_y;
            spatial_x = output.spatial_x;
            tuncert_frcal = output.frcal_uncert;
            ssf = 0;
        end
        output = model_gui(model_to_run,jclips);
        if isempty(output)
            return
        end
        if strcmp(model_to_run, 'model_psnr')
            model_to_run = output.model_to_run;
        elseif(strcmp(model_to_run, 'model_psnr_vfd'))
            model_additional_options.which_psnr=output.model_additional_options.which_psnr;
            model_additional_options.t_uncert=output.model_additional_options.t_uncert;
        elseif(strcmp(model_to_run, 'model_vqm_vfd'))
            model_additional_options.t_uncert=output.model_additional_options.t_uncert ;
            jtests.viewing_distance = output.jtests.viewing_distance;
        end
    end
    
    %Close Calibration window, begin cal
    pause(0.25);
%    delete (handles.bvqm_pc_cal);
    
    wb = waitbar(0,'BVQM is running.  Please wait...');
    pause(0.25); delete (handles.bvqm_pc_cal);
    waitbar(1/7); pause(0.25);
    
    if ~exist('jclips_sr')
        % no status output available...
        [jclips_sr, status.spatial.failed, jclips_sr_unfilt] = spatial_registration(jtests, jclips, 'quiet');
        save (status_file, 'status', 'model_selected','jclips','jclips_sr','jtests','cal_type','model_to_run','working_dir','ssf','jclips_sr_unfilt','tuncert_psnr','frac_samp','spatial_y','spatial_x','tuncert_frcal','model_additional_options');
    end

    waitbar(2/7); pause(0.25);

    if ~exist('jclips_vr')
        [jclips_vr, status.valid, jclips_vr_unfilt] = valid_region(jtests, jclips_sr, 'Frequency', 0.5, 'quiet');
        errors_frame_based.valid_region=status.valid;
        if status.valid, % error, continue without valid
            jclips_vr = jclips_sr;
        end
        save (status_file, 'status', 'model_selected','jclips','jclips_sr','jclips_vr','jtests','cal_type','model_to_run','working_dir','ssf', 'jclips_sr_unfilt', 'jclips_vr_unfilt','tuncert_psnr','frac_samp','spatial_y','spatial_x','tuncert_frcal','model_additional_options');
    end

    waitbar(3/7); pause(0.25);

    if ~exist('jclips_tr')
        % no status output available
       [jclips_tr,status.temporal] = temporal_registration(jtests, jclips_vr, 'Uncertainty', tuncert_frcal, 'quiet');
       if status.temporal.error,
           jclips_tr = jclips_vr;
       end
       % Per 12 Oct 06 e-mail:
       [jclips_tr] = fix_temporal(jtests, jclips_tr, 'StartPoint'); 
        save (status_file, 'status', 'model_selected','jclips','jclips_sr','jclips_vr','jclips_tr','jtests','cal_type','model_to_run','working_dir','ssf', 'jclips_sr_unfilt', 'jclips_vr_unfilt','tuncert_psnr','frac_samp','spatial_y','spatial_x','tuncert_frcal','model_additional_options');
    end

    waitbar(4/7); pause(0.25);

    % update jclips with the final calibration report.
    jclips = jclips_tr;
    
    if ~exist('crash')
        %Either this is the first time this function is running on the
        %clipset or the user wants to "start fresh" meaning that if PSNR
        %crashed then we need to make sure it removes all of its files
        %before running again.
        %Only calibrated version of PSNR will run in BVQM...
        warning off
        try
            delete([working_dir '\' 'temp_gclips.mat']);
        catch
            %File does not exist!
        end
        try
            delete([working_dir '\' 'temp_struct_orig.mat']);
        catch
            %File does not exist!
        end
        try
            delete([working_dir '\' 'tshifts.mat']);
        catch
            %File does not exist!
        end
        try
            delete([working_dir '\' 'results_psnr.mat']);
        catch
            %File does not exist!
        end
        try
            delete([working_dir '\' 'results_psnr.csv']);
        catch
            %File does not exist!
        end
        warning on
    end
    
    %Run PSNR Search and update jclips.
    if ~exist('jclips_temp')
        jclips_temp = bvqm_psnr_search_gclips(jtests,jclips,'results_psnr.csv',working_dir,'results_psnr',...
            'fraction_sampled',frac_samp,'spatial_uncertainty',spatial_x,spatial_y,...
            'temporal_uncertainty',tuncert_psnr,'silent','calibrated');
        delete([working_dir '\' 'results_psnr.mat']);
        delete([working_dir '\' 'results_psnr.csv']);
        save (status_file, 'status', 'model_selected','jclips','jclips_sr','jclips_vr','jclips_tr','jtests','cal_type','model_to_run','working_dir','ssf', 'jclips_sr_unfilt', 'jclips_vr_unfilt','tuncert_psnr','frac_samp','spatial_y','spatial_x','tuncert_frcal','jclips_temp','model_additional_options');
    end
    
    waitbar(5/7); pause(0.25);
    
    %Valid Region after search is WRONG.  Needs to be recalculated.
    if ~exist('jclips_vr')
        [jclips_vr, status.valid, jclips_vr_unfilt] = valid_region(jtests, jclips_temp,'quiet');
        if status.valid, % error, continue without valid
            jclips_vr = jclips_temp;
        end
        save (status_file, 'status', 'model_selected','jclips','jclips_sr','jclips_vr','jclips_tr','jtests','cal_type','model_to_run','working_dir','ssf', 'jclips_sr_unfilt', 'jclips_vr_unfilt','tuncert_psnr','frac_samp','spatial_y','spatial_x','tuncert_frcal','jclips_temp','model_additional_options');
    end
    
    jclips = jclips_vr;
    waitbar(6/7); pause(0.25);
    
        if ~exist('model_op')
            [model_op, jclips, status] = run_model(model_to_run, jtests, jclips, working_dir, temp_file, cal_type, status, model_additional_options);
            save (status_file, 'status','model_op', 'model_selected','jclips','jclips_sr','jclips_vr','jclips_tr','jtests','cal_type','model_to_run','working_dir','ssf', 'jclips_sr_unfilt', 'jclips_vr_unfilt','tuncert_psnr','frac_samp','spatial_y','spatial_x','tuncert_frcal','jclips_temp','model_additional_options');
    end

    waitbar(7/7); pause(0.25);

    % if variable 'model_op' is empty, an error occurred runing the model.
    % Warn user about this.
    if isempty(model_op)
        model_crash(working_dir, path_sep, clip_dir);
         model_op_err=1;
%         warnstr=sprintf('A fatal error occurred while running the model. \n You may have run out of memory.');
%         h=warndlg(warnstr, 'Error Occurred');
%         waitfor(h);
    else
        model_op_err=0;
        [s_rpt, d_rpt] = rpt_names(working_dir, path_sep, model_to_run);
        bvqm_pc_shortrpt(model_selected, jclips_tr, jtests, model_op, cal_type,working_dir, model_to_run,s_rpt,ssf,tuncert_psnr,status, [], [], [], []);
        bvqm_pc_longrpt (model_selected, jclips_tr, jtests, model_op, cal_type,working_dir, model_to_run,d_rpt,ssf,tuncert_psnr, status, [], jclips_vr_unfilt, jclips_sr_unfilt, []);
%        close(wb);
    end
    
end


    
%---------------------------
if ~exist('model_op')
    model_op_err=1;
end
    
if model_op_err==0
%     bvqm_pc_report('bvqm_pc_report_passer',working_dir, s_rpt, d_rpt, clip_dir);
%     bvqm_pc_report(working_dir, s_rpt, d_rpt, clip_dir);
     bvqm_pc_report([],working_dir,[], s_rpt,[], d_rpt,[], clip_dir);
     pause(0.25);
     close(wb);
end


%------------------------
% clean-up when the model crashes:
%------------------------
function model_crash(working_dir, path_sep, clip_dir)
        model_op_err=1;
        warnstr=sprintf('A fatal error occurred while running the model. \n You may have run out of memory. \n\n BVQM will exit.');
        h=warndlg(warnstr, 'Error Occurred');
        waitfor(h);

        warning off all;
        try  rmdir (strcat(working_dir, path_sep,'feature_*'), 's') ;   catch; err=1; end
        try delete (strcat(working_dir, path_sep,'bvqm-status.mat')) ;  catch; err=1; end
        try delete (strcat(working_dir, path_sep,'bvqm-opt_*.mat'));    catch; err=1; end
        try delete (strcat(working_dir, path_sep,'jclips.mat'));        catch; err=1; end

        fclose('all');
%        delete (handles.bvqm_pc_cal);
        close all;

%------------------------
% clean-up when the Excel crashes:
%------------------------
function excel_crash(working_dir, path_sep, clip_dir)

        warning off all;
        try  rmdir (strcat(working_dir, path_sep,'feature_*'), 's') ;   catch; err=1; end
        try delete (strcat(working_dir, path_sep,'bvqm-status.mat')) ;  catch; err=1; end
        try delete (strcat(working_dir, path_sep,'bvqm-opt_*.mat'));    catch; err=1; end
        try delete (strcat(working_dir, path_sep,'jclips.mat'));        catch; err=1; end

        fclose('all');
%        delete (handles.bvqm_pc_cal);
        close all;

%----------------------------------------------
%---% set-up file names for report files
%-----------------------------------------------
function [s_rpt, d_rpt] = rpt_names(working_dir, path_sep, model_to_run)

clk=fix(clock);
yr=sprintf('%4.0f',clk(1)) ; mo=sprintf('%02.0f',clk(2)); da=sprintf('%02.0f',clk(3));
hr=sprintf('%02.0f',clk(4)); mn=sprintf('%02.0f',clk(5)); sc=sprintf('%02.0f',clk(6));

date=strcat(yr,'-',mo,'-',da);
time=strcat(hr,'_',mn);

s_rpt=strcat(working_dir,path_sep,'bvqm_summary_report-',model_to_run,'-',date,'@',time, '.txt');
d_rpt=strcat(working_dir,path_sep,'bvqm_detail_report-',model_to_run,'-',date,'@',time, '.txt');

%---------------------------------------------
% --- run the model function that was selected:
% update 'status' variable, with field 'model'
%----------------------------------------------
function [model_op, jclips, status] = run_model(model_to_run, jtests, jclips, working_dir, temp_file, cal_type, status, add_ops)

if strcmp(model_to_run, 'model_general')
    [model_op, status.model] = model_general(jtests, jclips, working_dir, 'quiet');
%     if status.rm_mg ==1
%         warnmsg='An error has occurred running the model. \n Please check calibration settings.';
%         h=warndlg (sprintf(warnmsg), 'Error');
%         waitfor(h);
%         return
%     end
end

if strcmp(model_to_run, 'model_lowbw')
    [model_op, status.model] = model_lowbw(jtests, jclips, working_dir, temp_file, 'quiet');
end

if strcmp(model_to_run, 'model_fastlowbw')
    [model_op, status.model] = model_fastlowbw(jtests, jclips, working_dir, temp_file, 'quiet');
end

if strcmp(model_to_run, 'model_developers')
    [model_op, status.model] = model_developers(jtests, jclips, working_dir, 'quiet');
end

if strcmp(model_to_run, 'model_television')
    [model_op, status.model] = model_television(jtests, jclips, working_dir, 'quiet');
end

if strcmp(model_to_run, 'model_videoconferencing')
    [model_op, status.model] = model_videoconferencing(jtests, jclips, working_dir, 'quiet');
end

if strcmp(model_to_run, 'psnr_235')
    [model_op] = model_psnr(jtests, jclips, 'quiet', '235');
    status.model = 0;
end

if strcmp(model_to_run, 'psnr_255')
    [model_op] = model_psnr(jtests, jclips, 'quiet', '255');
    status.model = 0;
end

if strcmp(model_to_run, 'model_psnr_vfd')
    [model_op jclips_temp] = bvqm_model_psnr_vfd(jtests, jclips, working_dir,  ...
        'quiet', 't_uncert', add_ops.t_uncert, 'peak', add_ops.which_psnr);
    jclips = jclips_temp;
    status.gain_offset = 0;
    status.model = 0;
end

if strcmp(model_to_run, 'model_vqm_vfd')
    try
        [model_op] = bvqm_model_vqm_vfd(jtests, jclips, working_dir, 't_uncert', add_ops.t_uncert);
        status.model = 0;
    catch err
        errordlg('Model calculation failed. Neural Net toolbox is required. BVQM exiting.');
        return;
    end
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% % --- Executes on button press in pushbutton2_SelectDir.
% function pushbutton2_SelectDir_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton2_SelectDir (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
%
% [pathname] = uigetdir('','Specify File Path');
% set(handles.text5_WorkingDirectory, 'String', pathname)


% --- Executes on button press in pushbutton3_Back.
function pushbutton3_Back_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3_Back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
working_dir = get(handles.text5_WorkingDirectory, 'String');
%clip_dir = mat2str(cell2mat(get(handles.text7_clip_dir, 'String')));
clip_dir = cell2mat(get(handles.text7_clip_dir, 'String'));

if ispc ; path_sep = '\'; else; path_sep='/'; end
jcpf=strcat(working_dir, path_sep, 'jclips.mat');  %jclips path/file name

% want to query, then close window:
button = questdlg('Discard all previous data?', 'BVQM Question:', 'No');
waitfor(button);

if (strcmp(button, 'No'))  % keep data!
    % recreate bvqm-vcd.mat if it doesn't exist
    % (bvqm-vcd.mat doesn't exist if pressing back after crash recovery)
    vcd_fn = fullfile(working_dir, 'bvqm-vcd.mat'); % video clip data filename
    if ~exist(vcd_fn, 'file')
        load (fullfile(working_dir, 'jclips'));
        for jjj=1:size(jclips,2)
        vcd(jjj).scene=cell2mat(jclips(jjj).scene);
        vcd(jjj).hrc=cell2mat(jclips(jjj).hrc);
        vcd(jjj).filename=cell2mat(jclips(jjj).file_name);
        save(vcd_fn, 'vcd');
        end
    end
    % now exit this screen
    bvqm_pc(working_dir, clip_dir, button);
    pause(0.25);
    delete(handles.bvqm_pc_cal); % close;
end

if (strcmp(button, 'Cancel'))
    return
end

if (strcmp(button, 'Yes'))
    % silence any warning messages that occur (from deleting
    % nonexistant files e.g. bvqm-status.mat
    warning off all
    try  rmdir (strcat(working_dir, path_sep,'feature_*'), 's') ;   catch; err=1; end
    try delete (strcat(working_dir, path_sep,'bvqm-status.mat')) ;  catch; err=1; end
    try delete (strcat(working_dir, path_sep,'bvqm-opt_*.mat'));    catch; err=1; end
    try delete (strcat(working_dir, path_sep,'jclips.mat'));        catch; err=1; end
    try delete (fullfile(working_dir, 'bvqm-vcd.mat'));             catch; err=1 ; end
%    bvqm_pc('bvqm_pc_OpeningFcn', working_dir, clip_dir, button);
    bvqm_pc(working_dir, clip_dir, button);
    pause(0.25);
    delete(handles.bvqm_pc_cal); % close;
end

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over radiobutton1.
function radiobutton1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radiobutton9.
function radiobutton9_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton9


%  ---- Uncertainty function - querys for a value for uncert:
%function uncert = get_uncert
function [output] = get_uncert(qss)
% ask for the value for uncertainty. must be >= 0
uncert = -1.0;
[output]=bvqm_pc_uncertss(qss);
    

% --- Executes on button press in pushbutton4_exit.
function pushbutton4_exit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

warning off all;
if ispc ; path_sep = '\'; else; path_sep='/'; end
working_dir = get(handles.text5_WorkingDirectory, 'String');
clip_dir = mat2str(cell2mat(get(handles.text7_clip_dir, 'String')));

b1='Save & Exit'; b2='Exit'; b3='Cancel';

really= questdlg('Save intermediate structures before quit?', 'Confirm Save/Exit', b1,b2,b3,b2);
waitfor(really);
 
warning off all;
 if strcmp(really, b2)  
    try delete (fullfile(working_dir, 'bvqm-vcd.mat')); catch ; err=1 ; end  % video clip data filename
    try  rmdir (strcat(working_dir, path_sep,'feature_*'), 's') ;   catch; err=1; end
    try delete (strcat(working_dir, path_sep,'jclips.mat'));        catch; err=1; end
    try delete (strcat(working_dir, path_sep,'bvqm-status.mat')) ;  catch; err=1; end
    try delete (strcat(working_dir, path_sep,'bvqm-opt_*.mat'));    catch; err=1; end
    
    delete(handles.bvqm_pc_cal); % close all; %('bvqm_pc_report');
    fclose('all');
    pause(0.25);
 end
 
 if strcmp(really, b1)  % don't delete anything -- just exit.
    msg={'Structures are saved in the working directory:'; working_dir};
    h=msgbox(msg, 'Structures Saved', 'help');
    waitfor(h);
 
    delete(handles.bvqm_pc_cal); % close all; %('bvqm_pc_report');
    fclose('all');
    pause(0.25);
 end
 
function output_new = model_gui(model_to_run,jclips)
% if chose PSNR, determine peak white value.
output_new = model_to_run;
if strcmp(model_to_run, 'model_psnr')
    [peak_value] = questdlg('What would you like to use for the peak signal:  255 (computer white) or 235 (ITU-R Rec. BT-601 white)?', 'PSNR Model', '255', '235', []);

    if isempty(peak_value)
        return
    elseif strcmp(peak_value,'235'),
        output_new.model_to_run = 'psnr_235';
    else
        output_new.model_to_run = 'psnr_255';
    end
end

if(strcmp(model_to_run, 'model_psnr_vfd'))
    [output]=bvqm_pc_psnr_vfd_model();
    if isempty(output)
        return
    else
        output_new.model_additional_options.which_psnr=output.which_psnr;
        output_new.model_additional_options.t_uncert=output.uncert ;
    end
end

if(strcmp(model_to_run, 'model_vqm_vfd'))
    [output]=bvqm_pc_vqm_vfd_model(jclips(1).image_size);
    if isempty(output)
        return
    else
        viewing_distance=output.vdist;
        output_new.model_additional_options.t_uncert=output.uncert ;
        output_new.jtests.viewing_distance = viewing_distance;
    end
end
