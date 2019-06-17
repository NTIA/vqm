function varargout = bvqm_pc(varargin)
% BVQM_PC Application M-file for bvqm_pc.fig
%   BVQM_PC, by itself, creates a new BVQM_PC or raises the existing
%   singleton*.
%
%   H = BVQM_PC returns the handle to a new BVQM_PC or the handle to
%   the existing singleton*.
%
%   BVQM_PC('CALLBACK',hObject,eventData,handles,...) calls the local
%   function named CALLBACK in BVQM_PC.M with the given input arguments.
%
%   BVQM_PC('Property','Value',...) creates a new BVQM_PC or raises the
%   existing singleton*.  Starting from the left, property value pairs are
%   applied to the GUI before bvqm_pc_OpeningFunction gets called.  An
%   unrecognized property name or invalid value makes property application
%   stop.  All inputs are passed to bvqm_pc_OpeningFcn via varargin.
%
%   *See GUI Options - GUI allows only one instance to run (singleton).
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2000-2002 The MathWorks, Inc.


% Edit the above text to modify the response to help bvqm_pc

%
% WRITTEN BY : Mark A. McFarland, P.E.    303/497-4132
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',          mfilename, ...
    'gui_Singleton',     gui_Singleton, ...
    'gui_OpeningFcn',    @bvqm_pc_OpeningFcn, ...
    'gui_OutputFcn',     @bvqm_pc_OutputFcn, ...
    'gui_LayoutFcn',     [], ...
    'gui_Callback',      []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    varargout{1:nargout} = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before bvqm_pc is made visible.
function bvqm_pc_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bvqm_pc (see VARARGIN)


% Choose default command line output for bvqm_pc
handles.output = hObject;

% set first time entereing nonstandard filename variable to 0:
set(handles.text41, 'Value', 0);

%select the working directory, but only on first run:
wp_check_res=1;
wp_check_clp=1;
exit_flag=0;  % sets exit flag - bvqm will exit if 'cancel' is selected during directory selections.
% if cell2mat(varargin(4))
if isequal(size(varargin,2), 3)
    working_dir=(varargin{1});
    set(handles.text40, 'String', working_dir);
    clip_dir=(varargin{2});
    set(handles.clip_directory_text3, 'String', clip_dir);
    button=(varargin{3});   % no=keep data, yes=don't keep data

    % Populate the listbox1
    load_listbox(clip_dir,handles);
    % Remove any blank lines lingering in lb2:
    set(handles.listbox2, 'String', '')
    
    % Repopulate test data info: fps, rows, cols, vid std. 
%    if exist(fullfile(working_dir, 'jclips.mat', 'file'))  % jclips alread exists when going 'back' from cal screen
    if strcmp(button, 'No')  %If we keep data returning from cal screen.
    load(fullfile(working_dir, 'jclips.mat'));
        fps=jclips(1).fps;
        rows=jclips(1).image_size.rows;
        cols=jclips(1).image_size.cols;
        vid_std=jclips(1).video_standard;
        if strcmp(vid_std, 'interlace_lower_field_first'); vs_idx=1; end
        if strcmp(vid_std, 'interlace_upper_field_first'); vs_idx=2; end
        if strcmp(vid_std, 'progressive'); vs_idx=3; end
        set(handles.fps_edit, 'String', fps);
        set(handles.image_size_rows_edit, 'String', rows);
        set(handles.image_size_cols_edit, 'String', cols);
        set(handles.video_std_popupmenu1, 'Value', vs_idx);

        %Repoupulate selected files:
        %    set(handles.listbox2, 'String', jclips.file_name)
        cnt=1;
        num_clips=size(jclips, 2);
        while cnt <= num_clips
            lb2_fn(cnt) = (jclips(cnt).file_name);
            cnt=cnt+1;
            set(handles.listbox2, 'String', lb2_fn);
        end
    end
else   
    while isequal(wp_check_res, 1)
        [working_dir] = uigetdir('','BVQM: Select a RESULTS directory:');
        if isequal(working_dir,0);
            working_dir=pwd;
            exit_flag=1;
        end
        
%         [stat, msg, msg_id] = fileattrib(working_dir);
%         msg;

        try
            cannot_save=0;
            sfile=fullfile(working_dir,'jjj.mat');
            save(sfile, 'working_dir')
        catch
            cannot_save=1;
%             warnstr={'Results directory must have write permission.';'Please select new directory.'};
%             h=warndlg(warnstr, 'Write Permission Required');
%             waitfor(h);
%             working_dir=0;
%             wp_check_res=1;
        end
                
%         %    if isequal(wp_check_res, 0); % only display warning once
%        if  isequal(msg.UserWrite, 0)
        if  isequal(cannot_save, 1)
            warnstr={'Results directory must have write permission.';'Please select new directory.'};
            h=warndlg(warnstr, 'Write Permission Required');
            waitfor(h);
            working_dir=0;
            wp_check_res=1;
        else
            delete (sfile)
             if ~isempty(regexp(working_dir, ' '))  % no spaces in path--XX- Sept06- spaces allowed...
%                 warnstr={'Results directory may not have spaces in path.';'Please select new directory.'};
%                 h=warndlg(warnstr, 'No Spaces allowed in Path');
%                 waitfor(h);
%                 working_dir=0;
%                 wp_check_res=1;
                 wp_check_res=0;
             else
                 wp_check_res=0;
             end
        end
    end
    
        
    % NOTE: text40 not visible - just used to store data...
    if working_dir==0
        working_dir=pwd; 
        exit_flag=1 ; % set exit flag b/c cancel was pressed.
%         warnstr={'No results directory was selected'; 'Results directory set to current working directory.'};
%         h=warndlg(warnstr, 'Results directory not entered!');
%         waitfor(h);
%         fclose('all');
%         delete(handles.bvqm_pc);
    end
%   else
        set(handles.text40, 'String', working_dir);
        load_listbox(working_dir,handles);

        % set-up clip directory:
        % if cancel was selected already, don't ask for clip directory
        if isequal(exit_flag, 1)
            wp_check_clp=0;
            clip_dir=pwd;
        end
        
        while isequal(wp_check_clp, 1)
            [clip_dir] = uigetdir('','BVQM: Select a CLIP directory:');
            if isequal(clip_dir,0);
                clip_dir=pwd;
                exit_flag=1;  % set exit flag -- cancel was selected
            end
            [stat, msg, msg_id] = fileattrib(clip_dir);
            msg;
            
            % Just set msg.UserWrite=1 - no write permission needed now
            msg.UserWrite=1;

            %    if isequal(wp_check_res, 0); % only display warning once
            if  isequal(msg.UserWrite, 0)
                warnstr={'Clip directory must have write permission.';'Please select new directory with write permission.'};
                h=warndlg(warnstr, 'Write Permission Required');
                waitfor(h);
                working_dir=0;
                wp_check_clp=1;
            else
                if ~isempty(regexp(clip_dir, ' '))  % no spaces in path
%                     warnstr={'Clip directory may not have spaces in path.';'Please select new directory.'};
%                     h=warndlg(warnstr, 'No Spaces allowed in Path');
%                     waitfor(h);
%                     clip_dir=0;
%                     wp_check_clp=1;
                     wp_check_clp=0;
                 else
                    wp_check_clp=0;
                    set(handles.clip_directory_text3,'String', clip_dir);
                end
            end
        end
 %  end


    if clip_dir==0
         clip_dir=pwd;
         exit_flag=1;  % set exit flag b/c cancel was selected.
%          warnstr={'No clip directory was selected'; 'Clip directory set to current working directory.'};
%          h=warndlg(warnstr, 'Clip directory not entered!');
%          waitfor(h);
    end
    
% Populate the listbox
  load_listbox(clip_dir,handles);

% Remove any blank lines lingering in lb2:
   set(handles.listbox2, 'String', '')

    % Check to see if we are recovering from a crash - but only if cancel
    % wasn't selected during the directory selection.
    if isequal(exit_flag,0)
        %    Assume we are recovering from crash if 'jclips.mat' already exists in the working directory:
        if ispc ; path_sep = '\'; else; path_sep='/'; end
        jc_file= strcat(working_dir, path_sep, 'jclips.mat');
        if exist(jc_file, 'file')
            qstring={'A clips structure was found in the results directory.';...
                'Continue using recovered structure, or start fresh?'};

            button = questdlg(qstring, 'Clip structure found!','Start Fresh', 'Continue','Continue');
            waitfor(button);
            if strcmp(button, 'Continue')
                cr=1; % 1 signifies crash recovery
                pass={working_dir, cr};
%                bvqm_pc_cal('bvqm_pc_cal_Opening_Fcn', pass);
                bvqm_pc_cal([], pass);
                pause(0.5);
                delete(handles.bvqm_pc);
                pause(0.5);
            end
            if strcmp(button, 'Start Fresh')
                % silence any warning messages that occur (from deleting nonexistant files)
                warning off all;
                try  rmdir (strcat(working_dir, path_sep,'feature_*'), 's') ;   catch; err=1; end
                try  rmdir (strcat(working_dir, path_sep,'vfd_results_*'), 's');catch; err=1; end
                try delete (strcat(working_dir, path_sep,'vfd_*.mat'));         catch; err=1; end
                try delete (strcat(working_dir, path_sep,'bvqm-opt_*.mat'));    catch; err=1; end
                try delete (strcat(working_dir, path_sep,'jclips.mat'));        catch; err=1; end
                try delete (strcat(working_dir, path_sep,'bvqm-status.mat')) ;  catch; err=1; end
                try delete (strcat(working_dir, path_sep,'bvqm-linked.mat')) ;  catch; err=1; end
                try delete (fullfile(working_dir, 'bvqm-vcd.mat')); catch; err=1 ; end % video clip data filename

                try delete (strcat(working_dir, path_sep,'vqm_*.mat'));         catch; err=1; end
                try delete (strcat(working_dir, path_sep,'vfd_results.csv'));  catch; err=1; end
                try delete (strcat(working_dir, path_sep,'vfd.mat'));           catch; err=1; end
                
                % Display main screen, logo:
                x = vqm_logo;
                %x = imread('vqm_logo.jpg');
                axes(handles.axes1); % makes axis1 current
                imagesc(x);
                set(gca,'Visible','off')
            end
        end
    else % exit_flag
        fclose('all');
        delete(handles.bvqm_pc);
        return
    end
%     % Display main screen, logo:
%     x = vqm_logo;
%     %x = imread('vqm_logo.jpg');
%     axes(handles.axes1); % makes axis1 current
%     imagesc(x);
%     set(gca,'Visible','off')
end

try
    % Display main screen, logo:
    x = vqm_logo;
    %x = imread('vqm_logo.jpg');
    axes(handles.axes1); % makes axis1 current
    imagesc(x);
    set(gca,'Visible','off')
catch
end

% % set version number:
% set(handles.Name, 'String', 'BVQM_V1.4')


% --- Outputs from this function are returned to the command line.
 function varargout = bvqm_pc_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function axis1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axis1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% ------------------------------------------------------------
% Callback for listbox - populate listbox1:
% ------------------------------------------------------------
function varargout = listbox1_Callback(h, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
% Disable Change Directory ability after a video file is selected


% [clip_dir] = uigetdir('','BVQM: Select a Clip directory:');
% set(handles.clip_directory_text3,'String', clip_dir);
     
clip_dir = get(handles.clip_directory_text3,'String');
% get(handles.bvqm_pc,'SelectionType');

if ispc ; path_sep = '\'; else; path_sep='/'; end

if strcmp(get(handles.bvqm_pc,'SelectionType'),'open')
    listbox1_index = get(handles.listbox1,'Value');
    file_list = get(handles.listbox1,'String');
    clipdir_item = file_list{listbox1_index};
    clip_pathfn=strcat(clip_dir, path_sep, clipdir_item);
    %    fileattrib(clipdir_item)

    % if we're dealing with a directory....
    if  handles.is_dir(handles.sorted_index(listbox1_index))
        %   if  is_dir(sorted_index(listbox1_index))
        chosen_list=get(handles.listbox2,'String');
        %                     dirrrrrrrrrrrrr=1
        if length(chosen_list) ~=0
            message = {'               Cannot Change Directory!';...
                'Please choose video files from only one directory.';...
                'To choose a new directory, click "remove all" first.'};
            warndlg(message);
        else % a directory was chosen (. or ..)
            % cd (clipdir_item)
            % load_listbox(clipdir_item,handles);
            load_listbox(clip_pathfn,handles);
        end
    else
 %       [path,name,ext,ver] = fileparts(clipdir_item);
        set(handles.text6, 'String', clipdir_item) % located in clip data box..

    end
end

% ------------------------------------------------------------
% Read the current directory and sort the names; display in listbox1
% This occurs when a directory is selected in lb1.
% ------------------------------------------------------------
function load_listbox(dir_path,handles)

% Here - if there is a '..' or '.' at the end of the path, need to manually (!)
% reset path without trailing '..' or '.' .
%  e.g. path name of '/home/mark/projects/3141/bvqm/scratchdir14/..' should
%  be displayed as '/home/mark/projects/3141/bvqm'    and
%  '/home/mark/projects/3141/bvqm/.' should be '/home/mark/projects/3141/bvqm'
% 
% NOTE : this section not currently in use -
%        '.' and '..' problem fixed by forcing the selection of a clip dir
%        at same time that the working_dir is selected.

if ispc ; path_sep = '\'; else; path_sep='/'; end

[path, dd] = strtok(dir_path, '..');

cda=0;  % directory name not rebuilt...yet
clip_dir='';
if strcmp(dd, '..')  % fix path name - '..' selected
    tmp_dir_path = dir_path;
    num_slashes = size(strfind(tmp_dir_path, path_sep), 2);
    slash_count=1;
%     clip_dir='';
    while (slash_count <= num_slashes)
        [token, rem] = strtok(tmp_dir_path, path_sep);
        t(slash_count) = {token};
        r(slash_count) = {rem};
        tmp_dir_path = rem;

        % rebuild clip dir path:
        if slash_count <= num_slashes -2
            clip_dir = strcat(clip_dir, path_sep, t(slash_count));
            cda=1; % marker if a new clip directory was assigned
        end
        slash_count = slash_count +1;
    end
else
    if strcmp(dd, '.') % fix path name -  '.' selected
        if iscell(dir_path); dir_path=cell2mat(dir_path); end
        [tok,rem] = strtok(dir_path, '.');
        tok_size=size(tok,2) - 1;
        clip_dir = tok(1:tok_size)
        cda=1;
    end

end


if iscell(dir_path); dir_path=cell2mat(dir_path); end

% Only display avi any yuv files:
dir_struct_avi = dir(fullfile(dir_path,'*.avi'));
dir_struct_yuv = dir(fullfile(dir_path,'*.yuv'));
dir_struct=cat(1,dir_struct_avi,dir_struct_yuv);


% REMOVE '.' and '..' entries:
num_els = size(dir_struct, 1);
count=1; dir_struct_tmp = dir_struct;
while count <= num_els
    is_dir = dir_struct(count).isdir;
    if is_dir
        dir_struct(count)='';
        count=0;
        num_els=num_els-1;
    end
    count=count+1;
end

 [status, fa] = fileattrib(dir_path);
 
 if ~cda  % if clip directory name didn't have to be rebuilt
     clip_dir=fa.Name;
     %clip_dir=strcat(fa.Name,dir_path)
 end

 [sorted_names,sorted_index] = sortrows({dir_struct.name}');
 % odd use of handles, but it works well ... may want to rewrite later
 handles.file_names = sorted_names;
 handles.is_dir = [dir_struct.isdir];
 handles.sorted_index = [sorted_index];
guidata(handles.bvqm_pc,handles)
% dot=strcmp(handles.file_names, '.')
% dotdot=strcmp(handles.file_names, '..')

set(handles.listbox1,'String',handles.file_names, 'Value',1)
set(handles.clip_directory_text3,'String',clip_dir)


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function bvqm_pc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bvqm_pc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Add the current directory to the path, as the pwd might change thru' the
% gui. Remove the directory from the path when gui is closed
% (See bvqm_pc_DeleteFcn)
setappdata(hObject, 'StartPath', pwd);
addpath(pwd);
% set(bvqm_pc, 'CloseRequestFcn', 'bvqm_pc_CloseRequestFcn')

% --- Executes when user attempts to close bvqm_pc.
function bvqm_pc_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to bvqm_pc_imgsz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Hint: delete(hObject) closes the figure
warning off all;
% Remove any lingering files...
working_dir=get(handles.text40, 'String');
clip_dir=get(handles.clip_directory_text3, 'String');
if ispc ; path_sep = '\'; else; path_sep='/'; end

    try  rmdir (strcat(working_dir, path_sep,'feature_*'), 's') ;   catch; err=1; end
    try delete (strcat(working_dir, path_sep,'bvqm-status.mat')) ;  catch; err=1; end 
    try delete (strcat(working_dir, path_sep,'bvqm-opt_*.mat'));    catch; err=1; end
    try delete (strcat(working_dir, path_sep, 'bvqm-parsed-*.mat')); catch; err=1; end
%    try delete (strcat(working_dir, path_sep, 'bvqm-test_data.mat')); catch ; err=1 ; end
    try delete (strcat(working_dir, path_sep,'jclips.mat'));        catch; err=1; end
    try delete (strcat(clip_dir, path_sep, '%bvqm_*_*.*'));      catch; err=1; end
    try delete (strcat(working_dir, path_sep, 'bvqm-ptrunk.mat')); catch; err=1; end
    
     try load (fullfile(working_dir, 'bvqm-linked.mat'))
        num_linked=size(linked_file,2);
        count=1;
        while count <= num_linked
            try delete (linked_file{count}) ; catch; err=1 ; end
            count=count+1;
        end
        delete(fullfile(working_dir, 'bvqm-linked.mat'));
    catch
        nosuchfile=1;
     end
   
     fclose('all');
     delete(hObject);



% --- Executes during object deletion, before destroying properties.
function bvqm_pc_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to bvqm_pc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Remove the directory added to the path in the bvqm_pc_CreateFcn.
if isappdata(hObject, 'StartPath')
    rmpath(getappdata(hObject, 'StartPath'));
end
fclose('all');
%warning off all

% % Remove any lingering files...
%     try  rmdir (strcat(working_dir, path_sep,'feature_*'), 's') ;   catch; err=1; end
%     try delete (strcat(working_dir, path_sep,'bvqm-status.mat')) ;  catch; err=1; end 
%     try delete (strcat(working_dir, path_sep,'bvqm-opt_*.mat'));    catch; err=1; end
%     try delete (strcat(working_dir, path_sep, 'bvqm-parsed-*.mat')); catch; err=1; end
% %    try delete (strcat(working_dir, path_sep, 'bvqm-test_data.mat')); catch ; err=1 ; end
%     try delete (strcat(working_dir, path_sep,'jclips.mat'));        catch; err=1; end
%     try delete (strcat(clip_dir, path_sep, '%bvqm_*_*.*'));      catch; err=1; end
%     try delete (strcat(working_dir, path_sep, 'bvqm-ptrunk.mat')); catch; err=1; end
%     
%      working_dir=get(handles.text40, 'String');
%      try load (fullfile(working_dir, 'bvqm-linked.mat'))
%         num_linked=size(linked_file,2);
%         count=1;
%         while count <= num_linked
%             try delete (linked_file{count}) ; catch; err=1 ; end
%             count=count+1;
%         end
%         delete(fullfile(working_dir, 'bvqm-linked.mat'));
%     catch
%         nosuchfile=1;
%     end

% --- Executes on button press in add_pushbutton2.
function add_pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to add_pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ispc ; path_sep = '\'; else; path_sep='/'; end
working_dir = get(handles.text40, 'String');
clip_dir = get(handles.clip_directory_text3, 'String');

noadd=0;
file_list = get(handles.listbox1, 'String');
lb1_index = get(handles.listbox1,'Value');
num_files=length(lb1_index); %number of files selected

% Don't allow user to "add" a file if the listbox is empty
if isempty(file_list)
    return
end

for file_count=1:num_files  %for each selected file
    item=char(file_list{ lb1_index(file_count) }); %current selection
    hold_lb2 = get(handles.listbox2, 'String');
    [rows,cols] = size(hold_lb2);
    
    % is the selected file "readable"?
    [fid,message]=fopen(fullfile(clip_dir,item));
    if fid ==-1
        errorstr={'Error reading selected file:';' '; item;message;' '; 'File will not be added.'};
        h=errordlg(errorstr, 'File Read Error');
        waitfor(h);
        noadd=1;
        set(handles.listbox1, 'Value', lb1_index+1);
        return
    end       

    % logic to identify YUV file errors: file size not div. by 2:
    [pre,ext] = strtok(item, '.');
    if strcmp(ext, '.yuv')
        [fid,message] = fopen(strcat(clip_dir, path_sep, item));
        fseek(fid,0,'eof');
        fsize = ftell(fid);
        remm=rem(fsize,2);
        if remm % if ther eis  a remainder~= 0
            errstr='YUV file is not of a valid size.  It will not be added.';
            h=errordlg(errstr, 'YUV file error', 'on');
            waitfor(h)
            noadd=1;
            return
        end
    end

    % What type of file ?
    %   std:xxx_xxx_xxx.[avi,yuv]
    nsn=0;  % assume std file name to begin:
    % Step 1: find nested multiple _s
    double_=size(findstr(item,'__'),2);
    triple_=size(findstr(item,'___'),2);
    quad_  =size(findstr(item,'____'),2);
    quin_  =size(findstr(item,'_____'),2);
    tmp=size(findstr(item,'_'),2);
    if isequal(tmp,0); zero_ = 1; else; zero_=0; end

    if  (zero_ | double_ | triple_ | quad_ | quin_ )
        nsn=1;   % nsn = non standard filename
        std_match=0;
    end

    % if more than 2 underscores, not a std filename
    if tmp>2
        nsn=1;
        std_match=0;
    end
    
    %takes care of _xxx.xxx and xxx_.xxx
    if strncmp(item, '_', 1) | regexp(item, '_\.')
        nsn=1;
        std_match=0;
    end

    % what if there are multiple (or no) '.'s in the filename:
    num_dots=size(findstr(item,'.'),2);
    if ~isequal(size(findstr(item,'.'),2),1)
        nsn=1;
        std_match=0;
    end

    % if there are spaces in the clip name:
    if size(findstr(item, ' '),2) >1
        nsn=1;
        std_match=0;
    end

    % % Step 2: is the extention .avi or .yuv?
    if strcmp(ext, '.yuv'); ext_yuv=1;else; ext_yuv=0; end
    if strcmp(ext, '.avi');  ext_avi=1; else; ext_avi=0; end

    if ~(ext_yuv | ext_avi)
        noadd=1;
        message={'Only yuv and avi files are allowed.';'Please ensure that clip files have the proper extension.';...
            '';'This clip will not be added.'};
        title='Invalid file type';
        h=msgbox(message,title, 'error');
        waitfor(h);
            if lb1_index < size(get(handles.listbox1, 'String'),1)
                set(handles.listbox1, 'Value', lb1_index+1);
            end
    end
    
    % create new linked files for "standard" filenames (x_x_x.x or x_x.x)
    if ~nsn     % only proceed if std
        std_match=1;
        num_=size(findstr(item, '_'),2);

        if isequal(num_,2) %  xxx_xxx_xxx.xxx
            [test_name, mmm] = strtok(item, '_');
            %Don't add to linked file list if test name is already 'bvqm'
            if strcmp(test_name, 'bvqm')
                add_link=0;
            else
                add_link=1;
                test_name='bvqm';
            end

            [scene_name, mmm] = strtok(mmm, '_');
            [hrc_name_temp, mmm] = strtok(mmm, '_');   % rem should be empty at this point.
            [hrc_name,xxt] = strtok(hrc_name_temp,'.');
        end
        
        if isequal(num_,1)  % xxx_xxx.xxx
            test_name='bvqm';
            [scene_name,mmm] = strtok(item,'_');
            [hrc_name_temp, mmm] = strtok(mmm, '_');
            [hrc_name, xxt] = strtok(hrc_name_temp, '.');
         end
    end % ~nsn
     
    if ~(std_match) % if NOT std filename for avi or yuv file
        % only display this message one time:
        first_time=get(handles.text41, 'Value');

        if isequal(first_time,0)
            message={'Nonstandard clip filename will appear as:';'bvqm_SceneName_HRCName.ext';...
                'in the Selected Files listbox.';'';'Please click OK to enter Scene and HRC names.'};
            title='Nonstandard clip filename';
            h=msgbox(message,title, 'help');
            waitfor(h);
            set(handles.text41, 'Value', 1);
        end
        scene_name='' ; clip_name='';
        test_name='bvqm';

        while (isempty(cell2mat(scene_name)) | isempty(cell2mat(hrc_name)))
            prompt={'Enter a SCENE name:'; 'Enter a HRC name:'};
            title='Enter Scene, HRC';
            
            [ns]= inputdlg(prompt, title, 1, {pre, 'original'});
            if isempty(ns)
                return
            end

            scene_name=ns(1);
            hrc_name=ns(2);
        end

        if iscell(scene_name); scene_name=cell2mat(scene_name);end
        if iscell(hrc_name); hrc_name=cell2mat(hrc_name);end

        % in case 2 or 3 '.'s in the filename, just take the last as
        % the extention
        if ~(ext_yuv | ext_avi)
            if isequal(num_dots,2)
                [p,ext]= strtok(ext,'.');
            end

            if isequal(num_dots,3)
                [p,pp]= strtok(ext,'.');
                [p,ext]= strtok(pp,'.');
            end
        end
    end  % if ~std_match

        %only add file to listbox2 if there is no error
        if noadd ~= 1
            if ~exist('item_new')
                hold_lb2{rows+1} = item;
            else
                hold_lb2{rows+1} = item_new;
            end
            %        next_file_pushbutton7_Callback(hObject, eventdata, handles)
            % advance highlighted item in lb1:
            if lb1_index < size(get(handles.listbox1, 'String'),1)
                set(handles.listbox1, 'Value', lb1_index+1);
            end

            set(handles.listbox2, 'String', hold_lb2);
            set(handles.listbox2, 'Value', length(hold_lb2));

            set(handles.remove_pushbutton3, 'Enable', 'on');

            % are image size and fps data entered?
            check_imagesize(hObject, eventdata, handles);
        end  % noadd ~= 1


%     % Populate the "Auto File Parsing" box if the 1st item in lb2 is parsable.
%     lb2_index = get(handles.listbox2,'Value');
%     if lb2_index == 1
%         prev_file_pushbutton6_Callback(hObject, eventdata, handles);
%     end

    %next_file... fcn will update the video clip data fields...
    % note: nsn=non std name - passed
    % next_file_pushbutton7_Callback(hObject, eventdata, handles)

    set(handles.hrc_edit2, 'String', hrc_name)
    set(handles.scene_edit1, 'String', scene_name)

    if (strcmp(hrc_name, 'original'))
        set(handles.original_clip_checkbox1, 'Enable', 'inactive')
        set(handles.original_clip_checkbox1, 'Value', 1)
        set(handles.hrc_edit2, 'Enable', 'inactive')
        set(handles.scene_edit1, 'Enable', 'inactive');
    else
        set(handles.original_clip_checkbox1, 'Enable', 'inactive')
        set(handles.original_clip_checkbox1, 'Value', 0)
        set(handles.hrc_edit2, 'Enable', 'inactive')
        set(handles.scene_edit1, 'Enable', 'inactive');
    end

    % Store entered scene name, hrc name for later retrieval:
    vcd_fn = fullfile(working_dir, 'bvqm-vcd.mat'); % video clip data filename
    if exist(vcd_fn, 'file')
        load (vcd_fn);
        lb2_index = size(get(handles.listbox2, 'String'),1);
        vcd(lb2_index).scene=scene_name;
        vcd(lb2_index).hrc=hrc_name;
        vcd(lb2_index).filename=item;
        save(vcd_fn, 'vcd');
    else
        lb2_index = size(get(handles.listbox2, 'String'),1);
        vcd(lb2_index).scene=scene_name;
        vcd(lb2_index).hrc=hrc_name;
        vcd(lb2_index).filename=item;
        save(vcd_fn, 'vcd');
    end
fclose(fid);
end  % file loop


% --- Executes on button press in remove_pushbutton3.
function remove_pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to remove_pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

current_selection = get(handles.listbox2, 'Value');
resultsStr = get(handles.listbox2, 'String');
numResults = size(resultsStr,1);
if isempty(resultsStr) ; return; end

resultsStr(current_selection) = [];

% if removing last item, reset data fields
if isequal(numResults,length(current_selection)),
    %if isequal(numResults,current_selection),
    % resultsStr = {'<empty>'};
    current_selection = 1;
    set([handles.remove_pushbutton3], 'Enable', 'off');
    
    %clear test data boxes
    set(handles.image_size_rows_edit, 'String', '');
    set(handles.image_size_cols_edit, 'String', '');
    set(handles.fps_edit, 'String', '');
    set(handles.video_std_popupmenu1, 'Value', 3.0); % 3=prog.
    set(handles.length_text33, 'String', '');

    % clear video clip data boxes
    set(handles.test_edit, 'String', '');
    set(handles.scene_edit1, 'String', '');
    set(handles.hrc_edit2, 'String', '');
    set(handles.original_clip_checkbox1, 'Value', 0);
    set(handles.text6, 'String', ''); % clip file name

    % clear parsing panel:
%     set(handles.auto_parse_uipanel3, 'Visible', 'off')
end

%current_selection = min(current_selection, size(resultsStr,1));
% set(handles.listbox2, 'Value', current_selection-1, 'String', resultsStr)
set(handles.listbox2, 'Value', size(resultsStr,1), 'String', resultsStr)

% remove selected item from vcd list:
    % Store entered scene name, hrc name for later retrieval:
    working_dir=get(handles.text40, 'String');
    vcd_fn = fullfile(working_dir, 'bvqm-vcd.mat');
    load(vcd_fn);
    vcd(current_selection)=[];  % hmmm? but it works...
    save(vcd_fn, 'vcd');
    
    
    % video clip data filename
%      if exist (vcd_fn, 'file')
%         load (vcd_fn);
%         
%         
%         
%         lb2_index = size(get(handles.listbox2, 'String'),1)
%         vcd(lb2_index).scene=scene_name;
%         vcd(lb2_index).hrc=hrc_name;
%         vcd(lb2_index).filename=item
%         save(vcd_fn, 'vcd');
%      else  % first entry in vcd:
%         vcd(1).scene=scene_name;
%         vcd(1).hrc=hrc_name;
%         vcd(1).filename=item
%         save(vcd_fn, 'vcd');
%      end


%guidata(hObject, handles)

%-----------------------------------------------------------------------
% This function will update the test, scene, HRC, image size, and video
% standard fields associated with the clip when the file is double-clicked in listbox2.
%-----------------------------------------------------------------------
% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2

get(handles.bvqm_pc,'SelectionType');
if strcmp(get(handles.bvqm_pc,'SelectionType'),'open')
    lb2_index = get(handles.listbox2,'Value');
    file_list = get(handles.listbox2,'String');
    lb2_filename = file_list{lb2_index};
    [path,name,ext,ver] = fileparts(lb2_filename);

   % Parse filename to update video clip data frame with highlighted file info:
    [test_name, scene_name, hrc_name, ext] = parse_clip(lb2_filename,hObject,eventdata,handles);
    set(handles.text6, 'String', lb2_filename)
end


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%"Advanced Settings" Button - aka "optional settings", aka 'manual settings'
% --- Executes on button press in optional_settings_pushbutton5.
function optional_settings_pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to optional_settings_pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Gather data needed to fill in optional data
lb2_index = get(handles.listbox2,'Value');
file_list = get(handles.listbox2,'String');

%C heck to make sure the user has chosen some video at all
if isempty(file_list)
    warndlg('You must first select video clips.', 'No clips to process');
    return
end
if ispc ; path_sep = '\'; else; path_sep='/'; end
lb2_filename = file_list{lb2_index};
working_dir = get(handles.text40, 'String');
clip_dir = get(handles.clip_directory_text3, 'String');
lb2_pathfn = strcat(clip_dir, path_sep, lb2_filename);




% Don't allow users to alter settings for an 'original' clip:
% [test_name, rem] = strtok(lb2_filename, '_');
% [scene_name, rem] = strtok(rem, '_');
% [hrc_name_temp, rem] = strtok(rem, '_');   % rem should be empty at this point.
% [hrc_name,ext] = strtok(hrc_name_temp,'.');
% 
% if strcmp(hrc_name, 'original')
%     h=warndlg('Can not change advanced settings on original clips.', 'Advanced settings not available');
%     waitfor(h);
%     return
% end



% This statement will prevent errors from occuring when 'options' is
% pressed before files are browsed through:
test_field = get(handles.test_edit, 'String');
if isempty(test_field)
    prev_file_pushbutton6_Callback(hObject, eventdata, handles)
end

% Are we dealing with a paresd file?
[lb2_filename, ext] = strtok(lb2_filename, '@');
if strncmp(ext, '@', 1)
    parsed_file=1;
    file_size = NaN;
else
    parsed_file=0;
end

[fid,message] = fopen(lb2_pathfn);

fseek(fid,0,'eof');
file_size = ftell(fid);

% Pass the bvqm_pc_opt function these inputs:
rows = str2num(get(handles.image_size_rows_edit, 'String'));
cols = str2num(get(handles.image_size_cols_edit, 'String'));

if parsed_file
    if ispc ; path_sep = '\'; else; path_sep='/'; end
    load_file=strcat(working_dir, path_sep, 'bvqm-parsed-', lb2_filename, '.mat');
    load (load_file)
    %     loc_start = parsed(lb2_index).start_frame;
    %     loc_stop = parsed(lb2_index).stop_frame;
    %     num_frames = parsed(lb2_index).stop_frame - parsed(lb2_index).start_frame;

    % This takes care of assigning rows & cols to parsed clips in order to
    % computer the loc_stop.  ??? fill in rows,cols fields??
    [f_name,f_type]=strtok(lb2_filename, '.');
    if strcmp(f_type, '.yuv')
        rows=486 ; cols=720;
    end
    if strcmp(f_type, '.avi')
        avi_info=aviinfo(lb2_filename);
        rows=avi_info.Height;  cols=avi_info.Width;
    end
end

[f_name,f_type]=strtok(lb2_filename, '.');
loc_start = 1;
 
if strcmp(f_type, '.yuv')
    num_frames  = floor(file_size/(rows*cols*2));  % XXX - NOT ALWAYS ACCURATE)!
    loc_stop = num_frames;
end

if strcmp(f_type, '.avi')
%   avi_info = aviinfo(lb2_filename);
   avi_info = aviinfo(lb2_pathfn);
   num_frames=avi_info.NumFrames;
   loc_stop=num_frames;
end

fclose(fid);

working_dir = get(handles.text40, 'String');
[varargout] = bvqm_pc_opt(lb2_index, lb2_filename, num_frames, loc_start, loc_stop, rows, cols, parsed_file, working_dir);


% --- Executes on selection change in video_std_popupmenu1.
function video_std_popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to video_std_popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns video_std_popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from video_std_popupmenu1

% query to set default image size values based upon video std:
video_std_idx=get(handles.video_std_popupmenu1, 'Value');
video_std_lst=get(handles.video_std_popupmenu1, 'String');
video_std=video_std_lst(video_std_idx);





% --- Executes during object creation, after setting all properties.
function video_std_popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to video_std_popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in original_clip_checkbox1.
function original_clip_checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to original_clip_checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of original_clip_checkbox1


% --- Executes on button press in prev_file_pushbutton6.
function prev_file_pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to prev_file_pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Update highlighted selection when 'prev' button is pushed
lb2_filelist = get(handles.listbox2,'String');
lb2_filelist_size = size(lb2_filelist,1);
lb2_index = get(handles.listbox2, 'Value');

% save whatever is entered in rows, cols, fps boxes:
%  lb2_filename = lb2_filelist{lb2_index};
%  save_tdf(lb2_filename, hObject, eventdata, handles);

if lb2_index == 1
    lb2_index = lb2_index;
else
    lb2_index = lb2_index - 1;
end

if lb2_filelist_size==0; return; end   % in case no clips selected yet...

lb2_filename = lb2_filelist{lb2_index};
set(handles.listbox2, 'Value', lb2_index)
set(handles.text6, 'String', lb2_filename)

% Make sure image size, fps is filled in; if not, prompt for values:
fps=get(handles.fps_edit, 'String');
cols=get(handles.image_size_cols_edit, 'String');
rows=get(handles.image_size_rows_edit, 'String');
% if (strcmp(fps,'') | strcmp(cols,'') | strcmp(rows,''))


% can assume image size data is already valid (checked in 'check_imagesize' fcn
% if (isempty(fps) | isempty(cols) | isempty(rows))
%     check_imagesize(hObject, eventdata, handles);
% end

% Now, parse filename to update video clip data frame with highlighted file info:
[test_name, scene_name, hrc_name, ext] = parse_clip(lb2_filename,hObject, eventdata,handles);

% --- Executes on button press in next_file_pushbutton7.
function next_file_pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to next_file_pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% Update highlighted selection when 'next' button is pushed
lb2_filelist = get(handles.listbox2,'String');
lb2_filelist_size = size(lb2_filelist,1);
lb2_index = get(handles.listbox2, 'Value');

% save whatever is entered in rows, cols, fps boxes:
%  lb2_filename = lb2_filelist{lb2_index};
%  save_tdf(lb2_filename, hObject, eventdata, handles);


if lb2_index == lb2_filelist_size
    lb2_index = lb2_index;
else
    lb2_index = lb2_index + 1;
end

if lb2_filelist_size==0; return; end   % in case no clips selected yet...

lb2_filename = lb2_filelist{lb2_index};
set(handles.listbox2, 'Value', lb2_index);
set(handles.text6, 'String', lb2_filename);

% Make sure image size, fps is filled in; if not, prompt for values:
fps=get(handles.fps_edit, 'String');
cols=get(handles.image_size_cols_edit, 'String');
rows=get(handles.image_size_rows_edit, 'String');
if (strcmp(fps,'') | strcmp(cols,'') | strcmp(rows,''))
    check_imagesize(hObject, eventdata, handles);
end




% Parse filename to update video clip data frame with highlighted file info:
 [test_name, scene_name, hrc_name, ext] = parse_clip(lb2_filename,hObject,eventdata,handles);

% --- Executes on button press in parse_clip_pushbutton.
function parse_clip_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to parse_clip_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of parse_clip_pushbutton

% Get file info...
lb2_index = get(handles.listbox2, 'Value');
file_list=get(handles.listbox2,'String');
filename = file_list(lb2_index);
[test_name, rem] = strtok(filename, '_');
[scene_name, rem] = strtok(rem, '_');
[hrc_name_temp, rem] = strtok(rem, '_');  % rem should be empty at this point.
[hrc_name,ext] = strtok(hrc_name_temp,'.');


% Parse clip and create parsed clip filenames
clip_length = str2num(strtok(get(handles.length_text33, 'String'), ' '));
parse_length = str2num(get(handles.parse_cl_edit3, 'String'));
parse_shift = str2num(get(handles.parse_ts_edit4, 'String'));
fps = str2num(get(handles.fps_edit, 'String'));

start_time = 0 ; stop_time = 0;
start_frame = 1 ; stop_frame = 0;
parse_count = 1;
while stop_time <= clip_length
    stop_time = start_time + parse_length;
%    parsed_name(parse_count) = strcat(test_name, '_', scene_name, num2str(stop_time), '_', hrc_name, ext);
     parsed_name(parse_count) = strcat(filename, '@', num2str(stop_time));
    stop_frame = ceil(start_time + (parse_count * parse_length * fps));

    % create 'parsed' structure to track parsed files:
    
    % parsed(parse_count).filename = parsed_name(parse_count);
    parsed(parse_count).filename = filename;
    parsed(parse_count).parsed_name = parsed_name(parse_count);
    parsed(parse_count).clip_length = parse_length;
    parsed(parse_count).time_shift  = parse_shift;
    parsed(parse_count).start_time = start_time;
    parsed(parse_count).stop_time = stop_time;
    parsed(parse_count).start_frame = start_frame;
    parsed(parse_count).stop_frame = stop_frame;

    start_time = start_time + parse_shift;
    start_frame = stop_frame + 1;
    parse_count = parse_count +1;
end
set(handles.auto_parse_uipanel3, 'Visible', 'off')


% store files that appear AFTER the file which will be parsed so it can be redisplayed:
hold_lb2 = get(handles.listbox2,'String');
num_postfiles = size(hold_lb2,1) - lb2_index;
for postfile_count = 1:num_postfiles
    stored_filename{postfile_count} = file_list(lb2_index +postfile_count);
    postfile_count = postfile_count +1;
end  % now, files that appeared after parsed file are in stored in stored_filename{}

% update LB2 with parsed clips:
file_count = 1;

while file_count < parse_count
%    displayed_filename(file_count) = strcat(parsed_name(file_count),'@');
    displayed_filename(file_count) = parsed_name(file_count);
    % note: '-1' will prevent original file name from displaying in list.
    hold_lb2(lb2_index-1 + file_count) = displayed_filename(file_count);
    set(handles.listbox2, 'String', hold_lb2);
    hold_lb2 = get(handles.listbox2, 'String');

    file_count = file_count +1;
end

% add back the stored_filenames that came after the parsed file:
hold_lb2 = get(handles.listbox2, 'String');
maxx = size(hold_lb2,1);
for reload_count = 1:postfile_count-1
    hold_lb2(maxx +reload_count) = stored_filename{reload_count};
    set(handles.listbox2, 'String', hold_lb2);
end

working_dir = get(handles.text40, 'String');
if ispc ; path_sep = '\'; else; path_sep='/'; end
save_file=strcat(working_dir, path_sep, 'bvqm-parsed-',filename, '.mat');
%save_file=strcat(working_dir, path_sep, 'bvqm-parsed');
save_file = cell2mat(save_file);
save (save_file, 'parsed');

function parse_ts_edit4_Callback(hObject, eventdata, handles)
% hObject    handle to parse_ts_edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of parse_ts_edit4 as text
%        str2double(get(hObject,'String')) returns contents of parse_ts_edit4 as a double


% --- Executes during object creation, after setting all properties.
function parse_ts_edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parse_ts_edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function parse_cl_edit3_Callback(hObject, eventdata, handles)
% hObject    handle to parse_cl_edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of parse_cl_edit3 as text
%        str2double(get(hObject,'String')) returns contents of parse_cl_edit3 as a double


% --- Executes during object creation, after setting all properties.
function parse_cl_edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parse_cl_edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function image_size_rows_edit_Callback(hObject, eventdata, handles)
% hObject    handle to image_size_rows_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of image_size_rows_edit as text
%        str2double(get(hObject,'String')) returns contents of image_size_rows_edit as a double


% --- Executes during object creation, after setting all properties.
function image_size_rows_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_size_rows_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in select_all_pushbutton11.
function select_all_pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to select_all_pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% New Method: loop through each file, calling add_file fcn for each file
file_list = get(handles.listbox1, 'String');
num_files = size(file_list,1);
set(handles.listbox1,'Value', 1);

for file_count = 1:num_files
    add_pushbutton2_Callback(hObject, eventdata, handles);
end



% % Add all files from listbox1 to listbox2. do NOT include directories.
% file_list = get(handles.listbox1, 'String');
% num_files = size(file_list,1);
% for file_count = 1:num_files
%     item=(file_list(file_count));
%     mat_item = cell2mat(item);
%     [mat_item_prefix, mat_item_ext] = strtok(mat_item,'.');
%     
%     % if item is a directory, skip it; if not, add it to LB2
%     if (~isdir(mat_item))
%         % if item is not an '.avi' or '.yuv' file, then skip it
%         if (strcmp(mat_item_ext, '.yuv')) || (strcmp(mat_item_ext, '.avi'))
%         
%         % place item in listbox2 if NOT a directory
%         hold_lb2 = get(handles.listbox2, 'String');
%         [rows,cols]=size(hold_lb2);
%         hold_lb2{rows+1} = mat_item;
%         set(handles.listbox2, 'String', hold_lb2);
%         hold_lb2 = get(handles.listbox2, 'String');
%         hold_lb2 = unique(hold_lb2); %?XXXX
%         set([handles.remove_pushbutton3], 'Enable', 'on');
%         end
%         
%     end
% end
% 
% % Make sure image size, fps is filled in; if not, prompt for values:
% fps=get(handles.fps_edit, 'String');
% cols=get(handles.image_size_cols_edit, 'String');
% rows=get(handles.image_size_rows_edit, 'String');
% 
% %if (strcmp(fps,'') | strcmp(cols,'') | strcmp(rows,''))
% if (isempty(fps) | isempty(cols) | isempty(rows))
%     check_imagesize(hObject, eventdata, handles);
% end
% 
% % Populate the "Auto File Parsing" box.
% lb2_index = get(handles.listbox2,'Value');
% if lb2_index == 1
%         prev_file_pushbutton6_Callback(hObject, eventdata, handles);
% end
% %guidata(hObject, handles)
% 

% --- Executes on button press in remove_all_pushbutton12.
function remove_all_pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to remove_all_pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

working_dir=get(handles.text40, 'String');

% clear lb2, test data, and vide clip data
set(handles.listbox2, 'String', '');
set(handles.listbox1, 'Value', 1);

%clear test data boxes
set(handles.image_size_rows_edit, 'String', '');
set(handles.image_size_cols_edit, 'String', '');
set(handles.fps_edit, 'String', '');
set(handles.video_std_popupmenu1, 'Value', 3.0); % 3=prog.
set(handles.length_text33, 'String', '');

% clear video clip data boxes
set(handles.test_edit, 'String', '');
set(handles.scene_edit1, 'String', '');
set(handles.hrc_edit2, 'String', '');
set(handles.original_clip_checkbox1, 'Value', 0);
set(handles.text6, 'String', ''); % clip file name

% clear parsing panel:
% set(handles.auto_parse_uipanel3, 'Visible', 'off')

% TEMP FIX:
    warning off all; % no usefull info contained therein

% delete vcd file
    delete (fullfile(working_dir, 'bvqm-vcd.mat')); % video clip data filename
    delete('bvqm-ptrunk.mat');

% -- NO LONGER BEING USED.... but hold on to    
%     try load (fullfile(working_dir, 'bvqm-linked.mat'))
%         num_linked=size(linked_file,2);
%         count=1;
%         while count <= num_linked
%             try delete (linked_file{count}) ; catch; err=1 ; end
%             count=count+1;
%         end
%         delete(fullfile(working_dir, 'bvqm-linked.mat'));
%     catch
%         nosuchfile=1;
%     end

    
% --- This function will parse file names and update the Test Data & Video Clip Data info boxes.
% --- If the filename is not standard format, no parsing will occur.
function [test_name, scene_name, hrc_name, ext, rows, cols,fps] = parse_clip(lb2_filename,hObject, eventdata, handles)

% 03Jan06: change so that image size, fps is not automatically set by the
% file type.  It will be set by the video standard selection or user
% entered data.
set(handles.scene_edit1, 'Visible', 'on')
set(handles.hrc_edit2, 'Visible', 'on')
set(handles.original_clip_checkbox1, 'Visible', 'on')


% info needed to handle rows, cols, fps edit boxes:
if ispc ; path_sep = '\'; else; path_sep='/'; end
working_dir = get(handles.text40, 'String');
test_data_file = strcat(working_dir, path_sep, 'bvqm-test_data.mat');
%clip_dir = cell2mat(get(handles.clip_directory_text3, 'String'));
clip_dir = get(handles.clip_directory_text3, 'String');

%lb2_pathfn = strcat(clip_dir, path_sep, lb2_filename);
lb2_pathfn = fullfile(clip_dir,lb2_filename);
trunc_save_file=strcat(working_dir, path_sep, 'bvqm-ptrunk.mat');

% Modify lb2_filename if it is a parsed clip (to remove '@*' suffix:)
% Also, inactivate the auto file parsing box 
[temp_lb2_filename, rem] = strtok(lb2_filename, '@');
%[lb2_filename, rem] = strtok(lb2_filename, '@');

if strncmp(rem, '@', 1)
   set(handles.auto_parse_uipanel3, 'Visible', 'off')
 
 %   lb2_filename = temp_lb2_filename;
%          set(handles.parse_cl_edit3, 'Enable', 'off');
%          set(handles.parse_ts_edit4, 'Enable', 'off');
%          set(handles.parse_clip_pushbutton, 'Enable', 'off');
end

% Check for standard filename format
std_yuv_match = regexp(lb2_filename, '\w*_\w*_\w*.yuv\w*');
if isempty(std_yuv_match); std_yuv_match=0; end

std_avi_match = regexp(lb2_filename, '\w*_\w*_\w*.avi\w*');
if isempty(std_avi_match); std_avi_match=0; end

% TEMP FIX:  KNOWN=1 b/c was set when file was added.
%std_yuv_match=1; stc_avi_match=1;

% for names like *_*_*.*
if  ( (std_yuv_match) | (std_avi_match) ) % if std filename for avi or yuv file
    % Parse filename to update video clip data frame with highlighted file info:
    [test_name, rem] = strtok(lb2_filename, '_');
    [scene_name, rem] = strtok(rem, '_');
    [hrc_name_temp, rem] = strtok(rem, '_');   % rem should be empty at this point.
    [hrc_name,ext] = strtok(hrc_name_temp,'.');

else  % file name form : *_*.*
    test_name='bvqm';
    [scene_name, rem] = strtok(lb2_filename, '_');
    [hrc_name_temp, rem] = strtok(rem, '_');   % rem should be empty at this point.
    [hrc_name,ext] = strtok(hrc_name_temp,'.');

    set(handles.test_edit, 'String', test_name); 
    set(handles.test_edit, 'Enable', 'off');
    set(handles.scene_edit1, 'String', scene_name);
    set(handles.scene_edit1, 'Enable', 'inactive');
    set(handles.hrc_edit2, 'String', hrc_name);
    set(handles.hrc_edit2, 'Enable', 'inactive');

    if (strcmp(hrc_name, 'original'))
        set(handles.original_clip_checkbox1, 'Enable', 'inactive')
        set(handles.original_clip_checkbox1, 'Value', 1)
        set(handles.hrc_edit2, 'Enable', 'inactive')
    else
        set(handles.original_clip_checkbox1, 'Enable', 'inactive')
        set(handles.original_clip_checkbox1, 'Value', 0)
        set(handles.hrc_edit2, 'Enable', 'inactive')
    end
%else % if not standard filename format
    % nothing now...
end

%     set(handles.original_clip_checkbox1, 'Enable', 'on')
%     set(handles.original_clip_checkbox1, 'Value', 0)
%     set(handles.hrc_edit2, 'Enable', 'on')
%     set(handles.length_text33, 'String', '')
%     set(handles.auto_parse_text34, 'Enable', 'off')
%     set(handles.fps_edit, 'Enable', 'on')
%     set(handles.image_size_rows_edit, 'Enable', 'on')
%     set(handles.image_size_cols_edit, 'Enable', 'on')
     [prefix,ext] = strtok(lb2_filename,'.');
    % set output vars to '' if not std file name
%    scene_name=''; test_name=''; hrc_name='';

    %     set(handles.auto_parse_uipanel3, 'Visible', 'off')
    %     set(handles.hrc_edit2, 'String', '')
    %     set(handles.scene_edit1, 'String', '')
    %     set(handles.original_clip_checkbox1, 'Enable', 'on')
    %     set(handles.original_clip_checkbox1, 'Value', 0)
    %     set(handles.hrc_edit2, 'Enable', 'on')
    %     set(handles.length_text33, 'String', '')
    %     set(handles.auto_parse_text34, 'Enable', 'off')
    %     set(handles.fps_edit, 'Enable', 'on')
    %     set(handles.image_size_rows_edit, 'Enable', 'on')
    %     set(handles.image_size_cols_edit, 'Enable', 'on')
    %     [prefix,ext] = strtok(lb2_filename,'.');
    %     % set output vars to '' if not std file name
    %     scene_name=''; test_name=''; hrc_name='';
%end

if strcmp(ext, '.avi') % unparsed avi file (i.e. <15 secs)
%     avi_info=aviinfo(lb2_filename);
%     % load saved rows, cols, fps data:
%     if exist (test_data_file, 'file')
%         load (test_data_file);
%         set(handles.image_size_rows_edit, 'String', tdf_rows)
%         set(handles.image_size_cols_edit, 'String', tdf_cols)
%         set(handles.fps_edit, 'String', tdf_fps)
%     else
%         %if tdf doesn't exist, set to appropriate data from MATLAB fcn aviinfo
%         set(handles.image_size_rows_edit, 'String', avi_info.Height)
%         set(handles.image_size_rows_edit, 'Enable', 'on')
%         set(handles.image_size_cols_edit, 'String', avi_info.Width)
%         set(handles.image_size_cols_edit, 'Enable', 'on')
%         % NOTE: FramesPerSecond may return wrong value- set fps to aviinfo.Rate
%         set(handles.fps_edit, 'String', avi_info.Rate)
%         set(handles.fps_edit, 'Enable', 'on')
%     end
%     % save whatever is entered in test data fields:
%     save_tdf(hObject, eventdata,handles);


    % Compute length (in seconds) of the clip:
    avi_info=aviinfo(lb2_pathfn);

    length=avi_info.NumFrames/avi_info.FramesPerSecond;
    
    %length=avi_info.NumFrames/avi_info.Rate; 
    tlength = sprintf('%3.1f seconds', length);
    set(handles.length_text33, 'String', tlength);
    if length >= 15 % if parsing required
%         set(handles.auto_parse_text34, 'Enable', 'on')
%         set(handles.auto_parse_uipanel3, 'Visible', 'on')
%         set(handles.parse_clip_pushbutton, 'Visible', 'on')
%         set(handles.scene_edit1, 'Visible', 'off')
%         set(handles.hrc_edit2, 'Visible', 'off')
%         set(handles.original_clip_checkbox1, 'Visible', 'off')
% % %      TEMP FIX: files longer than 15 secs get truncated:
% %             if exist(trunc_save_file, 'file')
% %                 load(trunc_save_file);
% %                 tc=size(trunc_save_file,1) +1;
% %             else
% %                 tc=1;
% %             end
% % 
% %            % NOTE: need to display below values when manual settings are called and when jclips is written.     
% %             trunc_filename(tc)={lb2_filename};
% %             trunc_stop_frame(tc)=fps*15;
% %             save (trunc_save_file, 'trunc_filename', 'trunc_stop_frame')
% %         
    else


    end
elseif strcmp(ext, '.avi@*') % parsed avi file
    % load saved rows, cols, fps data:
     if exist (test_data_file, 'file')
%          load (test_data_file);
%          set(handles.image_size_rows_edit, 'String', tdf_rows)
%          set(handles.image_size_cols_edit, 'String', tdf_cols)
%          set(handles.fps_edit, 'String', tdf_fps)
    else
        %if tdf doesn't exist, set to default:
        set(handles.image_size_rows_edit, 'String', '')
        set(handles.image_size_rows_edit, 'Enable', 'on')
        set(handles.image_size_cols_edit, 'String', '')
        set(handles.image_size_cols_edit, 'Enable', 'on')
        set(handles.fps_edit, 'String', '')
        set(handles.fps_edit, 'Enable', 'on')
        set(handles.length_text33, 'String', 'parsed clip')
        set(handles.auto_parse_uipanel3, 'Visible', 'off')
    end

    % save whatever is entered in test data fields:
%     save_tdf(hObject, eventdata,handles);

end

if strcmp(ext,'.yuv') % unparsed yuv file (i.e. <15 secs)
    % load saved rows, cols, fps data:
%     if exist (test_data_file, 'file')
%         load (test_data_file);
%         set(handles.image_size_rows_edit, 'String', tdf_rows)
%         set(handles.image_size_cols_edit, 'String', tdf_cols)
%         set(handles.fps_edit, 'String', tdf_fps)
%     else
%         %if tdf doesn't exist, set to default:
%         set(handles.image_size_rows_edit, 'String', 486)
%         set(handles.image_size_rows_edit, 'Enable', 'on')
%         set(handles.image_size_cols_edit, 'String', 720)
%         set(handles.image_size_cols_edit, 'Enable', 'on')
%         set(handles.fps_edit, 'String', 29.97)
%         set(handles.fps_edit, 'Enable', 'on')
%     end
% 
%      % save whatever is entered in test data fields:
%      save_tdf(hObject, eventdata,handles);
    
    % Compute length (in seconds) of the clip, determine if parsing req'd:
    % skip if invalid fid - allows for exit and sa ???????????
    
    [fid,message] = fopen(lb2_pathfn);   % PPPPP path

%    if fid ~= -1
        fseek(fid,0,'eof');
        file_size = ftell(fid);


        rows = str2num(get(handles.image_size_rows_edit, 'String'));
        cols = str2num(get(handles.image_size_cols_edit, 'String'));
%         num_frames = round(file_size/(rows*cols*2));
        num_frames = file_size/(rows*cols*2);
        num_frames = floor(num_frames);
        fps = str2num(get(handles.fps_edit, 'String'));
        length = (num_frames/fps);
        tlength = sprintf('%3.1f seconds', num_frames/fps);
        set(handles.length_text33, 'String', tlength);
        if length >= 15 % if parsing required
%             set(handles.auto_parse_text34, 'Enable', 'on')
%             set(handles.auto_parse_uipanel3, 'Visible', 'on')
%             set(handles.parse_clip_pushbutton, 'Visible', 'on')
%             set(handles.scene_edit1, 'Visible', 'off')
%             set(handles.hrc_edit2, 'Visible', 'off')
%             set(handles.original_clip_checkbox1, 'Visible', 'off')
%      TEMP FIX: files longer than 15 secs get truncated:
% %             if exist(trunc_save_file, 'file')
% %                 load(trunc_save_file);
% %                 tc=size(trunc_save_file,1) +1;
% %             else
% %                 tc=1;
% %             end
% % 
% %            % NOTE: need to display below values when manual settings are called and when jclips is written.     
% %             trunc_filename(tc)={lb2_filename};
% %             trunc_stop_frame(tc)=fps*15;
% %             save (trunc_save_file, 'trunc_filename', 'trunc_stop_frame')
        else
%             set(handles.auto_parse_text34, 'Enable', 'off')
%             set(handles.auto_parse_uipanel3, 'Visible', 'off')
%             set(handles.parse_clip_pushbutton, 'Visible', 'off')
        end
        fclose(fid);
    elseif strcmp(ext, '.yuv@*') % parsed yuv file
        % load saved rows, cols, fps data:
        %     if exist (test_data_file, 'file')
        %         load (test_data_file);
        %         set(handles.image_size_rows_edit, 'String', tdf_rows)
        %         set(handles.image_size_cols_edit, 'String', tdf_cols)
        %         set(handles.fps_edit, 'String', tdf_fps)
        %     else
        %         %if tdf doesn't exist, set to default:
        %         set(handles.image_size_rows_edit, 'String', 486)
        %         set(handles.image_size_rows_edit, 'Enable', 'on')
        %         set(handles.image_size_cols_edit, 'String', 720)
        %         set(handles.image_size_cols_edit, 'Enable', 'on')
        %         set(handles.fps_edit, 'String', 29.97)
        %         set(handles.fps_edit, 'Enable', 'on')
        %         set(handles.length_text33, 'String', 'parsed clip')
        %         set(handles.auto_parse_uipanel3, 'Visible', 'off')
        %     end
        %
        %     % save whatever is entered in test data fields:
        %     save_tdf(hObject, eventdata,handles);

    end

    % if not a yuv or avi file, clear fields
    % 19June06 : this taken care of when file is added: 
    % allow ony .yuv or .avi files to be added to list.
    
%     if ~(strcmp(ext, '.yuv') | strcmp(ext, '.yuv@*') | strcmp(ext, '.avi') | strcmp(ext, '.avi@*'))
%         set(handles.image_size_rows_edit, 'String', '')
%         set(handles.image_size_cols_edit, 'String', '')
%         set(handles.fps_edit, 'String', '')
%         set(handles.length_text33, 'String', '');
%     end
%end

% populate video clip data box with info from bvqm-vcd.mat
vcd_fn=fullfile(working_dir, 'bvqm-vcd.mat');
lb2_index = get(handles.listbox2,'Value');


if exist (vcd_fn, 'file')
    load (vcd_fn);
    scene_name=vcd(lb2_index).scene;
    hrc_name=vcd(lb2_index).hrc;
    filename=vcd(lb2_index).filename;
       
    set(handles.scene_edit1, 'String', scene_name);
    set(handles.scene_edit1, 'Enable', 'inactive');
    set(handles.hrc_edit2, 'String', hrc_name);
    set(handles.hrc_edit2, 'Enable', 'inactive');

    if (strcmp(hrc_name, 'original'))
        set(handles.original_clip_checkbox1, 'Enable', 'inactive')
        set(handles.original_clip_checkbox1, 'Value', 1)
        set(handles.hrc_edit2, 'Enable', 'inactive')
    else
        set(handles.original_clip_checkbox1, 'Enable', 'inactive')
        set(handles.original_clip_checkbox1, 'Value', 0)
        set(handles.hrc_edit2, 'Enable', 'inactive')
    end
end

% --- Executes on button press in continue_pushbutton13.
function continue_pushbutton13_Callback(hObject, eventdata, handles, exit_pb)

% hObject    handle to continue_pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% check for exit_status variable (exit_pb); assign to 0 if nonexistent.
% NOTE: this variable is set when the 'EXIT' pb is pressed.
 ese=exist('exit_pb');
 if ~ese; exit_pb=0; end
lb2_filelist=get(handles.listbox2, 'String');
if isempty(lb2_filelist)
    warndlg('You must first select video clips.', 'No clips to process');
    return
end

% Initalize yuv file size error display marker
yuv_sz_err_dsp=0;
 
% % first, make sure image size, fps is filled in:
% check_imagesize(hObject, eventdata, handles);
%   no nlonger necessary - done at the beginning.

% If files not browsed through yet, force a browse in order to populate the
% test data boxes -- otherwise, errors occur...
test_field = get(handles.test_edit, 'String');
if isempty(test_field)
    prev_file_pushbutton6_Callback(hObject, eventdata, handles)
end


parsed_index =1;
lb2_index =1;
lb2_size=size(lb2_filelist,1);
% info needed to save entered rows, cols, fps data
if ispc ; path_sep = '\'; else; path_sep='/'; end
working_dir = get(handles.text40, 'String');
clip_dir = get(handles.clip_directory_text3, 'String');

        

%Loop through each item in listbox2 to update GClips structure:
% Assign clips default values, unless 'Optional Settings' were changed
while lb2_index <= lb2_size
    
    lb2_filename = lb2_filelist{lb2_index};
    test_data_file = strcat(working_dir, path_sep, 'bvqm-test_data.mat');

    % Was optional data already entered? Look for saved options .mat file:
    [fn_pre, fn_ext] = strtok(lb2_filename, '.');
    [fn_post] = strtok(fn_ext, '@');

    if ispc ; path_sep = '\'; else; path_sep='/'; end
    opt_filename = strcat(working_dir, path_sep, 'bvqm-opt_', fn_pre,  fn_post, '.mat');

    % set file_size, rows, cols variables:
    tmp_lb2_filename=strtok(lb2_filename, '@');
    open_f=strcat(clip_dir, path_sep, tmp_lb2_filename);
    [fid,message] = fopen(open_f);
    if fid ~= -1
        fseek(fid,0,'eof');
        file_size = ftell(fid);
    end
    rows= str2num(get(handles.image_size_rows_edit, 'String'));
    cols= str2num(get(handles.image_size_cols_edit, 'String'));

    if  exist(opt_filename, 'file') % Load Gclips optional data from optional settings
        load(opt_filename);
        [dum,opt_index] = size(oclips); % need to ref the last column...
        popupmenu1_string_list=get(handles.video_std_popupmenu1, 'String');
        [test_name, scene_name, hrc_name, ext] = parse_clip(lb2_filename,hObject,eventdata,handles);

        jclips(lb2_index).test  = {test_name};
        jclips(lb2_index).scene = {scene_name};
        jclips(lb2_index).hrc   = {hrc_name};

        % logic here to save user entered rows, cols, fps...
        jclips(lb2_index).image_size.rows = str2num(get(handles.image_size_rows_edit, 'String'));
        jclips(lb2_index).image_size.cols = str2num(get(handles.image_size_cols_edit, 'String'));
        jclips(lb2_index).spatial.horizontal = oclips(opt_index).spatial.horizontal;
        jclips(lb2_index).spatial.vertical   = oclips(opt_index).spatial.vertical;

        jclips(lb2_index).luminance_gain   = oclips(opt_index).luminance_gain;
        jclips(lb2_index).luminance_offset = oclips(opt_index).luminance_offset;
        jclips(lb2_index).video_standard = popupmenu1_string_list{get(handles.video_std_popupmenu1, 'Value')};
        jclips(lb2_index).fps            = str2num(get(handles.fps_edit, 'String'));
        jclips(lb2_index).subj_system    = {''};  % XXX don't know!

        jclips(lb2_index).mos = oclips(opt_index).mos;
        if isempty(jclips(lb2_index).mos)
            jclips(lb2_index).mos = NaN;
        end

        jclips(lb2_index).stdev = oclips(opt_index).stdev;
        if isempty(jclips(lb2_index).stdev);
            jclips(lb2_index).stdev = NaN;
        end

        jclips(lb2_index).inlsa_mos = oclips(opt_index).inlsa_mos;
        if isempty(jclips(lb2_index).inlsa_mos)
            jclips(lb2_index).inlsa_mos = NaN;
        end

        [jclips_filename, rem] = strtok(lb2_filename, '@');
        jclips(lb2_index).file_name     = {jclips_filename};


%         %Are we dealing with a parsed file?  If so, set the correct start/stop
%         %frame and align start/stop
%         parsed_mat = strcat(working_dir,path_sep, 'bvqm-parsed-', jclips_filename, '.mat');
%         if exist (parsed_mat, 'file')
%             load (parsed_mat);
%             jclips(lb2_index).loc_start = parsed(parsed_index).start_frame;
%             jclips(lb2_index).align_start = jclips(lb2_index).loc_start;
%             jclips(lb2_index).loc_stop  = parsed(parsed_index).stop_frame;
%             jclips(lb2_index).align_stop = jclips(lb2_index).loc_stop;
%             parsed_index = parsed_index + 1;
%         else % not a parsed file -set loc_start,stop
             jclips(lb2_index).loc_start = str2num(mat2str(cell2mat(oclips(opt_index).loc_start)));
             jclips(lb2_index).loc_stop  = str2num(mat2str(cell2mat(oclips(opt_index).loc_stop)));
%         end
 
% NO - leave entered values  05JULY06
        % Set align_start/stop to NaN- 26May06-
%         jclips(lb2_index).align_start = NaN;
%         jclips(lb2_index).align_stop = NaN;

                 jclips(lb2_index).align_start = str2num(oclips(opt_index).align_start);
                if isempty(jclips(lb2_index).align_start)
                    jclips(lb2_index).align_start = NaN;
                end
        %
                 jclips(lb2_index).align_stop  = str2num(oclips(opt_index).align_stop);
                if isempty(jclips(lb2_index).align_stop);
                    jclips(lb2_index).align_stop = NaN;
                end

        jclips(lb2_index).cvr = oclips(opt_index).cvr;

        %         if ischar(jclips(lb2_index).cvr.top)
        %             jclips(lb2_index).cvr.top = str2num(jclips(lb2_index).cvr.top);
        %         end
        %
        %         if iscellstr(jclips(lb2_index).cvr.bottom)
        %             jclips(lb2_index).cvr.bottom = str2num(cell2mat(jclips(lb2_index).cvr.bottom));
        %         end
        %         if ischar(jclips(lb2_index).cvr.left)
        %             jclips(lb2_index).cvr.left = str2num(jclips(lb2_index).cvr.left);
        %         end
        %
        %         if iscellstr(jclips(lb2_index).cvr.right)
        %             jclips(lb2_index).cvr.right = str2num(cell2mat(jclips(lb2_index).cvr.right));
        %         end

        jclips(lb2_index).viewers     = oclips(opt_index).viewers;

        jclips(lb2_index).hrc_definition   = {oclips(opt_index).hrc_definition};
        jclips(lb2_index).scene_definition = {oclips(opt_index).scene_definition};
        jclips(lb2_index).subj_system      = {oclips(opt_index).subj_system};

        jclips(lb2_index).scale.horizontal = oclips(opt_index).scale.horizontal;
        jclips(lb2_index).scale.vertical  = oclips(opt_index).scale.vertical;
        
        % set num_frames - will be needed later...
        num_frames = file_size/(rows*cols*2);
    else
        % no "advanced" or "optional" settings were inputted, so set defaults here...
%         tmp_lb2_filename=strtok(lb2_filename, '@');
%         open_f=strcat(clip_dir, path_sep, tmp_lb2_filename);
%         [fid,message] = fopen(open_f);

        if fid ~= -1
%             fseek(fid,0,'eof');
%             file_size = ftell(fid);  % set at beginning of loop...
            % set rows, cols values for yuv files:
            if strncmp(fn_ext, '.yuv', 4)
                % load saved rows, cols, fps data:
                if exist (test_data_file, 'file')
                    % load (test_data_file);
                    % rows=str2num(tdf_rows); cols=str2num(tdf_cols); fps =str2num(tdf_fps);
                else
                    %if tdf doesn't exist, set to default
%                     rows= str2num(get(handles.image_size_rows_edit, 'String'));
%                     cols= str2num(get(handles.image_size_cols_edit, 'String'));
                    fps = str2num(get(handles.fps_edit, 'String'));
                    num_frames = file_size/(rows*cols*2);
                    num_frames=floor(num_frames); 
                end
            end

            % set rows, cols values for avi files:
            if strncmp(fn_ext, '.avi*', 4)
                % load saved rows, cols, fps data:
                if exist (test_data_file, 'file')
                    % load (test_data_file);
                    % rows=str2num(tdf_rows); cols=str2num(tdf_cols); fps =str2num(tdf_fps);
                else
                    %if tdf doesn't exist, set to default
                    avi_info=aviinfo(fullfile(clip_dir,tmp_lb2_filename));
                    rows= avi_info.Height; cols = avi_info.Width;
                    num_frames=avi_info.NumFrames;
                end
            end
        else
            file_size = NaN;
            rows = NaN;
            cols = NaN;
            num_frames = NaN;
        end

        popupmenu1_string_list=get(handles.video_std_popupmenu1, 'String');
        [test_name, scene_name, hrc_name, ext] = parse_clip(lb2_filename,hObject,eventdata,handles);

        jclips(lb2_index).test  = {test_name};
        jclips(lb2_index).scene = {scene_name};
        jclips(lb2_index).hrc   = {hrc_name};
        %jclips(lb2_index).image_size.rows = str2num(get(handles.image_size_rows_edit, 'String'));
        %jlips(lb2_index).image_size.cols = str2num(get(handles.image_size_cols_edit, 'String'));
        jclips(lb2_index).image_size.rows = rows;
        jclips(lb2_index).image_size.cols = cols;

        jclips(lb2_index).spatial.horizontal = 0;
        jclips(lb2_index).spatial.vertical   = 0;
        jclips(lb2_index).luminance_gain   = 1.000;
        jclips(lb2_index).luminance_offset = 0.0000;
        jclips(lb2_index).video_standard = popupmenu1_string_list{get(handles.video_std_popupmenu1, 'Value')};
        jclips(lb2_index).fps            = str2num(get(handles.fps_edit, 'String'));
        jclips(lb2_index).subj_system   = {''};  % XXX don't know!

        jclips(lb2_index).mos         = NaN;
        jclips(lb2_index).stdev       = NaN;
        jclips(lb2_index).inlsa_mos = NaN;
        
        [jclips_filename, rem] = strtok(lb2_filename, '@');
        jclips(lb2_index).file_name     = {jclips_filename};  % CELL

        %Are we dealing with a parsed file?
        parsed_mat = strcat(working_dir,path_sep, 'bvqm-parsed-', jclips_filename, '.mat');
        if exist (parsed_mat, 'file')
            load (parsed_mat);
            jclips(lb2_index).loc_start = parsed(parsed_index).start_frame;
            jclips(lb2_index).loc_stop  = parsed(parsed_index).stop_frame;
            parsed_index = parsed_index + 1;
        else  % not a parsed file - but if bvqm-ptrunc file exists, set loc_stop to truncated value
% %             trunc_save_file=strcat(working_dir, path_sep, 'bvqm-ptrunk.mat');
% % 
% %             if exist(trunc_save_file, 'file')
% %                 % current clip=jclips_filename
% %                 load(trunc_save_file);
% %                  num_t=size(trunc_filename, 2);
% %                 tindex=1;
% %                 while tindex < num_t
% %                     trunc_fn=mat2str(cell2mat(trunc_filename(tindex)));
% %                     jclips_filename;
% %                     if strcmp(jclips_filename, trunc_fn)
% %                         match=1;
% %                         jclips(lb2_index).loc_start=1;
% %                         jclips(lb2_index).loc_stop=trunc_stop_frame(tindex);
% %                     else
% %                         jclips(lb2_index).loc_start = 1;   %?????
% %                         jclips(lb2_index).loc_stop  = num_frames;  %???????/
% %                     end
% %                     tindex=tindex+1;
% %                 end
% %             else
               
                jclips(lb2_index).loc_start = 1;
                jclips(lb2_index).loc_stop  = num_frames;

        end 
            % Set align_start/stop to NaN- 26May06-
            %         jclips(lb2_index).align_start = jclips(lb2_index).loc_start;
            %         jclips(lb2_index).align_stop  =
            %         jclips(lb2_index).loc_stop;
            jclips(lb2_index).align_start = NaN;
            jclips(lb2_index).align_stop = NaN;

            %jclips(lb2_index).image_size %%%%%%%%%
            [ds] = default_sroi(jclips(lb2_index).image_size);
            jclips(lb2_index).cvr = ds;

            jclips(lb2_index).viewers     = 0;
            jclips(lb2_index).hrc_definition   = {''};
            jclips(lb2_index).scene_definition = {''};
            jclips(lb2_index).subj_system = {''};
            jclips(lb2_index).inlsa_mos_2003 = NaN;
            jclips(lb2_index).scale.horizontal = 1000;
            jclips(lb2_index).scale.vertical   = 1000;


            if fid ~= -1; fclose(fid); end
        end
       
        % Check that YUV file sizes are legitimate:
        if isequal(fn_ext, '.yuv')
            exp_size=(num_frames*cols*rows*2);
            if exp_size ~= file_size
                if ~yuv_sz_err_dsp
                    % Display warning if there is an invalid YUV file
                    sl_line=sprintf('%8.f bytes ~= %dcols * %drows * %3.fframes', exp_size, cols, rows,num_frames);
                    qstnam='Unexpected YUV file size';
                    qststr={'YUV clips with an invalid file size exist:';' ';...
                        'Expected file size = 2 * #cols * #rows * #frames';...
%                         sl_line; ...
                        ' ';'Clip file sizes are not what was expected for the given image size.';...
                        'Please verify the image size.';
                        ' ';'Continuing with this indescrepancy will result in a fractional end frame.';...
                        'This fractional end frame will be rounded down, and the final results may be skewed.';...
                        ' ';'Press OK to continue with error.'
                        'Press CANCEL go back and to fix error.'};
                    do_this=questdlg(qststr, qstnam, 'OK', 'CANCEL', 'CANCEL');
                    waitfor(do_this);
                    yuv_sz_err_dsp=1;
                    if strcmp(do_this, 'CANCEL')
                        return
                    end
                end  % if strcmp(do_this..)
            end % if yuv_sz_err_dsp
        end % if YUV file
        % user was warned - fix loc_stop (whether it needs it or not):
        jclips(lb2_index).loc_stop=floor(jclips(lb2_index).loc_stop);
        lb2_index=lb2_index+1;
    end % END Create JClips
    
%
% TEMP: ANY clip longer than 15 secs: truncate to 15 secs 
%
% lb2_index=1;
% maxxx=0;
% while lb2_index <= lb2_size
%     fps=jclips(1).fps;
%     max_frame=fps*15;
%     if jclips(lb2_index).loc_stop > max_frame;
%         if maxxx==0 % first encounter -display msg:
%             msg={'Clips over 15 seconds in length were found.';' Only the first 15 seconds will be used.'};
%             h=helpdlg(msg, 'Clip(s) will be truncated');
%             waitfor(h);
%             maxxx=1;
%         end
%         jclips(lb2_index).loc_stop=max_frame;
%     end
%     lb2_index = lb2_index+1;
% end

% RESET test, scene, and HRC name to values entered by user (stored in vcd)
% loop through lb2 items and assign jclips.{test, scene, HRC} to values in  vcd
vcd_fn=fullfile(working_dir, 'bvqm-vcd.mat');
load (vcd_fn);

% 
% If all test names are the same, leave them - if different, return an
% error
if size(unique([jclips.test]),2) >1
    msg = sprintf('All video clips must have the same test name.\n\n');
    msg = sprintf('%sFile naming convention #1: <test>_<scene>_<HRC>\n', msg);
    msg = sprintf('%sFile naming convention #2: <scene>_<HRC>\n', msg);
    msg = sprintf('%sOtherwise scene and HRC are manually entered\n\n', msg);
    msg = sprintf('%sThe default test name is "bvqm".', msg);

    warndlg(msg, 'Multiple Tests Not Allowed');
    return
else
    tn=cell2mat(unique([jclips.test]));
end


% perform checks by comparing file contents with information in jclips
for cnt = 1:length(jclips),
    [~,ext]=strtok(jclips(cnt).file_name{1}, '.');
    
    if  strcmp(ext, '.avi'),
        avi_info=read_avi('info',fullfile(clip_dir,jclips(cnt).file_name{1}));
        
        % return an error if the frame rates differ significantly.
        fps = avi_info.FramesPerSecond;
        if fps < jclips(cnt).fps - 0.01 || fps > jclips(cnt).fps + 0.01,
            msg = sprintf('All video files must have the same frame rate.\n\n');
            msg = sprintf('%sFrame rate in header of AVI file "%s"\n is %f.', ...
                msg, jclips(cnt).file_name{1}, fps);
            msg = sprintf('%sExpected frame rate is %f', msg, jclips(cnt).fps);
            warndlg(msg, 'Incompatible Files');
            return
        end
    end

    % also return an error if images of 2+ sizes are selected.
    if jclips(1).image_size.cols ~= jclips(cnt).image_size.cols || ...
        jclips(1).image_size.rows ~= jclips(cnt).image_size.rows,
            msg = sprintf('All video files must have the same number of rows and columns\n');
            msg = sprintf('%sImage size in header of AVI file "%s"\n', ...
                msg, jclips(cnt).file_name{1});
            msg = sprintf('%sis (%d,%d) Expected value is (%d,%d)', msg,...
                jclips(cnt).image_size.rows, ...
                jclips(cnt).image_size.cols, ...
                jclips(1).image_size.rows, ...
                jclips(1).image_size.cols);
            warndlg(msg, 'Incompatible Files');
            return
    end
end




% Edit jclips based on data in vcd
lb2_idx=1;
while lb2_idx <= lb2_size
%    jclips(lb2_idx).test = {'bvqm'}; %  vcd(lb2_idx).test;
%   jclips(lb2_idx).test = {'ps2on'};
    jclips(lb2_idx).test = {tn}; %  vcd(lb2_idx).test;
    jclips(lb2_idx).scene={vcd(lb2_idx).scene};
    jclips(lb2_idx).hrc=  {vcd(lb2_idx).hrc};
    jclips(lb2_idx).file_name={vcd(lb2_idx).filename};
    lb2_idx=lb2_idx+1;
    
%     save_file= fullfile(working_dir, 'jclips');
%     save (save_file, 'jclips', 'working_dir');
end
    

%
% Create GTests from data in GClips
%
num_clips = size(jclips,2);
test_names = unique([jclips.test]);
num_unique_test_names = size(test_names,2);
lb2_filelist=get(handles.listbox2, 'String');
path = get(handles.clip_directory_text3, 'String');

if ispc ; path_sep = '\'; else; path_sep='/'; end
path = strcat(path, path_sep);

% This for loop obsoleted b/c all test names set to be the same
%       leave intact for future releases, though.
for test_index = 1:size(test_names,2)
    viewers(test_index) = 0 ; mos_worst(test_index) = [999] ; mos_best(test_index) = [-999];
    test_name = test_names(test_index);
    jtests(test_index).name = test_name;
    jtests(test_index).path = {path};
    jtests_scene_index = 1;
        
    for clip_index = 1:num_clips
        if strcmp(jclips(clip_index).test , test_name)
            % compute best and worst mos, num viewers for each test:
%              jtests(test_index).scenes(jtests_scene_index) = jclips(clip_index).scene;
%              jtests(test_index).hrc(jtests_scene_index) = jclips(clip_index).hrc;

              jtests(test_index).scenes = jclips(clip_index).scene;
              jtests(test_index).hrc = jclips(clip_index).hrc;


            jclips_viewers = jclips(clip_index).viewers;
            if isempty(jclips_viewers); jclips_viewers = 0 ; end
            
            viewers(test_index) = viewers(test_index) + jclips_viewers;
            jtests(test_index).viewers = viewers(test_index);
            % BUT- if it's a parsed clip, just count the # of viewers once:
            % to be implemented...
            
            mos_current = jclips(clip_index).mos;

            if mos_current < mos_worst(test_index)
                mos_worst(test_index) = mos_current;
                jtests(test_index).mos_worst = mos_worst(test_index);
            end

            if mos_current > mos_best(test_index)
                mos_best(test_index) = mos_current;
                jtests(test_index).mos_best = mos_best(test_index);
            end
            
            jtests_scene_index = jtests_scene_index + 1;
        end
       
        % fix jclips format errors:
%         if iscell(jclips(clip_index).scene)
%             jclips(clip_index).scene = cell2mat(jclips(clip_index).scene);
%         end
% 
%         if iscell(jclips(clip_index).hrc)
%             jclips(clip_index).hrc = cell2mat(jclips(clip_index).hrc);
%         end

        
        clip_index = clip_index + 1;
    end
%     jtests(test_index).scenes = unique(jtests(test_index).scenes);
%     jtests(test_index).hrc = unique(jtests(test_index).hrc);
    
    jtests.scenes = unique([jclips.scene]);
    jtests.hrc    = unique([jclips.hrc]);
    test_index = test_index + 1;
end  % END Create JTests



[clip_err, clip_msg]=check_clips(jclips, jtests, 'quiet');


if isequal(clip_err,1)
    error_str={'The clip structure created has errors.'; 'Please contact engineers at ITS to resolve this problem.'}
    h=errordlg(error_str, 'CLIP STRUCTURE ERROR');
    waitfor(h);
end

if isequal(clip_err, 2)
    error_str2='Error(s) must be fixed before continuing.';
    error_str={clip_msg; ''; error_str2};
    h=errordlg(error_str, 'Error');
    waitfor(h);
end

if isequal(clip_err, 3)
    b1='Ignore'; b2='Fix';
 
    button = questdlg(clip_msg,'Warning',b1, b2, b1);
    waitfor(button);
    
    if strcmp(button, b1)
        clip_err=0;
    end
    
    if strcmp(button,b2)
        clip_err=3;
    end
end

% sort jclips by scene then HRC.
 [order] = sort_clips_by ('scene', jclips, jtests);
 order2 = [];
 for loop=1:length(order),
     order2 = [order2, order{loop}];
 end
 jclips = jclips(order2);
 
if isequal(clip_err,0) % NO ERRORS - continue to cal window
    % parse if needed
    [jclips,jtests,status] = bvqm_pc_parse(jclips,jtests);
    
    % continue to cal window
    clip_dir = get(handles.clip_directory_text3, 'String');
    if ispc ; path_sep = '\'; else; path_sep='/'; end
    save_file= strcat(working_dir, path_sep, 'jclips');
    save (save_file, 'jclips', 'jtests', 'working_dir','status');

    %try delete('bvqm-test_data.mat'); catch ; err=1 ; end

    % IF 'exit' button was pressed, just delete window, don't run cal,
    if exit_pb == 1
        delete(handles.bvqm_pc);  % was ... ,1);   ??
        fclose('all');
    else
%        delete(handles.bvqm_pc);
        pause(0.25);
        cr=0; % 0 signifies no crash recovery
        pass={working_dir, cr};
%         bvqm_cal(working_dir,cr); 
%        bvqm_pc_cal(pass); 

% fix for compiler:
%        bvqm_pc_cal('bvqm_pc_cal_OpeningFcn',pass);
        bvqm_pc_cal([], pass);
        pause(0.25);
        delete(handles.bvqm_pc);
    end
end


% --- Executes on button press in pushbutton14_exit.
function pushbutton14_exit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

working_dir=get(handles.text40, 'String');
clip_dir=get(handles.clip_directory_text3, 'String');
if ispc ; path_sep = '\'; else; path_sep='/'; end

warning off all; % no usefull info contained therein

qstring='Would you like to save before exiting?';
b1='Save & Exit'; b2='Exit'; b3='Cancel';

button = questdlg(qstring, 'Confirm Save/Exit', b1, b2, b3, b2);

if (strcmp(button, b1))  % Save structures, then exit
    %   {'Message line 1';'Message line 2'}
    % Save jclips.mat, only if there is something to save...

    lb2_filelist = get(handles.listbox2,'String');
    if size(lb2_filelist) ~=0
        title='Clips structure saved';
        msg={'Results are saved in the working directory:'; working_dir;...
            'File name: jclips.mat'; };% ''; 'Press OK to exit.'};
        h=msgbox(msg, title, 'help');
        waitfor(h);
        continue_pushbutton13_Callback(hObject, eventdata, handles,1);
    else
       % close('bvqm_pc');
        delete(handles.bvqm_pc);
        fclose('all');
    end

    try  rmdir (strcat(working_dir, path_sep,'feature_*'), 's') ;   catch; err=1; end
    %    try delete (strcat(working_dir, path_sep,'bvqm-status.mat')) ;  catch; err=1; end
    try delete (strcat(working_dir, path_sep,'bvqm-opt_*.mat'));    catch; err=1; end
    try delete (strcat(working_dir, path_sep, 'bvqm-parsed-*.mat')); catch; err=1; end
    try delete (strcat(working_dir, path_sep, 'bvqm-ptrunk.mat')); catch; err=1; end

end

if (strcmp(button, b2)) % delete and just exit
    try  rmdir (strcat(working_dir, path_sep,'feature_*'), 's') ;   catch; err=1; end
    try delete (strcat(working_dir, path_sep,'bvqm-status.mat')) ;  catch; err=1; end 
    try delete (strcat(working_dir, path_sep,'bvqm-opt_*.mat'));    catch; err=1; end
    try delete (strcat(working_dir, path_sep, 'bvqm-parsed-*.mat')); catch; err=1; end
%    try delete (strcat(working_dir, path_sep, 'bvqm-test_data.mat')); catch ; err=1 ; end
    try delete (strcat(working_dir, path_sep,'jclips.mat'));        catch; err=1; end
%    try delete (strcat(working_dir, path_sep,'bvqm_*report*.txt')); catch; err=1; end
    try delete (strcat(clip_dir, path_sep, '%bvqm_*_*.*'));      catch; err=1; end
    try delete (strcat(working_dir, path_sep, 'bvqm-ptrunk.mat')); catch; err=1; end
    try delete (fullfile(working_dir, 'bvqm-vcd.mat')); catch; ;err=1; end
    
     try load (fullfile(working_dir, 'bvqm-linked.mat'))
        num_linked=size(linked_file,2);
        count=1;
        while count <= num_linked
            linked_file{count};
            try delete (linked_file{count}) 
            count=count+1;
            catch
                err=1;
            end
        end
        delete(fullfile(working_dir, 'bvqm-linked.mat'));
    catch
        nosuchfile=1;
    end

    delete(handles.bvqm_pc);
    pause(1);
end


% Check to see that image size (rows & cols), and fps are entered.
%   if not, offer suggestions
function check_imagesize(hObject, eventdata, handles)

clip_dir = get(handles.clip_directory_text3,'String');
working_dir = get(handles.text40, 'String');
lb2_filelist = get(handles.listbox2,'String');
lb2_index = get(handles.listbox2, 'Value');
lb2_filename=cell2mat(lb2_filelist(lb2_index));
[pre,ext]=strtok(lb2_filename, '.');

% if first file chosen is an avi file, get size data from aviinfo fcn:
if  (strcmp(ext, '.avi')) && (isequal(lb2_index,1))
    avi_info=aviinfo(fullfile(clip_dir,lb2_filename));
    cols=avi_info.Width;
    rows=avi_info.Height;
    fps =avi_info.FramesPerSecond;
    set(handles.image_size_cols_edit, 'String', cols);
    set(handles.image_size_rows_edit, 'String', rows);
    set(handles.fps_edit, 'String', fps);
    % assign proper video standard based upon image size:
    if isequal(rows,486) && isequal(cols,720)
        set(handles.video_std_popupmenu1, 'Value', 1.0); %1=lff interlaced lower field first
    elseif isequal(rows,480) && isequal(cols,720)
        set(handles.video_std_popupmenu1, 'Value', 1.0); %1=lff interlaced lower field first
    elseif isequal(rows,576) && isequal(cols,720)
        set(handles.video_std_popupmenu1, 'Value', 2.0); % 2=uff interlaced upper field first
    elseif isequal(rows,720) && isequal(cols,1280)
        set(handles.video_std_popupmenu1, 'Value', 3); %3=progressive
    elseif isequal(rows,1080) && isequal(cols,1920)
        set(handles.video_std_popupmenu1, 'Value', 2.0); %2=uff interlaced upper field first
    else
        set(handles.video_std_popupmenu1, 'Value', 3); %3=progressive
    end
end

rows = str2num(get(handles.image_size_rows_edit, 'String'));
cols = str2num(get(handles.image_size_cols_edit, 'String'));
fps = str2num(get(handles.fps_edit, 'String'));

% % Check to see if rows,cols,fps are valid data:
% is_str = (strcmp('',rows) | strcmp('',cols) | strcmp('',fps));
% is_emt = (isempty(rows) | isempty(cols) | isempty(fps));
% is_num = (isnumeric(rows) | isnumeric(cols) | isnumeric(fps));

%try
    
while (isempty(rows) | isempty(cols) | isempty(fps))
    button = bvqm_pc_sz;
    waitfor(button);

    if strcmp(button, 'NONE')  % window was killed
        return
    end
    
    if strcmp(button, 'NTSC')
        rows=486; cols=720; fps=29.97;
        set(handles.image_size_rows_edit, 'String', rows);
        set(handles.image_size_cols_edit, 'String', cols);
        set(handles.fps_edit, 'String', fps);
        set(handles.video_std_popupmenu1, 'Value', 1.0); % 1=lff
    end

    if strcmp(button, 'PAL')
        rows=576; cols=720; fps=25.00;
        set(handles.image_size_rows_edit, 'String', rows);
        set(handles.image_size_cols_edit, 'String', cols);
        set(handles.fps_edit, 'String', fps);
        set(handles.video_std_popupmenu1, 'Value', 2.0); % 2=uff
    end

    if strcmp(button, 'HDTV720p')
        rows=720; cols=1280; fps=29.97;
        set(handles.image_size_rows_edit, 'String', rows);
        set(handles.image_size_cols_edit, 'String', cols);
        set(handles.fps_edit, 'String', fps);
        set(handles.video_std_popupmenu1, 'Value', 2.0); % 2=uff
    end

    if strcmp(button, 'HDTV1080i')
        rows=1080; cols=1920; fps=29.97;
        set(handles.image_size_rows_edit, 'String', rows);
        set(handles.image_size_cols_edit, 'String', cols);
        set(handles.fps_edit, 'String', fps);
        set(handles.video_std_popupmenu1, 'Value', 2.0); % 2=uff
    end

    if strcmp(button, 'manual')

        prompt={'Image Size- COLS:', 'Image Size - ROWS:', 'Frames per Second:'};
        title= 'Enter Values';
        [answer] = inputdlg(prompt,title);

        num_els=size(answer,1);
        if num_els >2
            rows=str2num(cell2mat(answer(2)));
            cols=str2num(cell2mat(answer(1)));
            fps=str2num(cell2mat(answer(3)));
            answer=[rows cols fps];

            % for now, don't display error msg for bad input - just force
            % reentry...
            %             if ~isnumeric(answer)
            %                 h=errordlg('All fields must be numeric', 'Invalid data');
            %                 waitfor(h);
            %                 check_imagesize(hObject, eventdata, handles);
            %             else
            %                 set(handles.image_size_rows_edit, 'String', rows);
            %                 set(handles.image_size_cols_edit, 'String', cols);
            %                 set(handles.fps_edit, 'String', fps);
            %             end
            set(handles.image_size_rows_edit, 'String', rows);
            set(handles.image_size_cols_edit, 'String', cols);
            set(handles.fps_edit, 'String', fps);

        end
    end

    % update:
    rows = str2num(get(handles.image_size_rows_edit, 'String'));
    cols = str2num(get(handles.image_size_cols_edit, 'String'));
    fps = str2num(get(handles.fps_edit, 'String'));
    guidata(hObject, handles);
end


% % Get rows, cols, fps data from user
% function get_rcf
%     prompt={'Image Size- COLS:', 'Image Size - ROWS:', 'Frames per Second:'};
%     title= 'Enter Values';
%     [answer] = inputdlg(prompt,title);
%     if ~isnumeric(answer)
%          errordlg('All fields must be numeric', 'Invalid data');
%     end


%NOTE: NOT currently used - 02Feb06
    % This function will save the rows, cols, and fps fields for each clip -
% for later display on the screen and for the jclips structure:
function save_tdf(hObject, eventdata, handles)
% No matter what type of clip, if user enters image size data (rows or
% cols) or fps, store that data for the clip:
if ispc ; path_sep = '\'; else; path_sep='/'; end
working_dir = get(handles.text40, 'String');
test_data_file = strcat(working_dir, path_sep, 'bvqm-test_data.mat');

tdf_rows=get(handles.image_size_rows_edit, 'String');
tdf_cols=get(handles.image_size_cols_edit, 'String');
tdf_fps= get(handles.fps_edit, 'String');
save(test_data_file, 'tdf_rows', 'tdf_cols', 'tdf_fps');


% --- Executes on button press in pushbutton15_clipdir.
function pushbutton15_clipdir_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15_clipdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

wp_check_clp=1;
clip_dir_save = get(handles.clip_directory_text3,'String');

while isequal(wp_check_clp, 1)
            [clip_dir] = uigetdir('','BVQM: Select a CLIP directory:');
            if isequal(clip_dir,0);
                clip_dir=clip_dir_save;
            end
           [stat, msg, msg_id] = fileattrib(clip_dir);
            wp_check_clp=0;
            set(handles.clip_directory_text3,'String', clip_dir);
end

if [clip_dir] == 0
    clip_dir = clip_dir_save;
    return
end

set(handles.clip_directory_text3,'String', clip_dir);
% Populate the listbox
load_listbox(clip_dir,handles);
% Remove any lines in lb2:
set(handles.listbox2, 'String', '');

% clear test data:
set(handles.image_size_rows_edit, 'String', '');
set(handles.image_size_cols_edit, 'String', '');
set(handles.fps_edit, 'String', '');
set(handles.video_std_popupmenu1, 'Value', 3.0); % 3=progressive.
set(handles.length_text33, 'String', '');

% clear parsing panel:
% set(handles.auto_parse_uipanel3, 'Visible', 'off');

% clear video clip data box:
set(handles.test_edit, 'String', '');
set(handles.scene_edit1, 'String', '');
set(handles.hrc_edit2, 'String', '');
set(handles.original_clip_checkbox1, 'Value', 0);
set(handles.text6, 'String', '');  % clip file name



% --- Executes on button press in pushbutton16_test_data.
function pushbutton16_test_data_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16_test_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% save_rows=get(handles.image_size_rows_edit, 'String');
% save_cols=get(handles.image_size_cols_edit, 'String');
% save_fps= get(handles.fps_edit, 'String');
%save_ss=  get(handles.

% Gather data needed to fill in optional data
lb2_index = get(handles.listbox2,'Value');
file_list = get(handles.listbox2,'String');

% if no files selected first, do nothing:
if isempty(file_list)
    return
else
    % clear test data:
    set(handles.image_size_rows_edit, 'String', '');
    set(handles.image_size_cols_edit, 'String', '');
    set(handles.fps_edit, 'String', '');
    set(handles.video_std_popupmenu1, 'Value', 3.0); % 3=progressive.
    set(handles.length_text33, 'String', '');

    check_imagesize(hObject, eventdata, handles);
end




% --- Executes on key press with focus on add_pushbutton2 and none of its controls.
function add_pushbutton2_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to add_pushbutton2 (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
