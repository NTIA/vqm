function varargout = bvqm_pc_report(varargin)
% BVQM_PC_REPORT M-file for bvqm_pc_report.fig
%      BVQM_PC_REPORT, by itself, creates a new BVQM_PC_REPORT or raises the existing
%      singleton*.
%
%      H = BVQM_PC_REPORT returns the handle to a new BVQM_PC_REPORT or the handle to
%      the existing singleton*.
%
%      BVQM_PC_REPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BVQM_PC_REPORT.M with the given input arguments.
%
%      BVQM_PC_REPORT('Property','Value',...) creates a new BVQM_PC_REPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bvqm_pc_report_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bvqm_pc_report_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bvqm_pc_report

%
% WRITTEN BY : Mark A. McFarland, P.E.    303/497-4132
%


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bvqm_pc_report_OpeningFcn, ...
                   'gui_OutputFcn',  @bvqm_pc_report_OutputFcn, ...
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

function bvqm_pc_report_passer(working_dir, s_rpt, d_rpt, clip_dir)
% use command 'findobj' to find object handles....

handles=guidata('bvqm_pc_report');

handles.working_dir=(varargin{2});
set(handles.text9, 'String', working_dir);

handles.clip_dir=(varargin{5});
set(handles.text11, 'String', clip_dir);

handles.s_rpt=(varargin{3}); % summary rpt file name
handles.d_rpt=(varargin{4}); % detail rpt file name
set(handles.text10, 'String', s_rpt);
set(handles.text11, 'String', d_rpt);



%refresh_report(hObject, eventdata, handles, model_to_run);
%refresh_report(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function bvqm_pc_report_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bvqm_pc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Add the current directory to the path, as the pwd might change thru' the
% gui. Remove the directory from the path when gui is closed
% (See bvqm_pc_DeleteFcn)
% setappdata(hObject, 'StartPath', pwd);
% addpath(pwd);
% set(bvqm_pc, 'CloseRequestFcn', 'bvqm_pc_CloseRequestFcn')

% % Display logo:
% x = imread('vqm_logo.jpg');
% axes(handles.axes1); % makes axis1 current
% %set(handles.axes1, 'CurrentAxes', axes1)
% imagesc(x);
% set(gca,'Visible','off')

% --- Executes during object deletion, before destroying properties.
function bvqm_pc_report_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to bvqm_pc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Remove the directory added to the path in the bvqm_pc_CreateFcn.
% if isappdata(hObject, 'StartPath')
%     rmpath(getappdata(hObject, 'StartPath'));
% end
%pb4_quit_Callback(hObject, eventdata, handles)



% --- Executes just before bvqm_pc_report is made visible.
function bvqm_pc_report_OpeningFcn(hObject, eventdata, handles,varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bvqm_pc_report (see VARARGIN)

% Store working_dir for future use:

% working_dir=(varargin{1});
% set(handles.text9, 'String', working_dir);
% 
% clip_dir=(varargin{4});
% set(handles.text11, 'String', clip_dir);
% 
% s_rpt=(varargin{2}); % summary rpt file name
% d_rpt=(varargin{3}); % detail rpt file name
% set(handles.text10, 'String', s_rpt);
% set(handles.text11, 'String', d_rpt);

%refresh_report(hObject, eventdata, handles, model_to_run);

% handles=guidata(bvqm_pc_report);

% handles.working_dir=(varargin{2});
% handles.clip_dir=(varargin{5});
% handles.s_rpt=(varargin{3}); % summary rpt file name
% handles.d_rpt=(varargin{4}); % detail rpt file name
% % guidata(bvqm_pc_report, handles)

% handles.working_dir=(varargin{1});
% handles.clip_dir=(varargin{4});
% handles.s_rpt=(varargin{2}); % summary rpt file name
% handles.d_rpt=(varargin{3}); % detail rpt file name
% % guidata(bvqm_pc_report, handles)

handles.working_dir=(varargin{2});
handles.clip_dir=(varargin{8});
handles.s_rpt=(varargin{4}); % summary rpt file name
handles.d_rpt=(varargin{6}); % detail rpt file name

set(handles.text9, 'String', handles.working_dir);
set(handles.text11, 'String', handles.clip_dir);
set(handles.text10, 'String', handles.s_rpt);
set(handles.text11, 'String', handles.d_rpt);
%%%%

refresh_report(hObject, eventdata, handles);

% Choose default command line output for bvqm_pc_report
handles.output = hObject;

%Set the close function callback here.
set(handles.output,'CloseRequestFcn',{@my_close_fcn, handles})
%set(handles.output,'CloseRequestFcn',{@pb4_quit_Callback, handles.output})

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes bvqm_pc_report wait for user response (see UIRESUME)
% uiwait(handles.bvqm_pc_report);

%set up export to spreadsheet buttons for linux or windows:
if ~ispc
    set(handles.pushbutton7, 'String', 'View in OpenOffice');
    set(handles.pushbutton7, 'FontSize', 9);
%    set(handles.pushbutton7, 'Visible', 'off');
    set(handles.pushbutton12, 'String', 'Create CSV Table');
    set(handles.pushbutton12, 'TooltipString', 'Create a comma separated value (CSV) table');
end
    
set(handles.axes1,'visible','off');


function my_close_fcn(hObject, eventdata, handles)

if ispc ; path_sep = '\'; else; path_sep='/'; end
working_dir=get(handles.text9, 'String');
clip_dir=get(handles.text11, 'String');


b1='Save & Exit'; b2='Exit'; b3='Cancel';

% button = questdlg(qstring, 'Confirm Save/Exit', b1, b2, b3, b2);

really= questdlg('Save intermediate results before quit?', 'Confirm Save/Exit', b1,b2,b3,b2);
waitfor(really);
 
warning off all;
 if strcmp(really, b2)  
    try delete (fullfile(working_dir, 'bvqm-vcd.mat')); catch ; err=1 ; end  % video clip data filename
    try  rmdir (strcat(working_dir, path_sep,'feature_*'), 's') ;   catch; err=1; end
    try delete (strcat(working_dir, path_sep,'jclips.mat'));        catch; err=1; end
    try delete (strcat(working_dir, path_sep,'bvqm-status.mat')) ;  catch; err=1; end
    try delete (strcat(working_dir, path_sep,'bvqm-opt_*.mat'));    catch; err=1; end
    try delete (strcat(working_dir, path_sep,'vqm_*.mat'));         catch; err=1; end
    try  rmdir (strcat(working_dir, path_sep,'vfd_results_psnr'), 's');  catch; err=1; end
    try  rmdir (strcat(working_dir, path_sep,'vfd_results_vqm'), 's');   catch; err=1; end
    try delete (strcat(working_dir, path_sep,'vfd_results.csv'));  catch; err=1; end
    try delete (strcat(working_dir, path_sep,'vfd.mat'));           catch; err=1; end
    
    delete(handles.bvqm_pc_report); % close all; %('bvqm_pc_report');
    fclose('all');
    pause(1);
 end
 
 if strcmp(really, b1)  % don't delete anything -- just exit.
    msg={'Results are saved in the working directory:'; working_dir};
    h=msgbox(msg, 'Results Saved', 'help');
    waitfor(h);
 
    delete(handles.bvqm_pc_report); % close all; %('bvqm_pc_report');
    fclose('all');
    pause(1);
 end

%function refresh_report(hobject, eventdata, handles, model_to_run)
function refresh_report(hobject, eventdata, handles)

s_rpt=get(handles.text10, 'String'); % summary rpt file name
d_rpt=get(handles.text11, 'String'); % detail rpt file name

working_dir=get(handles.text9, 'String');
lfile=fullfile(working_dir, 'bvqm-status');

load (lfile);
% bvqm-status;
summary_bool=get(handles.radiobutton1,'Value');
detailed_bool=get(handles.radiobutton1,'Value');

if summary_bool==1
    filename=s_rpt;
else
    filename=d_rpt;
end

fid = fopen(filename,'r');

%Read in Report from Text File.
file = textscan(fid, '%s', 'delimiter', '\n');  %1by1 cell array of cell arrays
report=char(file{1,1});  % report is a character array
fclose(fid);
if size(report,2) > 83,
    report = report(:,1:83);
end

%Make this ASCII text look pretty by indenting everything and
%truncating it to make it fit in the box
% width=63;
% leftmargin(1:length(report(:,1)),1)=' ';
% head_foot(1:width+1)=' ';
% report_formatted=  [head_foot; leftmargin report(:,1:width); head_foot ];
% set(handles.text1, 'String', report_formatted); 
set(handles.text1, 'String', report);

if ispc ; path_sep = '\'; else; path_sep='/'; end
% set(handles.txt_location, 'String', strcat(strcat(pwd,path_sep),filename));
set(handles.txt_location, 'String', filename);


% --- Outputs from this function are returned to the command line.
function varargout = bvqm_pc_report_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pb3_back.
function pb3_back_Callback(hObject, eventdata, handles)
% hObject    handle to pb3_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pb4_quit.
function pb4_quit_Callback(hObject, eventdata, handles)
% hObject    handle to pb4_quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ispc ; path_sep = '\'; else; path_sep='/'; end
working_dir=get(handles.text9, 'String');
clip_dir=get(handles.text11, 'String');


b1='Save & Exit'; b2='Exit'; b3='Cancel';

% button = questdlg(qstring, 'Confirm Save/Exit', b1, b2, b3, b2);

really= questdlg('Save intermediate results before quit?', 'Confirm Save/Exit', b1,b2,b3,b2);
waitfor(really);
 
warning off all;
 if strcmp(really, b2)  
    try delete (fullfile(working_dir, 'bvqm-vcd.mat')); catch ; err=1 ; end  % video clip data filename
    try  rmdir (strcat(working_dir, path_sep,'feature_*'), 's') ;   catch; err=1; end
    try delete (strcat(working_dir, path_sep,'jclips.mat'));        catch; err=1; end
    
    try delete (strcat(working_dir, path_sep,'bvqm-status.mat')) ;  catch; err=1; end
    try delete (strcat(working_dir, path_sep,'bvqm-opt_*.mat'));    catch; err=1; end
    try delete (strcat(working_dir, path_sep,'vqm_*.mat'));         catch; err=1; end
    try  rmdir (strcat(working_dir, path_sep,'vfd_results_psnr'), 's');  catch; err=1; end
    try  rmdir (strcat(working_dir, path_sep,'vfd_results_vqm'), 's');   catch; err=1; end
    try delete (strcat(working_dir, path_sep,'vfd_results.csv'));  catch; err=1; end
    try delete (strcat(working_dir, path_sep,'vfd.mat'));           catch; err=1; end
    
    delete(handles.bvqm_pc_report); % close all; %('bvqm_pc_report');
    fclose('all');
    pause(1);
 
 end
 
 if strcmp(really, b1)  % don't delete anything -- just exit.
     
    msg={'Results are saved in the working directory:'; working_dir};
    h=msgbox(msg, 'Results Saved', 'help');
    waitfor(h);

    delete(handles.bvqm_pc_report); % close all; %('bvqm_pc_report');
    fclose('all');
    pause(1);
 end
     

% --- Executes on selection change in popupmenu1. - SELECT GRAPH
%function popupmenu1_Callback(hObject, eventdata, handles,model_to_run)
function popupmenu1_Callback(hObject, eventdata, handles, model_op)

% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

working_dir=get(handles.text9, 'String');
load (fullfile(working_dir, 'bvqm-status'));   % to extract model_op structure

% This gets called whenever a user chooses a graph.
menu_index=get(handles.popupmenu1,'Value');
set(handles.listbox2,'Visible','off');
set(handles.text5,'Visible','off');
set(handles.pushbutton8,'Enable','on');
axes(handles.axes1);
if (menu_index==1)
    cla;
    title('');
    xlabel('');
    ylabel('');
    set(handles.pushbutton8,'Enable','off');
    set(handles.axes1,'visible','off');

elseif menu_index==2
    if length(model_op.data(1,:)) == 1,
        two_wide = [model_op.data(1,:), nan];
        barh(two_wide,'facecolor',[0.3 0.3 1]);
    else
        barh(model_op.data(1,:),'facecolor',[0.3 0.3 1]);
    end
    if length(model_op.clip_name) < 6,
        tic_labels = model_op.clip_name(1:length(model_op.clip_name));
        tic_nums = 1:length(model_op.clip_name);
        ylabel('Clip','FontSize',10);
    else
        width = ceil(length(model_op.clip_name)/6);
        tic_nums = 1:width:length(model_op.clip_name);
        tic_labels = model_op.clip_name(tic_nums);
        ylabel('Clips (alphabetically)','FontSize',10);
    end
    set(handles.axes1,'YTick',tic_nums);
    set(handles.axes1,'YTickLabel',tic_labels);
    set(handles.axes1,'FontSize',8);
    xlabel('VQM Score');
    title(strcat(model_op.par_name{1},': Quality by Clip'), 'Interpreter','none','FontSize',10);
elseif menu_index==3
    %plot one bar of average VQM per HRC
    hrc_stats=ave_par_values(model_op,'scene');
    if length(hrc_stats.data(1,:)) == 1,
        two_wide = [hrc_stats.data(1,:), nan];
        barh(two_wide,'facecolor',[0.3 0.3 1]);
    else
        barh(hrc_stats.data(1,:),'facecolor',[0.3 0.3 1]);
    end
    [a,b]=strtok(hrc_stats.clip_name, '_');
    [c,d]=strtok(b,'_'); [hrc_name,e]=strtok(d, '_');
     hrc_stats.hrc_name=hrc_name;
    if length(hrc_stats.hrc_name) < 6,
        tic_labels = hrc_stats.hrc_name(1:length(hrc_stats.hrc_name));
        tic_nums = 1:length(hrc_stats.hrc_name);
        ylabel('HRC','FontSize',10);
    else
        width = ceil(length(hrc_stats.hrc_name)/6);
        tic_nums = 1:width:length(hrc_stats.hrc_name);
        tic_labels = hrc_stats.hrc_name(tic_nums);
        ylabel('HRC (alphabetically)','FontSize',10);
    end
    set(handles.axes1,'YTick',tic_nums);
    set(handles.axes1,'YTickLabel',tic_labels);
    set(handles.axes1,'FontSize',8);
    xlabel('Average VQM Score');
    title(strcat(model_op.par_name{1},': Quality by HRC'), 'Interpreter','none','FontSize',10);    
elseif menu_index==4
    %plot one bar of average VQM per scene across HRC's
    scene_stats=ave_par_values(model_op,'hrc');
    if length(scene_stats.data(1,:)) == 1,
        two_wide = [scene_stats.data(1,:), nan];
        barh(two_wide,'facecolor',[0.3 0.3 1]);
    else
        barh(scene_stats.data(1,:),'facecolor',[0.3 0.3 1]);
    end
    [a,b]=strtok(scene_stats.clip_name, '_');
    [scn_name,d]=strtok(b,'_'); 
    scene_stats.scene_name=scn_name;
    if length(scene_stats.scene_name) < 6,
        tic_labels = scene_stats.scene_name(1:length(scene_stats.scene_name));
        tic_nums = 1:length(scene_stats.scene_name);
        ylabel('Scene','FontSize',10);
    else
        width = ceil(length(scene_stats.scene_name)/6);
        tic_nums = 1:width:length(scene_stats.scene_name);
        tic_labels = scene_stats.scene_name(tic_nums);
        ylabel('Scene (alphabetically)','FontSize',10);
    end
    set(handles.axes1,'YTick',tic_nums);
    set(handles.axes1,'YTickLabel',tic_labels);
    set(handles.axes1,'FontSize',8);
    xlabel('Average VQM Score');
    title(strcat(model_op.par_name{1},': Quality by Scene'), 'Interpreter','none','FontSize',10);    
elseif menu_index==5
    %plot all the paramaters for a single scene
    %enable the listbox to select the scene
    set(handles.listbox2,'Visible','on');
    set(handles.text5,'Visible','on');
    set(handles.listbox2, 'String',  char(model_op.clip_name));
    refresh_scene_detail(hObject, eventdata, handles,model_op);
else
    cla;
    title('');
    xlabel('');
    ylabel('');
    set(handles.pushbutton8,'Enable','off');
    set(handles.axes1,'visible','off');
end

%
function refresh_scene_detail(hObject, eventdata, handles, model_op)

working_dir=get(handles.text9, 'String');
load (fullfile(working_dir, 'bvqm-status'));   % to extract model_op structure

%This function looks at the selected scene and plots its details
selected=get(handles.listbox2,'Value');
bar(model_op.data(:,selected),'facecolor',[0.3 0.3 1]);
parlabels=model_op.par_name;
parlabels{1}='VQM';  %Use "VQM" instead of model name (Model Name Appears in Title)
set(handles.axes1,'XTickLabel',parlabels);
set(handles.axes1,'FontSize',8);
xlabel('Parameter');
ylabel('Value');
title(strcat(model_op.par_name{1}, strcat(' Scene Details: ',model_op.clip_name{selected})), 'Interpreter','none','FontSize',10);    


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radiobutton1. - SUMMARY RPT
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

refresh_report(hObject, eventdata, handles);

% --- Executes on button press in radiobutton2. - DETAIL REPORT.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

refresh_report(hObject, eventdata, handles);


%View in Excel
% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)

% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%edisplay jclips and model_op into their own spreadsheets 

working_dir=get(handles.text9, 'String');
load (fullfile(working_dir, 'bvqm-status'));
if ispc ; path_sep = '\'; else; path_sep='/'; end

clk=fix(clock);
hr=sprintf('%02.0f',clk(4)); mn=sprintf('%02.0f',clk(5)); sc=sprintf('%02.0f',clk(6));
time=strcat(hr,'_',mn);

jclips_file_out=strcat(working_dir,path_sep,char(jtests.name),'_clips_',date,'@',time, '.csv');
model_file_out=strcat(working_dir,path_sep,char(jtests.name),'_model_',date,'@', time, '.csv');

[jclipsdata jclipshead modeldata modelhead unitedheader]= prepexport(jclips, jtests, model_op, cal_type,hObject, eventdata, handles);

% write out CSV files.
make_csv(jclips_file_out, jclipsdata,unitedheader, jclipshead);
make_csv(model_file_out, modeldata,unitedheader, modelhead);
if ispc
    uiwait(msgbox(sprintf('MS-Excel not available.  Use another spreadsheet application open comma-delimited files "%s" and "%s"', jclips_file_out, model_file_out),...
        'Open Spreadsheet', 'none'));
else
    sys_cmd=['./oo_view.sh ' jclips_file_out ' ' model_file_out];

    [status,result] = system(sys_cmd);

    % Catch errors from system command:
    if ~(isequal(status, 0))
        h=errordlg({'Error opening OpenOffice Spreadsheet.'; ' ';'Is OpenOffice installed on this computer?'},'Export Error');
        waitfor(h);
    end           
end

 
function [jclipsdata  jclipshead modeldata modelhead unitedheader]= prepexport(jclips, jtests, model_op, cal_type,hObject, eventdata, handles);
%Function prepexport prepares the data about jclips and model_op for
%export to excel or csv

working_dir=get(handles.text9, 'String');
load (fullfile(working_dir, 'bvqm-status'));

jclipshead= {'test','scene','hrc','filename','location start','location stop',...
    'align start','align stop','video standard','fps','image size: rows',...
    'image size: columns', 'spatial: horizontal', 'spatial: vertical',...
    'scale: horizontal', 'scale: vertical', 'luminance gain',...
    'luminance offset', 'cvr: top', 'cvr: bottom', 'cvr: left', 'cvr: right',...
    'mos' ...
    };

%unwrap some of the multi-level structs
img_size=[jclips.image_size];
spatial=[jclips.spatial];
cvr=[jclips.cvr];
scale=[jclips.scale];

% Replace NaNs with [] for jclips.mos,inlsa_mos,stdev for display in excel:
for count=1:size(jclips,2)
    if isnan(jclips(count).mos)
        jclips(count).mos=[];
    end
    if isnan(jclips(count).stdev)
        jclips(count).stdev=[];
    end
    if isnan(jclips(count).inlsa_mos)
        jclips(count).inlsa_mos=[];
    end
end


% jclipsdata= [ [jclips.test]' [jclips.scene]' [jclips.hrc]' [jclips.file_name]' ...
%     {jclips.loc_start}'   {jclips.loc_stop}'  {jclips.align_start}'  ...
%     {jclips.align_stop}'  {jclips.video_standard}' {jclips.fps}'...
%     {img_size.rows}' {img_size.cols}' {spatial.horizontal}'...
%     {spatial.vertical}'  {scale.horizontal}' {scale.vertical}' ...
%     {jclips.luminance_gain}' , {jclips.luminance_offset}'...
%      {cvr.top}' {cvr.bottom}' {cvr.left}' {cvr.right}' [jclips.subj_system]'...
%      {jclips.mos}' {jclips.stdev}' {jclips.inlsa_mos}' {jclips.viewers}' ...
%      [jclips.hrc_definition]' [jclips.scene_definition]'];

% Changed June06: disp align_start,stop instead of loc_start,stop for orig
% clips.
jclipsdata= [ [jclips.test]' [jclips.scene]' [jclips.hrc]' [jclips.file_name]' ...
    {jclips.loc_start}'   {jclips.loc_stop}'  {jclips.align_start}'  ...
    {jclips.align_stop}'  {jclips.video_standard}' {jclips.fps}'...
    {img_size.rows}' {img_size.cols}' {spatial.horizontal}'...
    {spatial.vertical}'  {scale.horizontal}' {scale.vertical}' ...
    {jclips.luminance_gain}' , {jclips.luminance_offset}'...
     {cvr.top}' {cvr.bottom}' {cvr.left}' {cvr.right}' ...
     {jclips.mos}'  ...
     ];
 
%  modeldata=[ [model_op.clip_name]'  num2cell(model_op.data')];
%  modelhead=['clip name' model_op.par_name];
% This section to make display test, scene, hrc, & fn (in the spreadsheet)for each entry in the
% model_op structure
for cnt = 1:size(jclips,2)
    jc_cn=strcat(jclips(cnt).test,'_', jclips(cnt).scene, '_', jclips(cnt).hrc);
    for mcnt=1:size(model_op.clip_name,2)
        mo_cn=model_op.clip_name(mcnt);
        if strcmp(jc_cn, mo_cn)
            tmo(mcnt).test=jclips(cnt).test;
            tmo(mcnt).scene= jclips(cnt).scene;
            tmo(mcnt).hrc=jclips(cnt).hrc;
            tmo(mcnt).fn = jclips(cnt).file_name;
        end
        mcnt=mcnt+1;
    end
    cnt=cnt+1;
end

% modeldata=[ [model_op.clip_name]'  num2cell(model_op.data')];
modeldata=[ [tmo.test]', [tmo.scene]', [tmo.hrc]', [tmo.fn]',  num2cell(model_op.data')];
% modelhead=['clip name' model_op.par_name];
modelhead=['test', 'scene','hrc', 'filename', model_op.par_name];

unitedheader = strcat('Test: ', jclips(1).test, '    Video: ', jtests.path, '     Date: ', date);

 
% --- Executes on button press in pushbutton8.
% --- Export the graph - using the print function
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Export the Graph-  using the print function

%Let the user choose the file and format

[filename, pathname, filterindex] = uiputfile( ...
{'*.pdf', 'Portable Document Format (*.pdf)';...
'*.jpg;*.jpeg;','JPEG Image (*.jpg, *.jpeg)';...
 '*.bmp','Bitmap (*.bmp)'; ...
 '*.png','Portable Network Graphic (*.png)'; ...
 '*.tiff','TIFF (*.tiff)';...
 '*.ill','Adobe Illustrator (*.ill)';...
 '*.eps','Encapsulated Post Script (*.eps)';}, 'Save as');


%Unfortunately, cannot just export/print the plot; the wholefigure is printed.
%As a workaroudn we use the '-noui' option which hides all of the
%buttons and text. Then we hide the panels and turn the background white.
% Future implementations may want to create a new figure and
%regenerate the graph and then export from there. This should avoid
%exporting a bunch of blank space.

% if filename ~= 0 && pathname ~=0
%if length(filename) > 0 && length(pathname) > 0
if ~(isequal(filename,0) | isequal(pathname,0))
    set(gcf,'PaperPositionMode','auto');
    %make the panels invisible so that you don't see it when you print
    set(handles.uipanel1,'Visible','off');
    %change the background to white
    set(handles.uipanel2,'BackgroundColor',[1 1 1]);
    set(handles.uipanel3,'Visible','off');
    formats={'-dpdf','-djpeg','-dbmp','-dpng','-dtiff','-dill','-depsc'};

    print(formats{filterindex},'-noui',strcat(pathname,filename)); %print while hiding the gui

    %now make the panels visible again
    set(handles.uipanel1,'Visible','on');
    set(handles.uipanel2,'BackgroundColor',[.831 .816 .784]);
    set(handles.uipanel3,'Visible','on');
else
    return
end


% --- Executes on selection change in listbox2.
%function listbox2_Callback(hObject, eventdata, handles, model_op)
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2

working_dir=get(handles.text9, 'String');
lfile=fullfile(working_dir, 'bvqm-status');

load (lfile); % bvqm-status;

%This listbox consists of a lists of clips to choose from.
refresh_scene_detail(hObject, eventdata, handles,model_op);


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
%This listbox consists of a lists of clips to choose from.
%Let's populate the box with clips.



function text1_Callback(hObject, eventdata, handles)
% hObject    handle to text1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text1 as text
%        str2double(get(hObject,'String')) returns contents of text1 as a double


%Export to Excel
% --- Executes on button press in pushbutton12.
%function pushbutton12_Callback(hObject, eventdata, handles,model_selected, jclips, jtests, model_op, cal_type)
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


working_dir=get(handles.text9, 'String');
lfile=fullfile(working_dir, 'bvqm-status');

clk=fix(clock);
hr=sprintf('%02.0f',clk(4)); mn=sprintf('%02.0f',clk(5)); sc=sprintf('%02.0f',clk(6));
time=strcat(hr,'_',mn);

load (lfile); % bvqm-status;

directoryname=uigetdir('','Pick a directory to save to');
if directoryname ~= 0

    [jclipsdata  jclipshead modeldata modelhead unitedheader]= prepexport(jclips, jtests, model_op, cal_type,hObject, eventdata, handles);
    if ispc ; path_sep = '\'; else; path_sep='/'; end

    % workin in linux or computer without Excel:
    jclips_file_out=strcat(directoryname,path_sep,char(jtests.name),'_clips_',date,'@',time, '.csv');
    model_file_out=strcat(directoryname,path_sep,char(jtests.name),'_model_',date,'@', time, '.csv');
    calib_file_out=strcat(directoryname,path_sep,char(jtests.name),'_calibration_',date,'@',time);

    make_csv(jclips_file_out, jclipsdata,unitedheader, jclipshead);
    make_csv(model_file_out, modeldata,unitedheader, modelhead);
    [error_status] = bvqm_pc_calexport(jclips, jtests, calib_file_out);
    vfd_file_in = strcat(working_dir,path_sep,'vfd_results.csv');
    vfd_file_out = strcat(directoryname,path_sep,char(jtests.name),'_model_',date,'@', time, '_vfd.csv');
    if exist(vfd_file_in,'file'),
        copyfile(vfd_file_in,vfd_file_out);
        vfd_string = strcat('Model VFD data: ', char(jtests.name),'_model_',date,'@', time, '_vfd.csv');
    else
        vfd_string = '';
    end
    if ~error_status,
        msgbox( 'Data could not be exported to CSV.', 'Export Failed','Error');
    else
        msg='CSV files may be viewed using any spreadsheet application that can read comma delimited data.';
        message= {'CSV files saved:';''; 'Clip data:';jclips_file_out;'';'Model data:'; model_file_out; '';'Manual Calibration Output files:'; calib_file_out; '_sheet1.csv '; calib_file_out; '_sheet2.csv'; calib_file_out; '_sheet3.csv '; '' ; vfd_string ;'' ; msg};
        msgbox( message,'Files written', 'help');
    end

end

function make_csv(filename, data, title, headings)

% note: jclipsdata=data, jclipshead=headings, unitedheader=title

%disp('Started writing CSV ...');

% Open the file for writing
fid = fopen(filename, 'w');

% write out the title 
fprintf(fid, '%s\n\n', title{1});

% write out headers
for i = 1:length(headings)
    fprintf(fid, [headings{i}, ',']);
end;

% erase the last comma and go to the next line
fseek(fid, -1, 'eof');
fprintf(fid, '\n');

% write out the data
for (i = 1:size(data,1))
    for (j = 1:size(data,2))

        % write out one data element and place a comma
        a=data{i,j};
        if ischar(a)
            fprintf(fid, [a, ',']);
        else
            fprintf(fid, [mat2str(a), ',']);
        end
    end

    % erase the last comma and go to the next line
    fseek(fid, -1, 'eof');
    fprintf(fid, '\n');
end;

% erase the last comma and go to the next line
fseek(fid, -1, 'eof');
fprintf(fid, '\n');

fclose(fid);

%disp('Done writing CSV');

