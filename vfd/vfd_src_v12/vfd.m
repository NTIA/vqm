function vfd(proc_file, orig_file, scan_type, results_file, varargin)
% VFD 'proc_file' 'orig_file' 'scan_type' 'results_file' options
%
%   Estimate the variable frame delays (VFD) of all frames (or fields) in
%   the user specified processed video file ('proc_file') that is in either
%   uncompressed UYVY AVI file format (default) or raw big-YUV file format
%   (optional).  The original video file is given by 'orig_file'.  The user
%   must specify the 'scan_type' of the video files as either
%   'progressive', 'interlaced_uff' (interlaced upper field first), or
%   'interlaced_lff' (interlaced lower field first), since this information
%   is not available in the AVI file format.  VFD results are appended to 
%   the file 'results_file' (see the RESULTS section below for a complete
%   description of the results output file).  
%
%   This routine does not perform any spatial registration so the original
%   and processed clips are assumed to have been spatially registered
%   beforehand.  The original and processed video files can have a
%   different number of frames but ideally, there should be a matching
%   original frame (or field) for every processed frame (or field).
%
% SYNTAX
%   vfd 'proc_file' 'orig_file' 'scan_type' 'results_file' options
%
% DESCRIPTION
%    For each frame (or field) in the processed file, this algorithm finds
%    the frame (or field) in the original file that minimizes the mean
%    squared error, subject to the constraints imposed by the optional
%    arguments.  See the RESULTS section below for the results output file
%    format.
%
%   Any or all of the following optional properties may be requested (the
%   'yuv' option is required for yuv files, but not for avi files since
%   this information is read from the avi header).
%
%   'yuv' rows cols        Specifies the number of rows and cols of the 
%                          big-YUV files. 
%
%   'sroi' top left bottom right    Only use the specified spatial region 
%                                   of interest (sroi) for the vfd
%                                   calculation.  By default, all of the
%                                   image is used.  The sroi is inclusive,
%                                   where top/left start at 1. 
%
%   'troi' fstart fstop    Only calculate vfd for the specified temporal
%                          region of interest (troi) of the processed video
%                          clip, where fstart and fstop are included and
%                          given in frames.  By default, the vfd of the
%                          entire processed video clip is calculated
%                          (fstart = 1, fstop = number of frames in file).
%
%   'first_align' a           Specifies the best guess for the original
%                             file frame number that corresponds to the
%                             first frame in the processed troi.  By
%                             default, this is set to fstart, which is set
%                             to 1 when troi is not specified.
%
%   't_uncert' t              Specifies the temporal uncertainty (plus or
%                             minus t frames) over which to search.  The
%                             processed remains fixed and the original is
%                             shifted.  The center (zero shift) point for
%                             the first frame (or field) is given by
%                             first_align.  By default, temporal
%                             uncertainty is set to 30 frames.  It can have
%                             a minimum value of 1 frame.  When the
%                             original cannot be shifted by the temporal 
%                             uncertainty (e.g., perhaps near the ends of 
%                             the sequence), the original will be shifted
%                             up to the maximum extent possible.
%
%   'reframe'  Allow for the possibility that the processed video clip has
%              been reframing.  This option is only valid for a scan_type
%              of 'interlaced_uff' or 'interlaced_lff'.  Reframing can vary
%              throughout the processed clip, although this should be rare.
%              This option will increase the runtime substantially since
%              extra spatial shifts must be examined, but it should be used
%              if there is any possibility of reframing existing in the
%              processed video clip.  See Section 3.1.2 of NTIA Report
%              TR-02-392 for a definition of reframing.
%
%   'causal'   Impose causality constraint so that later frames (fields) in
%              the processed clip cannot align to original frames (fields)
%              that are earlier in time than found for the proceeding
%              processed frames (fields).  For interlaced video, a
%              one-field jump back in time is allowed since this is
%              indicative of a frozen frame.  By default, causality is
%              turned off (yes, codecs can output non-causal sequences).
%              But specifying the causal option is usually recommended. 
%
%   'verbose'   Display output and plots during processing.
%
% RESULTS
%   The output results file is in Comma-Separated Values (CSV format):
%   
%   'proc_file', proc_indices
%   'orig_file', orig_indices
%
%   proc_indices and orig_indices are row vectors of length equal to the
%   number of frames in the processed temporal region of interest (for
%   scan_type = 'progressive') or of length equal to 2 * the number of 
%   frames in the processed temporal region of interest (for scan_type =
%   'interlaced_uff' or 'interlaced_lff').  For each processed frame (or
%   field) index, the best matching original frame (or field) index is
%   given, where 1 is the first frame (or field) in the file.
%
% EXAMPLES
%   vdf 'proc.avi' 'orig.avi' 'progressive' 'results.csv'
%   vfd 'proc.yuv' 'orig.yuv' 'interlaced_lff' 'results.csv' 'yuv' 486 720 'reframe' 'causal' 'verbose' 
%

if nargin == 0,
    fprintf(' vfd ''proc_file'' ''orig_file'' ''scan_type'' ''results_file'' options\n');
    fprintf('\n');
    fprintf('   Estimate the variable frame delays (VFD) of all frames (or fields) in\n');
    fprintf('   the user specified processed video file (''proc_file'') that is in either\n');
    fprintf('   uncompressed UYVY AVI file format (default) or raw big-YUV file format\n');
    fprintf('   (optional).  The original video file is given by ''orig_file''.  The user\n');
    fprintf('   must specify the ''scan_type'' of the video files as either\n');
    fprintf('   ''progressive'', ''interlaced_uff'' (interlaced upper field first), or\n');
    fprintf('   ''interlaced_lff'' (interlaced lower field first), since this information\n');
    fprintf('   is not available in the AVI file format.  VFD results are appended to \n');
    fprintf('   the file ''results_file'' (see the RESULTS section below for a complete\n');
    fprintf('   description of the results output file).  \n');
    fprintf('\n');
    fprintf('   This routine does not perform any spatial registration so the original\n');
    fprintf('   and processed clips are assumed to have been spatially registered\n');
    fprintf('   beforehand.  The original and processed video files can have a\n');
    fprintf('   different number of frames but ideally, there should be a matching\n');
    fprintf('   original frame (or field) for every processed frame (or field).\n');
    fprintf('\n');
    fprintf(' SYNTAX\n');
    fprintf('   vfd ''proc_file'' ''orig_file'' ''scan_type'' ''results_file'' options\n');
    fprintf('\n');
    fprintf(' DESCRIPTION\n');
    fprintf('    For each frame (or field) in the processed file, this algorithm finds\n');
    fprintf('    the frame (or field) in the original file that minimizes the mean\n');
    fprintf('    squared error, subject to the constraints imposed by the optional\n');
    fprintf('    arguments.  See the RESULTS section below for the results output file\n');
    fprintf('    format.\n');
    fprintf('\n');
    fprintf('   Any or all of the following optional properties may be requested (the\n');
    fprintf('   ''yuv'' option is required for yuv files, but not for avi files since\n');
    fprintf('   this information is read from the avi header).\n');
    fprintf('\n');
    fprintf('   ''yuv'' rows cols        Specifies the number of rows and cols of the \n');
    fprintf('                          big-YUV files. \n');
    fprintf('\n');
    fprintf('   ''sroi'' top left bottom right    Only use the specified spatial region \n');
    fprintf('                                   of interest (sroi) for the vfd\n');
    fprintf('                                   calculation.  By default, all of the\n');
    fprintf('                                   image is used.  The sroi is inclusive,\n');
    fprintf('                                   where top/left start at 1. \n');
    fprintf('\n');
    fprintf('   ''troi'' fstart fstop    Only calculate vfd for the specified temporal\n');
    fprintf('                          region of interest (troi) of the processed video\n');
    fprintf('                          clip, where fstart and fstop are included and\n');
    fprintf('                          given in frames.  By default, the vfd of the\n');
    fprintf('                          entire processed video clip is calculated\n');
    fprintf('                          (fstart = 1, fstop = number of frames in file).\n');
    fprintf('\n');
    fprintf('   ''first_align'' a           Specifies the best guess for the original\n');
    fprintf('                             file frame number that corresponds to the\n');
    fprintf('                             first frame in the processed troi.  By\n');
    fprintf('                             default, this is set to fstart, which is set\n');
    fprintf('                             to 1 when troi is not specified.\n');
    fprintf('\n');
    fprintf('   ''t_uncert'' t              Specifies the temporal uncertainty (plus or\n');
    fprintf('                             minus t frames) over which to search.  The\n');
    fprintf('                             processed remains fixed and the original is\n');
    fprintf('                             shifted.  The center (zero shift) point for\n');
    fprintf('                             the first frame (or field) is given by\n');
    fprintf('                             first_align.  By default, temporal\n');
    fprintf('                             uncertainty is set to 30 frames.  It can have\n');
    fprintf('                             a minimum value of 1 frame.  When the\n');
    fprintf('                             original cannot be shifted by the temporal \n');
    fprintf('                             uncertainty (e.g., perhaps near the ends of \n');
    fprintf('                             the sequence), the original will be shifted\n');
    fprintf('                             up to the maximum extent possible.\n');
    fprintf('\n');
    fprintf('   ''reframe''  Allow for the possibility that the processed video clip has\n');
    fprintf('              been reframing.  This option is only valid for a scan_type\n');
    fprintf('              of ''interlaced_uff'' or ''interlaced_lff''.  Reframing can vary\n');
    fprintf('              throughout the processed clip, although this should be rare.\n');
    fprintf('              This option will increase the runtime substantially since\n');
    fprintf('              extra spatial shifts must be examined, but it should be used\n');
    fprintf('              if there is any possiblity of reframing existing in the\n');
    fprintf('              processed video clip.  See Section 3.1.2 of NTIA Report\n');
    fprintf('              TR-02-392 for a definition of reframing.\n');
    fprintf('\n');
    fprintf('   ''causal''   Impose causality constraint so that later frames (fields) in\n');
    fprintf('              the processed clip cannot align to original frames (fields)\n');
    fprintf('              that are earlier in time than found for the proceeding\n');
    fprintf('              processed frames (fields).  For interlaced video, a\n');
    fprintf('              one-field jump back in time is allowed since this is\n');
    fprintf('              indicative of a frozen frame.  By default, causality is\n');
    fprintf('              turned off (yes, codecs can output non-causal sequences).\n');
    fprintf('              But specifying the causal option is usually recommended. \n');
    fprintf('\n');
    fprintf('   ''verbose''   Display output and plots during processing.\n');
    fprintf('\n');
    fprintf(' RESULTS\n');
    fprintf('   The output results file is in Comma-Separated Values (CSV format):\n');
    fprintf('\n');
    fprintf('   ''proc_file'', proc_indices\n');
    fprintf('   ''orig_file'', orig_indices\n');
    fprintf('\n');
    fprintf('   proc_indices and orig_indices are row vectors of length equal to the\n');
    fprintf('   number of frames in the processed temporal region of interest (for\n');
    fprintf('   scan_type = ''progressive'') or of length equal to 2 * the number of \n');
    fprintf('   frames in the processed temporal region of interest (for scan_type =\n');
    fprintf('   ''interlaced_uff'' or ''interlaced_lff'').  For each processed frame (or\n');
    fprintf('   field) index, the best matching original frame (or field) index is\n');
    fprintf('   given, where 1 is the first frame (or field) in the file.\n');
    fprintf('\n');
    fprintf(' EXAMPLES\n');
    fprintf('   vdf ''proc.avi'' ''orig.avi'' ''progressive'' ''results.csv''\n');
    fprintf('   vfd ''proc.yuv'' ''orig.yuv'' ''interlaced_lff'' ''results.csv'' ''yuv'' 486 720 ''reframe'' ''causal'' ''verbose'' \n');
    fprintf('\n');
    return;
end

% strip off the extra single quotes ''
proc_file = eval(proc_file);
orig_file = eval(orig_file);
scan_type = eval(scan_type);
results_file = eval(results_file);

% Validate the scan_type
if (~strcmpi(scan_type,'progressive') && ~strcmpi(scan_type,'interlaced_lff') && ~strcmpi(scan_type,'interlaced_uff'))
    error('Invalid scan_type');
end

% Validate input arguments and set their defaults
is_yuv = 0;  % default file type, uncompressed UYVY AVI
is_whole_image = 1;
is_whole_time = 1;
first_align = 0;
t_uncert = 30;
reframe = 0;
causal = 0;
verbose = 0;
cnt=1;
while cnt <= length(varargin),
    if ~ischar(varargin{cnt}),
        error('Property value passed into vfd is not recognized');
    end
    if strcmpi(eval(char(varargin(cnt))),'yuv') == 1
        rows = str2double(varargin{cnt+1});
        cols = str2double(varargin{cnt+2});
        is_yuv = 1;
        cnt = cnt + 3;
    elseif strcmpi(eval(char(varargin(cnt))),'sroi') == 1
        top = str2double(varargin{cnt+1});
        left = str2double(varargin{cnt+2});
        bottom = str2double(varargin{cnt+3});
        right = str2double(varargin{cnt+4});
        is_whole_image = 0;
        cnt = cnt + 5;
    elseif strcmpi(eval(char(varargin(cnt))),'troi') == 1
        fstart = str2double(varargin{cnt+1});
        fstop = str2double(varargin{cnt+2});
        is_whole_time = 0;
        cnt = cnt + 3;
    elseif strcmpi(eval(char(varargin(cnt))),'first_align') == 1
        first_align = str2double(varargin{cnt+1});
        cnt = cnt + 2;
    elseif strcmpi(eval(char(varargin(cnt))),'t_uncert') == 1
        t_uncert = str2double(varargin{cnt+1});
        cnt = cnt + 2;
    elseif strcmpi(eval(char(varargin(cnt))),'reframe') == 1
        reframe = 1;
        cnt = cnt + 1;
    elseif strcmpi(eval(char(varargin(cnt))),'causal') == 1
        causal = 1;
        cnt = cnt + 1;
    elseif strcmpi(eval(char(varargin(cnt))),'verbose') == 1
        verbose = 1;
        cnt = cnt + 1;
    else
        error('Property value passed into vfd not recognized');
    end
end

% Validate reframing option
if (reframe && strcmpi(scan_type,'progressive'))
    error('Reframe option not allowed for progressive video');
end

% Get the processed and original file information
if (~is_yuv)  % AVI file
    
    % Get processed file information
    [avi_info] = read_avi('Info',proc_file);
    rows = avi_info.Height;
    cols = avi_info.Width;
    tframes = avi_info.NumFrames;  % total frames in processed file

    
    % Get original file information
    [avi_info_orig] = read_avi('Info',orig_file);
    rows_orig = avi_info_orig.Height;
    cols_orig = avi_info_orig.Width;
    tframes_orig = avi_info_orig.NumFrames;
    
else  % big-YUV file
    
    % Get the processed file information
    [fid, message] = fopen(proc_file, 'r');
    if fid == -1
        fprintf(message);
        error('Cannot open processed big-YUV file %s', proc_file);
    end
    % Find last frame.
    fseek(fid,0, 'eof');
    tframes = ftell(fid) / (2 * rows * cols);
    fclose(fid);
    
    % Get the original file information
    rows_orig = rows;
    cols_orig = cols;
    [fid, message] = fopen(orig_file, 'r');
    if fid == -1
        fprintf(message);
        error('Cannot open original big-YUV file %s', orig_file);
    end
    % Find last frame.
    fseek(fid,0, 'eof');
    tframes_orig = ftell(fid) / (2 * rows * cols);
    fclose(fid);
    
end

% Verify that the processed and original are the same resolution
if (rows ~= rows_orig || cols ~= cols_orig)
    error('Processed and original files have different image resolutions.');
end

% Set/Validate the SROI
if (is_whole_image)
    top = 1;
    left = 1;
    bottom = rows;
    right = cols;
elseif (top<1 || left<1 || bottom>rows || right>cols || top>bottom || left>right)
    error('Invalid spatial region of interest (SROI) for the processed file.');
end

% Set/Validate the TROI of the processed file
if (is_whole_time)
    fstart= 1;
    fstop = tframes;
elseif (fstart<1 || fstop>tframes || fstart>fstop)
    error('Invalid temporal region of interest (TROI) for the processed file.');
end

%  Assign the original first alignment point and validate
if (~first_align)  % a value for first_align was not input by the user so assign default
    first_align = fstart;
end
if (first_align<1 || first_align>tframes_orig)
    error('Invalid first_align for the original file.');
end

%  Validate t_uncert
if (t_uncert<1 || t_uncert>tframes_orig)
    error('Invalid search uncertainty t_uncert.');
end

%  Find the original time segment to read
offset_orig = first_align-fstart;  % search offset of the orig file with respect to the proc file
fstop_orig = min(tframes_orig, fstop+offset_orig+t_uncert);  % the last original frame to read

% Read extra orig frames at the beginning to allow for search range.
% Calculate a new first_align point that is referenced to the first
% original frame that is read rather than to the first frame in the original file.
fstart_orig = max(1, fstart+offset_orig-t_uncert);
new_first_align = first_align-fstart_orig+1;  % no change if I start reading the first frame of orig

% Read in the original and processed S-T segments
if (~is_yuv)  % AVI file
    
    % Read in video and clear color planes to free up memory
    [yp, cb, cr] = read_avi('YCbCr',proc_file,'frames',fstart,fstop,'sroi',top,left,bottom,right);
    clear cb cr;
    [yo, cb, cr] = read_avi('YCbCr',orig_file,'frames',fstart_orig,fstop_orig,'sroi',top,left,bottom,right);
    clear cb cr;
    
else  % big-YUV file
    
    % Read in video and clear color planes to free up memory
    [yp] = read_bigyuv(proc_file,'frames',fstart,fstop,'size',rows,cols,'sroi',top,left,bottom,right);
    [yo] = read_bigyuv(orig_file,'frames',fstart_orig,fstop_orig,'size',rows,cols,'sroi',top,left,bottom,right);
    
end

% Generate the function call with all the desired options
func_call = 'est_var_frame_delays(yp,yo,''normalize'',';  % always use the normalize option as it seems to work the best

if (reframe)
    func_call = strcat(func_call,'''reframe'',');
end

if (causal)
    func_call = strcat(func_call,'''causal'',');
end

if (verbose)
    func_call = strcat(func_call,'''verbose'',');
end

if (strcmpi(scan_type,'interlaced_lff'))
    func_call = strcat(func_call,'''interlaced'',1,');
    new_first_align = 2*new_first_align-1;  % convert to fields for est_var_frame_delays
end

if (strcmpi(scan_type,'interlaced_uff'))
    func_call = strcat(func_call,'''interlaced'',2,');
    new_first_align = 2*new_first_align-1;  % convert to fields for est_var_frame_delays
end

func_call = strcat(func_call,'''first_align'',',num2str(new_first_align),',','''t_uncert'',',num2str(t_uncert),')');

% Call the est_var_frame_delays function to get results
[results results_rmse results_fuzzy results_fuzzy_mse] = eval(func_call);

%  If the VFD alignment algorithm failed, then results==0.
%  In that case, use the default time alignment.
[nrows, ncols, nframes] = size(yp);
if (results == 0)
    if (~strcmpi(scan_type,'progressive'))  % interlaced
        npts = nframes*2;
    else  % progressive
        npts = nframes;
    end
    results = new_first_align:new_first_align+npts-1;
    fprintf('VFD algorithm failed, using default time alignment.\n');
end

% Translate results to use the orig and proc file indexing
if (strcmpi(scan_type,'progressive'))
    proc_indices = (fstart-1) + (1:length(results));
    orig_indices = (fstart_orig-1) + results;
else % interlaced
    proc_indices = 2*(fstart-1) + (1:length(results));
    orig_indices = 2*(fstart_orig-1) + results;
end

% Save results
fid_results = fopen(results_file,'a');  % open results file for appending

if (strcmpi(scan_type,'progressive'))
    fprintf(fid_results,'File Name, Matching Frame Indices\n');
else % interlaced
    fprintf(fid_results,'File Name, Matching Field Indices\n');
end

npts = length(proc_indices);

fprintf(fid_results,'%s, ',proc_file);
for i = 1:npts-1
    fprintf(fid_results,'%f, ',proc_indices(i));
end
fprintf(fid_results,'%f\n',proc_indices(i+1));

fprintf(fid_results,'%s, ',orig_file);
for i = 1:npts-1
    fprintf(fid_results,'%f, ',orig_indices(i));
end
fprintf(fid_results,'%f\n',orig_indices(i+1));
fclose(fid_results);

close all;

end
