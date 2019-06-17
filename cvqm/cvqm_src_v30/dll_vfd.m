function [results results_fuzzy] = dll_vfd(delay, varargin)
% VFD options
%
%   Estimate the variable frame delays (VFD) of all frames (or fields) in
%   the original and processed video clips already stored in dll_video.
%   These files must be in the UYVY AVI file format.  VFD results are
%   output in a row vector.
%
%   This routine does not perform any spatial registration so the original
%   and processed clips are assumed to have been spatially registered
%   beforehand.  The original and processed video files can have a
%   different number of frames but ideally, there should be a matching
%   original frame (or field) for every processed frame (or field).
%
% SYNTAX
%   results = dll_vfd(options)
%
% DESCRIPTION
%    For each frame (or field) in the processed file, this algorithm finds
%    the frame (or field) in the original file that minimizes the mean
%    squared error, subject to the constraints imposed by the optional
%    arguments.
%
%   Any or all of the following optional properties may be requested.
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
%
% RESULTS
%   The output results are in an array that includes the indices of the
%   frames (or fields) of the original file that is the best fit for each
%   processed file frame (or field) in sequential order.  These indices
%   start at the beginning of the processed section of the clip, not the
%   beginning of the entire clip.
%
%   These results may be processed by the dll_vfd_print function, along
%   with the temporal delay found by a previous calibration algorithm, to
%   output a file with the indices that correlate to the entire video file.
%
%
% EXAMPLES
%   vfd();
%   vfd('t_uncert', 30, 'reframe', 'causal');
%

aba_t = 8;

% Validate the scan_type
scan_type = dll_video('get_video_standard', 1);
if strcmp(scan_type, 'interlace_lower_field_first');
    scan_type = 'interlaced_lff';
elseif strcmp(scan_type, 'interlace_upper_field_first');
    scan_type = 'interlaced_uff';
end
if (~strcmpi(scan_type,'progressive') && ~strcmpi(scan_type,'interlaced_lff') && ~strcmpi(scan_type,'interlaced_uff'))
    error('Invalid scan_type');
end

% Set rewind points
dll_video('set_rewind', 1);
dll_video('set_rewind', 2);

% Validate input arguments and set their defaults
t_uncert = 30; %ceil(dll_video('fps', 1));
reframe = 0;
causal = 0;
cnt=1;
while cnt <= length(varargin),
    if ~ischar(varargin{cnt}),
        error('Property value passed into vfd is not recognized');
    end
    if strcmpi(varargin{cnt},'t_uncert') == 1
        t_uncert = varargin{cnt+1};
        cnt = cnt + 2;
    elseif strcmpi(varargin{cnt},'reframe') == 1
        reframe = 1;
        cnt = cnt + 1;
    elseif strcmpi(varargin{cnt},'causal') == 1
        causal = 1;
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
    
% Get processed file information
[rows, cols] = dll_video('size', 2);

% Get original file information
[rows_orig, cols_orig] = dll_video('size', 1);

% Verify that the processed and original are the same resolution
if (rows ~= rows_orig || cols ~= cols_orig)
    error('Processed and original files have different image resolutions.');
end

% Generate the length of the video files to be read (based on the shorter
% time slice)
tlength_original = dll_calib_video('total_sec', 1);
tlength_processed = dll_calib_video('total_sec', 2);
if delay < 0 % 8/5/11 if delay is -, take frames off the end of the processed
    tlength_processed = tlength_processed + delay/dll_video('fps', 2);
end

% Read the video files
dll_video('set_tslice', 1, tlength_original);
dll_video('set_tslice', 2, tlength_processed);
yo = dll_calib_video('tslice', 1);
yp = dll_calib_video('tslice', 2);

% Generate the function call with all the desired options
func_call = 'est_var_frame_delays(yp,yo,''normalize'',';  % always use the normalize option as it seems to work the best

if (reframe)
    func_call = strcat(func_call,'''reframe'',');
end

if (causal)
    func_call = strcat(func_call,'''causal'',');
end

if (strcmpi(scan_type,'interlaced_lff'))
    func_call = strcat(func_call,'''interlaced'',1,');
end

if (strcmpi(scan_type,'interlaced_uff'))
    func_call = strcat(func_call,'''interlaced'',2,');
end

func_call = strcat(func_call,'''t_uncert'',',num2str(t_uncert),')');

% Call the est_var_frame_delays function to get results
[results results_rmse results_fuzzy results_fuzzy_mse] = eval(func_call);

%  If the VFD alignment algorithm failed, then results==0.                 
%  In that case, use the default time alignment.                           
[~, ~, nframes] = size(yp);
new_first_align = 1;
if (results == 0)                                                          
    if (~strcmpi(scan_type,'progressive'))  % interlaced                         
        npts = nframes*2;                                                       
    else  % progressive                                                        
        npts = nframes;                                                          
    end                                                                          
    results = new_first_align:new_first_align+npts-1;                            
    fprintf('VFD algorithm failed, using default time alignment.\n');
else
    aba = ave_best_aligned(results_fuzzy_mse);  % results_fuzzy_mse is loaded from mat file
    if (aba > aba_t)
        if (is_interlaced)
            npts = nframes*2;
        else
            npts = nframes;
        end
        results = new_first_align:new_first_align+npts-1;  % Overwrite the VFD results with gclips alignment
    end
end                                                              

% Return to rewind points
dll_video('rewind', 1);
dll_video('rewind', 2);

end
