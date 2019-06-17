function psnr_vfd(clip_dir, test, scan_type, psnr_file, varargin)
% PSNR_VFD 'clip_dir' 'test' 'scan_type' 'psnr_file' options
%   Estimate the Variable Frame Delay (VFD) Y-channel PSNR (PSNR) of all
%   clips and HRCs (Hypothetical Reference Circuits) in a video test (input
%   argument "test") where the video clips are stored in the specified
%   directory ("clip_dir").  The video clips must have names that conform to
%   the standard naming conventions test_scene_hrc.avi (or optionally,
%   test_scene_hrc.yuv) with no extra '_' or '.' in the file names.  "test"
%   is the name of the test, "scene" is the name of the scene, and "hrc" is
%   the name of the HRC.  The name of the original reference clip for the
%   PSNR calculation must be "test_scene_original.avi".
%   
%   The user must specify the 'scan_type' of the video files as either
%   'progressive', 'interlaced_uff' (interlaced upper field first), or
%   'interlaced_lff' (interlaced lower field first), since this information
%   is not available in the AVI or the Big YUV file formats.
%   
%   PSNR_VFD uses the results file output by the PSNR_SEARCH program 
%   ("psnr_file") as a calibration starting point (i.e., Yshift, Xshift,
%   Tshift, Gain, Offset), but goes one step further by applying VFD
%   matching of the original video stream to the processed video stream
%   (i.e., the original video stream is modified so that it matches the
%   processed video stream frame-by-frame, or field-by-field for interlaced
%   video).  Output results from the PSNR_VFD program are stored in a file
%   with the same root name as "psnr_file" (i.e., name without the file
%   extension), appended with "_vfd.csv" (for Comma Separated Value).
%
% SYNTAX
%   psnr_vfd 'clip_dir' 'test' 'scan_type' 'psnr_file' options
%
% DESCRIPTION
%   This program uses the clip calibration information (i.e., Yshift,
%   Xshift, Tshift , Gain, Offset) from PSNR_SEARCH (in "psnr_file").  The 
%   original clip is shifted by (Yshift, Xshift, Tshift) with respect to
%   the processed clip.  The Y-image of the processed clip is multiplied
%   by Gain and then the Offset is added.  For speed, all calibration is
%   held constant for the VFD estimation.  The VFD estimation algorithm is
%   applied to find the best matching original frame (or field) for each
%   processed frame (or field).  Then, the original is VFD-matched to the 
%   processed.  A final gain and offset correction (called Gain_Adjust and
%   Offset_Adjust) is applied to the processed video (i.e.,
%   Gain_Adjust*y_proc + Offset_Adjust) before calculation of the final
%   VFD-corrected PSNR (PSNR_VFD).  Two perception-based VFD parameters are
%   also calculated (Par1_VFD and Par2_VFD) that attempt to capture the
%   perceptual distortions in the flow of scene motion (since these
%   distortions are removed from PSNR_VFD).  Par1_VFD is extracted from
%   only the VFD information while Par2_VFD uses both the VFD information
%   and the motion in the processed video clip.  For a description of these
%   two parameters, see the 2011 NTIA Technical Memorandum (TM) entitled
%   "Variable Frame Delay (VFD) Parameters for Video Quality Measurements."
%
%   The above procedure is repeated for each processed clip in the
%   "clip_dir" that belongs to the video test specified by "test". 
%
%   A peak signal of 255 is used for calculation of PSNR.  Double precision
%   calculations are used everywhere.  A 64-bit operating system with at
%   least 4 GB of free memory is recommended since the entire double
%   precision versions of the original and processed sequences must be held
%   in memory.
%
%   Any or all of the following optional properties may be requested (the
%   first option is required for yuv files, not avi files).
%
%   'yuv' rows cols    Specifies the number of rows and cols for the Big
%                      YUV files (if using Big YUV files).  The default is
%                      to assume AVI files (*.avi).  Big YUV format is a
%                      binary format for storing ITU-R Recommendation
%                      BT.601 video sequences.  The format can be used for
%                      any image size.  In the Big YUV format, all the
%                      frames are stored sequentially in one big binary
%                      file. The sampling is 4:2:2 and image pixels are
%                      stored sequentially by video scan line as bytes in
%                      the following order: Cb1 Y1 Cr1 Y2 Cb3 Y3 Cr3 Y4…,
%                      where Y is the luminance component, Cb is the blue
%                      chrominance component, Cr is the red chrominance
%                      component, and the subscript is the pixel number.
%
%   'sroi' top left bottom right    Only use the specified spatial region 
%                                   of interest (sroi) of the processed
%                                   video clip for all calculations. The
%                                   sroi is inclusive, where top/left start
%                                   at 1. By default, sroi is the entire
%                                   image reduced by the calibration shift
%                                   (Yshift, Xshift).  If the user inputs a
%                                   sroi, allowance must be made for the 
%                                   maximum spatial shift encountered in
%                                   the PSNR_SEARCH results ("psnr_file").
%                                   For interlaced video, top must be odd
%                                   while bottom must be even.
%
%   'troi' fstart fstop    Only calculate PSNR_VFD for the specified
%                          temporal region of interest (troi) of the
%                          processed video clip, where fstart and fstop are
%                          included and given in frames.  By default, the
%                          troi is the entire file reduced by the temporal
%                          calibration shift (Tshift).  If the user inputs
%                          an fstart and fstop, allowance must be made for
%                          the maximum temporal shift encountered in 
%                          "psnr_file".
%
%   't_uncert' t    Specifies the temporal uncertainty (plus or minus t
%                   frames) over which to perform the VFD search.  The
%                   processed remains fixed and the original is shifted.
%                   The center (zero shift) point for the temporal search
%                   assumes the temporal alignment given by Tshift in the 
%                   "psnr_file" results from PSNR_SEARCH. By default,
%                   temporal uncertainty is set to 30 frames.  It can have
%                   a minimum value of 1 frame.  When the original cannot
%                   be shifted by the temporal uncertainty (e.g., perhaps
%                   near the ends of the sequence), the original will be
%                   shifted up to the maximum extent possible.  To
%                   accomodate a full temporal search at the beginning and
%                   end of the sequence, increase fstart and decrease fstop
%                   in the processed troi accordingly.
%
%   'reframe'  Allow for the possibility that the processed video clip has
%              been reframing.  This option is only valid for a scan_type
%              of 'interlaced_uff' or 'interlaced_lff'.  Reframing can vary
%              throughout the processed clip, although this should be rare.
%              This option will increase the runtime substantially since
%              extra spatial shifts must be examined, but it should be used
%              if there is any possibility of time varying reframing
%              existing in the processed video clip.  See Section 3.1.2 of
%              NTIA Report TR-02-392 for a definition of reframing.  For
%              constant reframing that was properly detected by the
%              PSNR_SEARCH program, Tshift will have a half frame (0.5)
%              Tshift.  This condition will be detected/corrected by 
%              PSNR_VFD (the original video is reframed accordingly).
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
%   'verbose'   Display intermediate progress during processing.
%
% EXAMPLES
%   These three examples illustrate the Big YUV format for QCIF, CIF, and
%   VGA. 
%
%   psnr_vfd 'd:\q01\' 'q01' 'progressive' 'q01_psnr.csv' 'yuv' 144 176 'causal'
%   psnr_vfd 'd:\c01\' 'c01' 'progressive' 'c01_psnr.csv' 'yuv' 288 352 'causal'
%   psnr_vfd 'd:\v01\' 'v01' 'progressive' 'v01_psnr.csv' 'yuv' 480 640 'causal'
%
%   This example illustrates uncompressed UYVY AVI format for 525-line video.
%
%   psnr_vfd 'd:\rr525\' 'rr525' 'interlaced_lff' 'rr525_psnr.csv' 'causal' 'reframe'
%

% This prints out help if the user runs with no command line arguments
if nargin == 0,
    fprintf('PSNR_VFD ''clip_dir'' ''test'' ''scan_type'' ''psnr_file'' options\n');
    fprintf('\n');
    fprintf('  Estimate the Variable Frame Delay (VFD) Y-channel PSNR (PSNR) of all\n');
    fprintf('  clips and HRCs (Hypothetical Reference Circuits) in a video test (input\n');
    fprintf('  argument "test") where the video clips are stored in the specified\n');
    fprintf('  directory ("clip_dir").  The video clips must have names that conform to\n');
    fprintf('  the standard naming conventions test_scene_hrc.avi (or optionally,\n');
    fprintf('  test_scene_hrc.yuv) with no extra ''_'' or ''.'' in the file names.  "test"\n');
    fprintf('  is the name of the test, "scene" is the name of the scene, and "hrc" is\n');
    fprintf('  the name of the HRC.  The name of the original reference clip for the\n');
    fprintf('  PSNR calculation must be "test_scene_original.avi".\n');
    fprintf('\n');
    fprintf('  The user must specify the ''scan_type'' of the video files as either\n');
    fprintf('  ''progressive'', ''interlaced_uff'' (interlaced upper field first), or\n');
    fprintf('  ''interlaced_lff'' (interlaced lower field first), since this information\n');
    fprintf('  is not available in the AVI or the Big YUV file formats.\n');
    fprintf('\n');
    fprintf('  PSNR_VFD uses the results file output by the PSNR_SEARCH program\n');
    fprintf('  ("psnr_file") as a calibration starting point (i.e., Yshift, Xshift,\n');
    fprintf('  Tshift, Gain, Offset), but goes one step further by applying VFD\n');
    fprintf('  matching of the original video stream to the processed video stream\n');
    fprintf('  (i.e., the original video stream is modified so that it matches the\n');
    fprintf('  processed video stream frame-by-frame, or field-by-field for interlaced\n');
    fprintf('  video).  Output results from the PSNR_VFD program are stored in a file\n');
    fprintf('  with the same root name as "psnr_file" (i.e., name without the file\n');
    fprintf('  extension), appended with "_vfd.csv" (for Comma Separated Value).\n');
    fprintf('\n');
    fprintf('SYNTAX\n');
    fprintf('  psnr_vfd ''clip_dir'' ''test'' ''scan_type'' ''psnr_file'' options\n');
    fprintf('\n');
    fprintf('DESCRIPTION\n');
    fprintf('  This program uses the clip calibration information (i.e., Yshift,\n');
    fprintf('  Xshift, Tshift , Gain, Offset) from PSNR_SEARCH (in "psnr_file").  The\n');
    fprintf('  original clip is shifted by (Yshift, Xshift, Tshift) with respect to\n');
    fprintf('  the processed clip.  The Y-image of the processed clip is multiplied\n');
    fprintf('  by Gain and then the Offset is added.  For speed, all calibration is\n');
    fprintf('  held constant for the VFD estimation.  The VFD estimation algorithm is\n');
    fprintf('  applied to find the best matching original frame (or field) for each\n');
    fprintf('  processed frame (or field).  Then, the original is VFD-matched to the\n');
    fprintf('  processed.  A final gain and offset correction (called Gain_Adjust and\n');
    fprintf('  Offset_Adjust) is applied to the processed video (i.e.,\n');
    fprintf('  Gain_Adjust*y_proc + Offset_Adjust) before calculation of the final\n');
    fprintf('  VFD-corrected PSNR (PSNR_VFD).  Two perception-based VFD parameters are\n');
    fprintf('  also calculated (Par1_VFD and Par2_VFD) that attempt to capture the\n');
    fprintf('  perceptual distortions in the flow of scene motion (since these\n');
    fprintf('  distortions are removed from PSNR_VFD).  Par1_VFD is extracted from\n');
    fprintf('  only the VFD information while Par2_VFD uses both the VFD information\n');
    fprintf('  and the motion in the processed video clip.  For a description of these\n');
    fprintf('  two parameters, see the 2011 NTIA Technical Memorandum (TM) entitled\n');
    fprintf('  "Variable Frame Delay (VFD) Parameters for Video Quality Measurements."\n');
    fprintf('\n');
    fprintf('  The above procedure is repeated for each processed clip in the\n');
    fprintf('  "clip_dir" that belongs to the video test specified by "test".\n');
    fprintf('\n');
    fprintf('  A peak signal of 255 is used for calculation of PSNR.  Double precision\n');
    fprintf('  calculations are used everywhere.  A 64-bit operating system with at\n');
    fprintf('  least 4 GB of free memory is recommended since the entire double\n');
    fprintf('  precision versions of the original and processed sequences must be held\n');
    fprintf('  in memory.\n');
    fprintf('\n');
    fprintf('  Any or all of the following optional properties may be requested (the\n');
    fprintf('  first option is required for yuv files, not avi files).\n');
    fprintf('\n');
    fprintf('  ''yuv'' rows cols    Specifies the number of rows and cols for the Big\n');
    fprintf('                     YUV files (if using Big YUV files).  The default is\n');
    fprintf('                     to assume AVI files (*.avi).  Big YUV format is a\n');
    fprintf('                     binary format for storing ITU-R Recommendation\n');
    fprintf('                     BT.601 video sequences.  The format can be used for\n');
    fprintf('                     any image size.  In the Big YUV format, all the\n');
    fprintf('                     frames are stored sequentially in one big binary\n');
    fprintf('                     file. The sampling is 4:2:2 and image pixels are\n');
    fprintf('                     stored sequentially by video scan line as bytes in\n');
    fprintf('                     the following order: Cb1 Y1 Cr1 Y2 Cb3 Y3 Cr3 Y4…,\n');
    fprintf('                     where Y is the luminance component, Cb is the blue\n');
    fprintf('                     chrominance component, Cr is the red chrominance\n');
    fprintf('                     component, and the subscript is the pixel number.\n');
    fprintf('\n');
    fprintf('  ''sroi'' top left bottom right    Only use the specified spatial region\n');
    fprintf('                                  of interest (sroi) of the processed\n');
    fprintf('                                  video clip for all calculations. The\n');
    fprintf('                                  sroi is inclusive, where top/left start\n');
    fprintf('                                  at 1. By default, sroi is the entire\n');
    fprintf('                                  image reduced by the calibration shift\n');
    fprintf('                                  (Yshift, Xshift).  If the user inputs a\n');
    fprintf('                                  sroi, allowance must be made for the\n');
    fprintf('                                  maximum spatial shift encountered in\n');
    fprintf('                                  the PSNR_SEARCH results ("psnr_file").\n');
    fprintf('                                  For interlaced video, top must be odd\n');
    fprintf('                                  while bottom must be even.\n');
    fprintf('\n');
    fprintf('  ''troi'' fstart fstop    Only calculate PSNR_VFD for the specified\n');
    fprintf('                         temporal region of interest (troi) of the\n');
    fprintf('                         processed video clip, where fstart and fstop are\n');
    fprintf('                         included and given in frames.  By default, the\n');
    fprintf('                         troi is the entire file reduced by the temporal\n');
    fprintf('                         calibration shift (Tshift).  If the user inputs\n');
    fprintf('                         an fstart and fstop, allowance must be made for\n');
    fprintf('                         the maximum temporal shift encountered in\n');
    fprintf('                         "psnr_file".\n');
    fprintf('\n');
    fprintf('  ''t_uncert'' t    Specifies the temporal uncertainty (plus or minus t\n');
    fprintf('                  frames) over which to perform the VFD search.  The\n');
    fprintf('                  processed remains fixed and the original is shifted.\n');
    fprintf('                  The center (zero shift) point for the temporal search\n');
    fprintf('                  assumes the temporal alignment given by Tshift in the\n');
    fprintf('                  "psnr_file" results from PSNR_SEARCH. By default,\n');
    fprintf('                  temporal uncertainty is set to 30 frames.  It can have\n');
    fprintf('                  a minimum value of 1 frame.  When the original cannot\n');
    fprintf('                  be shifted by the temporal uncertainty (e.g., perhaps\n');
    fprintf('                  near the ends of the sequence), the original will be\n');
    fprintf('                  shifted up to the maximum extent possible.  To\n');
    fprintf('                  accomodate a full temporal search at the beginning and\n');
    fprintf('                  end of the sequence, increase fstart and decrease fstop\n');
    fprintf('                  in the processed troi accordingly.\n');
    fprintf('\n');
    fprintf('  ''reframe''  Allow for the possibility that the processed video clip has\n');
    fprintf('             been reframing.  This option is only valid for a scan_type\n');
    fprintf('             of ''interlaced_uff'' or ''interlaced_lff''.  Reframing can vary\n');
    fprintf('             throughout the processed clip, although this should be rare.\n');
    fprintf('             This option will increase the runtime substantially since\n');
    fprintf('             extra spatial shifts must be examined, but it should be used\n');
    fprintf('             if there is any possibility of time varying reframing\n');
    fprintf('             existing in the processed video clip.  See Section 3.1.2 of\n');
    fprintf('             NTIA Report TR-02-392 for a definition of reframing.  For\n');
    fprintf('             constant reframing that was properly detected by the\n');
    fprintf('             PSNR_SEARCH program, Tshift will have a half frame (0.5)\n');
    fprintf('             Tshift.  This condition will be detected/corrected by\n');
    fprintf('             PSNR_VFD (the original video is reframed accordingly).\n');
    fprintf('\n');
    fprintf('  ''causal''   Impose causality constraint so that later frames (fields) in\n');
    fprintf('             the processed clip cannot align to original frames (fields)\n');
    fprintf('             that are earlier in time than found for the proceeding\n');
    fprintf('             processed frames (fields).  For interlaced video, a\n');
    fprintf('             one-field jump back in time is allowed since this is\n');
    fprintf('             indicative of a frozen frame.  By default, causality is\n');
    fprintf('             turned off (yes, codecs can output non-causal sequences).\n');
    fprintf('             But specifying the causal option is usually recommended.\n');
    fprintf('\n');
    fprintf('  ''verbose''   Display intermediate progress during processing.\n');
    fprintf('\n');
    fprintf('EXAMPLES\n');
    fprintf('  These three examples illustrate the Big YUV format for QCIF, CIF, and VGA\n');
    fprintf('\n');
    fprintf('  psnr_vfd ''d:\\q01\\'' ''q01'' ''progressive'' ''q01_psnr.csv'' ''yuv'' 144 176 ''causal''\n');
    fprintf('  psnr_vfd ''d:\\c01\\'' ''c01'' ''progressive'' ''c01_psnr.csv'' ''yuv'' 288 352 ''causal''\n');
    fprintf('  psnr_vfd ''d:\\v01\\'' ''v01'' ''progressive'' ''v01_psnr.csv'' ''yuv'' 480 640 ''causal''\n');
    fprintf('\n');
    fprintf('  This example illustrates uncompressed UYVY AVI format for 525-line video.\n');
    fprintf('\n');
    fprintf('  psnr_vfd ''d:\\rr525\\'' ''rr525'' ''interlaced_lff'' ''rr525_psnr.csv'' ''causal'' ''reframe''\n');
    fprintf('\n');
    return;
end

% Strip off the extra single quotes '' on the required inputs.  This is
% required for command line arguments in standalone execuatables.
clip_dir = eval(clip_dir); 
test = eval(test);
scan_type = eval(scan_type);
psnr_file = eval(psnr_file);

% Generate the name of the file to store PSNR_VFD results
dot = strfind(psnr_file,'.');
vfd_file = strcat(psnr_file(1:dot(length(dot))-1), '_vfd.csv');

% Validate the scan_type
if (~strcmpi(scan_type,'progressive') && ~strcmpi(scan_type,'interlaced_lff') && ~strcmpi(scan_type,'interlaced_uff'))
    error('Invalid scan_type');
end

% Define the peak signal level
peak = 255.0;

% Define the sub-sampling factor on the pixels for the final gain and
% offset adjusting fit, which is performed right before calculation of
% PSNR_VFD.
fraction_sampled = 0.1;

% Add extra \ in clip_dir in case user did not
clip_dir = strcat(clip_dir,'\');

% Validate input arguments and set their defaults
file_type = 'avi';  % default file type, uncompressed UYVY AVI
is_yuv = 0;
is_whole_image = 1;
is_whole_time = 1;
t_uncert = 30;  % Default plus or minus temporal search, in frames
reframe = 0;
causal = 0;
verbose = 0;

cnt=1;
while cnt <= length(varargin),
    if strcmpi(eval(char(varargin(cnt))),'yuv') == 1
        rows = str2double(varargin{cnt+1});
        cols = str2double(varargin{cnt+2});
        is_yuv = 1;
        file_type = 'yuv';
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
    elseif strcmpi(eval(char(varargin(cnt))), 't_uncert') ==1
        t_uncert = str2double(varargin{cnt+1});
        cnt = cnt + 2;
    elseif strcmpi(eval(char(varargin(cnt))),'reframe') == 1
        reframe = 1;
        cnt = cnt + 1;
        % Make sure video is not progressive for this option
        if (strcmpi(scan_type,'progressive'))
            error('Reframe option not allowed for progressive video');
        end
    elseif strcmpi(eval(char(varargin(cnt))),'causal') == 1
        causal = 1;
        cnt = cnt + 1;
    elseif strcmpi(eval(char(varargin(cnt))),'verbose') == 1
        verbose = 1;
        cnt = cnt +1;
    else
        error('Property value passed into psnr_vfd not recognized');
    end
end

% If not progressive and user inputs an SROI, they must have an odd top and
% an even bottom.  Otherwise the field ordering will reverse.
if (~strcmpi(scan_type,'progressive') && ~is_whole_image && (~mod(top,2) || mod(bottom,2)))
    error('SROI top must be odd and bottom must be even for interlaced video.');
end

%  Get a directory listing
files = dir(clip_dir);  % first two files are '.' and '..'
num_files = size(files,1);

% Find the HRCs and their scenes for the specified video test
hrc_list = {};
scene_list = {};
for i=3:num_files
    this_file = files(i).name;
    und = strfind(this_file,'_'); % find underscores and period
    dot = strfind(this_file,'.');  % Will only use the last dot
    if(size(und,2)==2) % possible standard naming convention file found
        this_test = this_file(1:und(1)-1);  % pick off the test name
        if(~isempty(strmatch(test,this_test,'exact')) && ...
                ~isempty(strmatch(file_type,this_file(dot(length(dot))+1:length(this_file)),'exact')))  % test clip found
            this_scene = this_file(und(1)+1:und(2)-1);
            this_hrc = this_file(und(2)+1:dot(length(dot))-1);
            % See if this HRC already exists and find its list location
            loc = strmatch(this_hrc,hrc_list,'exact');
            if(loc)  % HRC already present, add to scene list for that HRC
                if(size(strmatch(this_scene,scene_list{loc},'exact'),1)==0)
                    scene_list{loc} = [scene_list{loc} this_scene];
                end
            else  % new HRC found
                hrc_list = [hrc_list;{this_hrc}];
                this_loc = size(hrc_list,1);
                scene_list(this_loc) = {{this_scene}};
            end
        end
    end
end

scene_list = scene_list';
num_hrcs = size(hrc_list,1);
if (num_hrcs == 0)
    error('No files with standard naming convention found.');
end

% Results struct to store psnr_vfd results
results_vfd = struct('test', {}, 'scene', {}, 'hrc', {}, 'gain_adjust', {}, 'offset_adjust', {}, ...
    'psnr_vfd', {}, 'par1_vfd', {}, 'par2_vfd', {}, 'orig_indices', {}, 'proc_indices', {});

% Read in the psnr results output by the PSNR_SEARCH program.  Test clips
% that do not have PSNR_SEARCH results will be skipped.
psnr_import = importdata(psnr_file);
%  Check to make sure that the imported structure has the correct
%  characteristics for the textdata and data arrays
if (size(psnr_import.textdata,2)~=9 || size(psnr_import.data,2)~=6 || size(psnr_import.textdata,1)-1 ~= size(psnr_import.data,1))
    error('Invalid psnr_file.');
end
nclips = size(psnr_import.data,1);  % The number of clips in the file
% Load the structure to hold the PSNR_SEARCH results
results_psnr = struct('test', {}, 'scene', {}, 'hrc', {}, 'yshift', {}, ...
    'xshift', {}, 'tshift', {}, 'gain', {}, 'offset', {}, 'psnr', {});
for i = 1:nclips
    results_psnr(i).test = psnr_import.textdata{i+1,1};
    results_psnr(i).scene = psnr_import.textdata{i+1,2};
    results_psnr(i).hrc = psnr_import.textdata{i+1,3};
    results_psnr(i).yshift = psnr_import.data(i,1);
    results_psnr(i).xshift = psnr_import.data(i,2);
    results_psnr(i).tshift = psnr_import.data(i,3);
    results_psnr(i).gain = psnr_import.data(i,4);
    results_psnr(i).offset = psnr_import.data(i,5);
    results_psnr(i).psnr = psnr_import.data(i,6);
end

% Process one HRC at a time to compute average PSNR_VFD for that HRC
index = 1;  % index used to store psnr_vfd results
fid_vfd = fopen(vfd_file,'a');  % open vfd_file for appending and write out the header
if (strcmpi(scan_type,'progressive'))
    fprintf(fid_vfd,'Test,Scene,HRC,Gain_Adjust,Offset_Adjust,PSNR_VFD,Par1_VFD,Par2_VFD,(Proc Orig) Matching Frame Indices\n');
else % interlaced
    fprintf(fid_vfd,'Test,Scene,HRC,Gain_Adjust,Offset_Adjust,PSNR_VFD,Par1_VFD,Par2_VFD,(Proc Orig) Matching Field Indices\n');
end
fclose(fid_vfd);
for i = 1:num_hrcs
    
    psnr_ave = 0;  % psnr_vfd average summer for this HRC
    par1_ave = 0;  % par1_vfd average summer for this HRC
    par2_ave = 0;  % par2_vfd average summer for this HRC
    this_hrc = hrc_list{i};
    if(strcmpi('original',this_hrc)) % Don't process original
        continue;
    end
    num_scenes = size(scene_list{i},2);  % Number of scenes in this HRC
    
    for j = 1:num_scenes
        
        this_scene = scene_list{i}{j};
        
        %  Find this clip's calibration information in results_psnr
        this_clip  = find((strcmpi({results_psnr.test},test) & strcmpi({results_psnr.scene},this_scene) & strcmpi({results_psnr.hrc},this_hrc)));
        if (isempty(this_clip) && verbose)
            fprintf('Skipping Clip %s_%s_%s, Calibration information not found in psnr_file.\n', test, this_scene, this_hrc);
            continue;
        end
        
        %  Assign the scene information to the results_vfd structure
        results_vfd(index).test = test;
        results_vfd(index).scene = this_scene;
        results_vfd(index).hrc = this_hrc;
        
        %  Pick off the calibration information for this clip
        this_yshift = results_psnr(this_clip).yshift;
        this_xshift = results_psnr(this_clip).xshift;
        this_tshift = results_psnr(this_clip).tshift;
        this_gain = results_psnr(this_clip).gain;
        this_offset = results_psnr(this_clip).offset;
        this_psnr = results_psnr(this_clip).psnr;  % This will be used as a check against psnr_vfd
        
        % Read original and processed video files
        if (~is_yuv)  % YUV file parameters not specified, AVI assummed
            % Re-generate the original and processed avi file names
            orig = strcat(clip_dir, test,'_', this_scene, '_', 'original', '.avi');
            proc = strcat(clip_dir, test,'_', this_scene, '_', this_hrc, '.avi');
            [avi_info] = read_avi('Info',orig);
            [avi_info_proc] = read_avi('Info',proc);
            rows = avi_info_proc.Height;
            cols = avi_info_proc.Width;
            % Check to make sure processed and original sizes match
            rows_orig = avi_info.Height;
            cols_orig = avi_info.Width;
            if (rows ~= rows_orig || cols ~= cols_orig)
                error('Original and processed image sizes do not match.')
            end
            
            % Set/Validate the SROI of the processed video
            if (is_whole_image) % make SROI the whole image less the calibration shift
                if (this_xshift <= 0)  % Original is shifted left or 0 wrt processed
                    left = 1-this_xshift;
                    right = cols;
                else  % Original is shifted right wrt processed
                    left = 1;
                    right = cols-this_xshift;
                end
                if (this_yshift <= 0)  % Original is shifted up or 0 wrt processed
                    top = 1-this_yshift;
                    if(~strcmpi(scan_type,'progressive') && ~mod(top,2))  % Must start on odd line for interlaced video
                        top = top + 1;
                    end
                    bottom = rows;
                else  % Original is shifted down wrt processed
                    top = 1;
                    bottom = rows-this_yshift;
                    if(~strcmpi(scan_type,'progressive') && mod(bottom,2))  % Must end on even line for interlaced video
                        bottom = bottom - 1;
                    end
                end
            end
            if (top<1 || left<1 || bottom>rows || right>cols)
                fprintf('Skipping Clip %s_%s_%s, invalid processed SROI, top=%f, left=%f, bottom=%f, right = %f.\n', ...
                    test, this_scene, this_hrc, top, left, bottom, right);
                continue;
            end
            
            %  Set the matching original SROI and validate
            left_orig = left + this_xshift;
            right_orig = right + this_xshift;
            top_orig = top + this_yshift;
            bottom_orig = bottom + this_yshift;
            % Odd y_shift, correct to preserve field ordering for interlaced, new top and bottom create two extra lines that 
            % will be eliminated in the reframe.
            if (mod(this_yshift,2) && ~strcmpi(scan_type,'progressive'))  
                top_orig = top_orig - 1;
                bottom_orig = bottom_orig + 1;
            end
            if (top_orig<1 || left_orig<1 || bottom_orig>rows || right_orig>cols)  % Original SROI wrt processed SROI
                fprintf('Skipping Clip %s_%s_%s, original xshift=%f and yshift=%f\n', ...
                    test, this_scene, this_hrc, this_xshift, this_yshift);
                fprintf('produces invalid SROI top_orig=%f, left_orig=%f, bottom_orig=%f, right_orig=%f.\n', ...
                    top_orig, left_orig, bottom_orig, right_orig);
                continue;
            end
            
            tframes_orig = avi_info.NumFrames;  % total frames in orig file
            tframes_proc = avi_info_proc.NumFrames;
            % Validate that orig and proc have the same number of frames
            if (tframes_orig ~= tframes_proc)
                fprintf('\n%s_%s_%s: orig & proc files have different number of frames; longer file will be truncated.\n', ...
                    test, this_scene, this_hrc);
                tframes = min(tframes_orig,tframes_proc); % the actual number of frames that will be used
            else
                tframes = tframes_proc;
            end
            
            % Set/Validate the time segment of the processed video
            if (is_whole_time) % use whole time segment less the calibration shift
                if (this_tshift <= 0)  % original is shifted left or 0 wrt processed
                    fstart = ceil(1-this_tshift);
                    fstop = tframes;
                else  % original is shifted right wrt processed
                    fstart = 1;
                    fstop = floor(tframes-this_tshift);
                end
            end
            if (fstart<1 || fstop>tframes)
                fprintf('Skipping Clip %s_%s_%s, invalid processed TROI, fstart=%f, fstop=%f.\n', ...
                    test, this_scene, this_hrc, fstart, fstop);
                continue;
            end
            
            % Set the matching original fstart and fstop and validate.  The
            % original will contain an extra frame when reframing is
            % required.
            fstart_orig = floor(fstart+this_tshift);
            fstop_orig = ceil(fstop+this_tshift);
            if (fstart_orig<1 || fstop_orig>tframes)  % Original TROI wrt Processed TROI
                fprintf('Skipping Clip %s_%s_%s, original tshift=%f produces\n', test, this_scene, this_hrc, this_tshift);
                fprintf('invalid original TROI: fstart_orig=%f, fstop_orig=%f.\n', fstart_orig, fstop_orig);
                continue;
            end
            
            %  Calculate to see how many extra original frames there are at
            %  the beginning and end of the sequence for VFD calculations.
            %  We would like at least t_uncert extra frames.
            beg_extra = min(fstart_orig-1, t_uncert);
            end_extra = min(tframes_orig-fstop_orig, t_uncert);
            first_align = 1 + beg_extra;  % This original frame best aligns to the first frame in the processed TROI
            
            % Read in video and clear color planes to free up memory
            [y_orig,cb,cr] = read_avi('YCbCr',orig,'frames',fstart_orig-beg_extra,...
                fstop_orig+end_extra, 'sroi',top_orig,left_orig,...
                bottom_orig,right_orig);
            clear cb cr;
            [y_proc,cb,cr] = read_avi('YCbCr',proc,'frames',fstart,fstop,...
                'sroi',top,left,bottom,right);
            clear cb cr;
            
        else  % YUV file
            % Re-generate the original and processed YUV file name
            orig = strcat(clip_dir, test,'_', this_scene, '_', 'original', '.yuv');
            proc = strcat(clip_dir, test,'_', this_scene, '_', this_hrc, '.yuv');
            
            % Set/Validate the SROI of the processed video
            if (is_whole_image) % make SROI the whole image less the calibration shift
                if (this_xshift <= 0)  % Original is shifted left or 0 wrt processed
                    left = 1-this_xshift;
                    right = cols;
                else  % Original is shifted right wrt processed
                    left = 1;
                    right = cols-this_xshift;
                end
                if (this_yshift <= 0)  % Original is shifted up or 0 wrt processed
                    top = 1-this_yshift;
                    if(~strcmpi(scan_type,'progressive') && ~mod(top,2))  % Must start on odd line for interlaced video
                        top = top + 1;
                    end
                    bottom = rows;
                else  % Original is shifted down wrt processed
                    top = 1;
                    bottom = rows-this_yshift;
                    if(~strcmpi(scan_type,'progressive') && mod(bottom,2))  % Must end on even line for interlaced video
                        bottom = bottom - 1;
                    end
                end
            end
            if (top<1 || left<1 || bottom>rows || right>cols)
                fprintf('Skipping Clip %s_%s_%s, invalid processed SROI, top=%f, left=%f, bottom=%f, right = %f.\n', ...
                    test, this_scene, this_hrc, top, left, bottom, right);
                continue;
            end
            
            %  Set the matching original SROI and validate
            left_orig = left + this_xshift;
            right_orig = right + this_xshift;
            top_orig = top + this_yshift;
            bottom_orig = bottom + this_yshift;
            % Odd y_shift, correct to preserve field ordering for interlaced, new top and bottom create two extra lines that 
            % will be eliminated in the reframe.
            if (mod(this_yshift,2) && ~strcmpi(scan_type,'progressive'))  
                top_orig = top_orig - 1;
                bottom_orig = bottom_orig + 1;
            end
            if (top_orig<1 || left_orig<1 || bottom_orig>rows || right_orig>cols)  % Original SROI wrt processed SROI
                fprintf('Skipping Clip %s_%s_%s, original xshift=%f and yshift=%f\n', ...
                    test, this_scene, this_hrc, this_xshift, this_yshift);
                fprintf('produces invalid SROI, top_orig=%f, left_orig=%f, bottom_orig=%f, right_orig=%f.\n', ...
                    top_orig, left_orig, bottom_orig, right_orig);
                continue;
            end
            
            % Find the total frames of the input original file
            [fid, message] = fopen(orig, 'r');
            if fid == -1
                fprintf(message);
                error('Cannot open this clip''s bigyuv file, %s', orig);
            end
            % Find last frame.
            fseek(fid,0, 'eof');
            tframes_orig = ftell(fid) / (2 * rows * cols);
            fclose(fid);
            % Find the total frames of the processed file
            [fid, message] = fopen(proc, 'r');
            if fid == -1
                fprintf(message);
                error('Cannot open this clip''s bigyuv file, %s', proc);
            end
            % Find last frame.
            fseek(fid,0, 'eof');
            tframes_proc = ftell(fid) / (2 * rows * cols);
            fclose(fid);
            % Validate that orig and proc have the same number of frames
            if (tframes_orig ~= tframes_proc)
                fprintf('\n%s_%s_%s: orig & proc files have different number of frames; longer file will be truncated.\n', ...
                    test, this_scene, this_hrc);
                tframes = min(tframes_orig,tframes_proc);
            else
                tframes = tframes_proc;
            end
            
            % Set/Validate the time segment of the processed video
            if (is_whole_time) % use whole time segment less the calibration shift
                if (this_tshift <= 0)  % original is shifted left or 0 wrt processed
                    fstart = ceil(1-this_tshift);
                    fstop = tframes;
                else  % original is shifted right wrt processed
                    fstart = 1;
                    fstop = floor(tframes-this_tshift);
                end
            end
            if (fstart<1 || fstop>tframes)
                fprintf('Skipping Clip %s_%s_%s, invalid processed TROI, fstart=%f, fstop=%f.\n', ...
                    test, this_scene, this_hrc, fstart, fstop);
                continue;
            end
            
            % Set the matching original fstart and fstop and validate.  The original will contain an extra frame 
            % when reframing is required.
            fstart_orig = floor(fstart+this_tshift);
            fstop_orig = ceil(fstop+this_tshift);
            if (fstart_orig<1 || fstop_orig>tframes)  % Original TROI wrt Processed TROI
                fprintf('Skipping Clip %s_%s_%s, original tshift=%f produces\n', test, this_scene, this_hrc, this_tshift);
                fprintf('invalid original TROI: fstart_orig=%f, fstop_orig=%f.\n', fstart_orig, fstop_orig);
                continue;
            end
            
            %  Calculate to see how many extra original frames there are at
            %  the beginning and end of the sequence for VFD calculations.
            %  We would like at least t_uncert extra frames.
            beg_extra = min(fstart_orig-1, t_uncert);
            end_extra = min(tframes_orig-fstop_orig, t_uncert);
            first_align = 1 + beg_extra;  % This original frame best aligns to the first frame in the processed TROI
            
            % Read in video and clear color planes to free up memory
            [y_orig,cb,cr] = read_bigyuv(orig,'frames',fstart_orig-beg_extra,...
                fstop_orig+end_extra,'size',rows,cols,'sroi',top_orig,...
                left_orig,bottom_orig,right_orig);
            clear cb cr;
            [y_proc,cb,cr] = read_bigyuv(proc,'frames',fstart,fstop,...
                'size',rows,cols,'sroi',top,left,bottom,right);
            clear cb cr;
            
        end
        
        %  Reframe the original if required
        if (rem(this_tshift,1))  % Non-integer tshift
            y_orig = reframe_video(y_orig, scan_type);
        end
        
        % Convert everything to double precision before any calculations
        % are performed.
        y_orig = double(y_orig);  
        y_proc = double(y_proc);
        
        %  Correct the processed for gain and offset as read in from the
        %  psnr_file.
        y_proc = this_gain*y_proc + this_offset;
        
        %  Generate the call to the VFD estimation function
        func_call = 'est_var_frame_delays(y_proc,y_orig,''normalize'',';  % use normalize option as it seems to work best
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
            first_align = 2*first_align-1;  % convert to fields for est_var_frame_delays
        end
        if (strcmpi(scan_type,'interlaced_uff'))
            func_call = strcat(func_call,'''interlaced'',2,');
            first_align = 2*first_align-1;  % convert to fields for est_var_frame_delays
        end
        func_call = strcat(func_call,'''first_align'',',num2str(first_align),',','''t_uncert'',',num2str(t_uncert),')');
        
        % Call the est_var_frame_delays function to get VFD results
        [results results_rmse results_fuzzy results_fuzzy_mse] = eval(func_call);
        
        %  If the VFD alignment algorithm failed, then results==0.
        %  In that case, use the time alignment given by the psnr_file and
        %  generate psuedo results where the alignment just increases by
        %  one frame (or field) at a time.
        [nrows, ncols, nframes] = size(y_proc);
        if (results == 0)
            vfd_failed = 1;  % Set a logical variable to record that the VFD algorithm failed.
            if (~strcmpi(scan_type,'progressive'))  % interlaced
                npts = nframes*2;
            else  % progressive
                npts = nframes;
            end
            results = first_align:first_align+npts-1;
            fprintf('WARNING: VFD algorithm failed for clip %s_%s_%s, using psnr_file time alignment.\n', ...
                test, this_scene, this_hrc);
        else
            vfd_failed = 0;
            npts = length(results);
        end
        
        % This code translates the VFD results to use the orig and proc FILE indexing
        if (strcmpi(scan_type,'progressive'))
            proc_indices = (fstart-1) + (1:length(results));
            orig_indices = (fstart_orig-beg_extra-1) + results;
        else % interlaced
            proc_indices = 2*(fstart-1) + (1:length(results));
            orig_indices = 2*(fstart_orig-beg_extra-1) + results;
            % Add one to orig_indices if the original was reframed
            if (rem(this_tshift,1))  % Non-integer tshift
                orig_indices = orig_indices+1;
            end
        end
        results_vfd(index).orig_indices = orig_indices;
        results_vfd(index).proc_indices = proc_indices;
        
        % Apply the VFD correction to the original clip to make it look
        % like the processed clip.
        y_orig = vfd_match(y_orig, scan_type, results);
        
        % Reshape for gain/offset fit and PSNR calculation
        y_proc = reshape(y_proc,nrows*ncols*nframes,1);
        y_orig = reshape(y_orig,nrows*ncols*nframes,1);
        
        % Perform the final gain and offset fit using randomly sub-sampled pixels
        rand_nums = round(randperm((nrows*ncols*nframes))); %Randomizes numbers from 1 to nrows*ncols*nframes
        this_fit = polyfit(y_proc(rand_nums(1:round(nrows*ncols*nframes*fraction_sampled))),...
                                y_orig(rand_nums(1:round(nrows*ncols*nframes*fraction_sampled))),1);
        clear rand_nums;
        results_vfd(index).gain_adjust = this_fit(1);
        results_vfd(index).offset_adjust = this_fit(2);
        
        %  Calculate the final PSNR_VFD
        this_psnr_vfd = 10*(log10(peak*peak)-log10(sum(((this_fit(1)*y_proc+this_fit(2))-y_orig).^2)/(nrows*ncols*nframes)));
        results_vfd(index).psnr_vfd = this_psnr_vfd;
        clear y_orig;  % Done with y_orig
        
        %  Print out a warning if the psnr_vfd is significantly lower than
        %  the psnr, which indicates that the VFD alignment might be
        %  suspect.  Here we are using a threshold of 1.5 db.
        if (this_psnr-this_psnr_vfd > 1.5)
            fprintf('WARNING: VFD alignment results may be unreliable for clip %s_%s_%s, psnr_vfd=%5.4f, psnr=%5.4f.\n', ...
                test, this_scene, this_hrc, this_psnr_vfd, this_psnr);
        end
        
        %  Reshape y_proc to 3D and calculate TI_RMS
        y_proc = reshape(y_proc,nrows,ncols,nframes);
        y_proc = cat(3,y_proc(:,:,1),y_proc);  % Add extra frame at the beginning for TI calculation
        y_proc = diff(y_proc,1,3);  % first order difference along 3rd dimension
        y_proc = reshape(y_proc,nrows*ncols,nframes);
        y_proc = y_proc.^2;
        y_proc = sqrt(sum(y_proc)./(nrows*ncols));
        if (~strcmpi(scan_type,'progressive'))  % Replicate every other sample for interlaced video
            y_proc = reshape(repmat(y_proc,2,1), 1, 2*nframes);
        end
        
        %  Calculate the diff of the VFD information, which forms the basis
        %  for both par1_vfd and par2_vfd.  The VFD information is
        %  converted to a vector that gives Abnormal Frame Jumps (AFJs).
        % Subtracting 1 and maxing with zero produces a parameter that (1)
        % does not penalize for normal field/frame delivery (where the VFD
        % field/frame indices increase by one from one field/frame to the
        % next), (2) does not penalize for frame/field repeats (e.g., where
        % the VFD frame indices stay fixed from one frame to the next), and
        % (3) does not penalize for interlaced frame repeats (where the VFD
        % field indices jump back one in time from one field to the next).
        % A non-impairment value of 0 is used for the first field/frame
        % (which must be padded since it's diff is not available).  
        if (vfd_failed)
            
            % This equation assumes that both the early and late time sides
            % used for the frame jump estimates are absolutely correct
            % (i.e., no fuzzy alignments).  The fuzzy alignment information
            % is not available here as the VFD algorithm failed.  Abnormal
            % Frame Jumps (AFJ) is then calculated as:
            afj = max([0 abs(diff(results))-1], 0);
            
        else
            
            % This code assumes fuzzy uncertainty on both the early and
            % late time sides when calculating frame jumps (the uncertainty
            % is given by results_fuzzy array).  Here, frame jumps are only
            % penalized when they are absolutely certain to be correct
            % (i.e., no fuzzy overlapping). 
            fuzzy_max_early = min(max(results_fuzzy(:,1:npts-1)),results(2:npts));
            fuzzy_min_late = max(min(results_fuzzy(:,2:npts)),fuzzy_max_early);
            afj = max([0 fuzzy_min_late-fuzzy_max_early-1], 0);
            
        end
        
        %  Calculate the pure VFD parameter (par1_vfd) and the TI weighted
        %  variant VFD parameter (par2_vfd).
        par1_vfd = log10(sqrt(mean(afj.^2))+1);
        results_vfd(index).par1_vfd = par1_vfd;
        par2_vfd = log10(sqrt(mean((afj.*log10(1+y_proc)).^2))+1);
        results_vfd(index).par2_vfd = par2_vfd;
        
        %  Output the clip information, current time, and the psnr_vfd
        t = clock;
        if (verbose)
            fprintf('Clip %s_%s_%s at %d:%d, psnr_vfd = %5.4f, par1_vfd = %5.4f, par2_vfd = %5.4f\n', ...
                test, this_scene, this_hrc, t(4), t(5), this_psnr_vfd, par1_vfd, par2_vfd);
        end
        
        % Write out the psnr_vfd results for this clip into the vfd_file
        fid_vfd = fopen(vfd_file,'a');
        fprintf(fid_vfd,'%s, %s, %s, %5.4f, %5.4f, %5.4f, %5.4f, %5.4f,', results_vfd(index).test, results_vfd(index).scene, ...
            results_vfd(index).hrc, results_vfd(index).gain_adjust, results_vfd(index).offset_adjust, ...
            results_vfd(index).psnr_vfd, results_vfd(index).par1_vfd, results_vfd(index).par2_vfd);
        for k = 1:npts-1
            fprintf(fid_vfd,'%d %d,',proc_indices(k), orig_indices(k));
        end
        fprintf(fid_vfd,'%d %d\n',proc_indices(k+1), orig_indices(k+1));
        fclose(fid_vfd);
        
        %  Add to the HRC summer
        psnr_ave = psnr_ave + this_psnr_vfd;
        par1_ave = par1_ave + par1_vfd;
        par2_ave = par2_ave + par2_vfd;
        
        %  Increment the results_vfd counter
        index = index+1;
        
    end
    
    % Compute average psnr_vfd for this HRC
        psnr_ave = psnr_ave/(num_scenes);
        par1_ave = par1_ave/(num_scenes);
        par2_ave = par2_ave/(num_scenes);
        if(verbose)
            fprintf('\nHRC = %s, psnr_vfd_ave = %5.4f, par1_vfd_ave = %5.4f, par2_vfd_ave = %5.4f\n\n', ...
                this_hrc, psnr_ave, par1_ave, par2_ave);
        end
        
        pause(1);
        close all;  % closes open figures
        fclose('all');  % closes open files
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [yout] = reframe_video(yin, scan_type)
%  function [yout] = reframe_video(yin, scan_type)
%  This function reframes a 3D input video array yin, with a scan_type of
%  either 'interlaced_uff' or 'interlaced_lff', and produces a 3D output
%  video array yout.  yout will have one less frame than yin and have its
%  number of rows reduced by two lines.  The number of rows in yin must be
%  even and yin must contain at least 2 video frames.
%
[nrows, ncols, nframes] = size(yin);
if (mod(nrows,2) || nframes<2)
    error ('reframe_video function requires an even number of rows and at least 2 video frames');
end

% Split_into_fields
yin = reshape(yin, 2, nrows/2, ncols, nframes);
if (strcmpi(scan_type,'interlaced_lff'))
    late_field = squeeze(yin(1,2:nrows/2,:,1:nframes-1));
    early_field = squeeze(yin(2,1:nrows/2-1,:,2:nframes));
elseif (strcmpi(scan_type,'interlaced_uff'))
    early_field = squeeze(yin(1,2:nrows/2,:,2:nframes));
    late_field = squeeze(yin(2,1:nrows/2-1,:,1:nframes-1));
else
    error('Unsupported scan_type in function reframe_video');
end
clear yin;

% Reframe video and return: 
% For interlaced_lff the lower_field is now the late_field and the
% upper_field is now the early_field.  
% For interlaced_uff the lower_field is now the early_field and the
% upper_field is now the late_field. 
if (strcmpi(scan_type,'interlaced_lff'))  %  Reframe code for interlaced_lff
    yout(1,:,:,:) = early_field;
    yout(2,:,:,:) = late_field;
else  %  Reframe code for interlaced_uff
    yout(1,:,:,:) = late_field;
    yout(2,:,:,:) = early_field;
end
clear late_field early_field;
yout = reshape(yout, nrows-2, ncols, nframes-1);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [yout] = vfd_match(yin, scan_type, results)
%  function [yout] = vfd_match(yin, scan_type, results)
%  This function converts 3D input video array 'yin' into another 3D video
%  array 'yout' with frames and/or fields ordered according to the Variable
%  Frame Delay (VFD) 'results'.  The scan_type must be either
%  'progressive', 'interlaced_lff', or 'interlaced_uff'.
%

[nrows, ncols, nframes_orig] = size(yin);

if (strcmpi(scan_type,'progressive'))
    nframes = length(results);
    is_interlaced = 0;
elseif (strcmpi(scan_type, 'interlaced_lff'))  
    nframes = length(results)/2;  % These are field results, so must half to get nframes
    is_interlaced = 1;
    field_first = 1;  % Use same field_first definition as est_var_frame_delays
elseif (strcmpi(scan_type, 'interlaced_uff'))
    nframes = length(results)/2;
    is_interlaced = 1;
    field_first = 2;
else
    error('Unsupported scan_type in function vfd_match.');
end

% yout will be the VFD-corrected original and will have the same number of
% frames as the processed clip.  Use the same 'single' or 'double' rule as
% read_tslice for the image precision of yout.
if (nrows > 650)
    yout = zeros(nrows,ncols,nframes,'single');
else
    yout = zeros(nrows,ncols,nframes,'double');
end

if(is_interlaced)
    
    for j = 1:nframes
        % Get matching original field for the early processed field
        orig_frame_num = ceil(results(2*j-1)/2);  % The frame number that contains the original field
        early_field = mod(results(2*j-1),2);  % =1 if early field, =0 if late field
        [yo1 yo2] = split_into_fields(squeeze(yin(:,:,orig_frame_num)));
        if (early_field)
            switch field_first
                case(1)
                    this_orig1 = yo1;
                case(2)
                    this_orig1 = yo2;
            end
        else % late field
            switch field_first
                case(1)
                    this_orig1 = yo2;
                case(2)
                    this_orig1 = yo1;
            end
        end
        % Get matching original field for the late processed field
        orig_frame_num = ceil(results(2*j)/2);  % The frame number that contains the original field
        early_field = mod(results(2*j),2);  % =1 if early field, =0 if late field
        [yo1 yo2] = split_into_fields(squeeze(yin(:,:,orig_frame_num)));
        if (early_field)
            switch field_first
                case(1)
                    this_orig2 = yo1;
                case(2)
                    this_orig2 = yo2;
            end
        else % late field
            switch field_first
                case(1)
                    this_orig2 = yo2;
                case(2)
                    this_orig2 = yo1;
            end
        end
        % Joint the two original fields into a frame
        switch field_first
            case(1)
                this_orig = join_into_frames(this_orig1,this_orig2);
            case(2)
                this_orig = join_into_frames(this_orig2,this_orig1);
        end
        yout(:,:,j) = this_orig;
    end
    clear this_orig this_orig1 this_orig2 yo1 yo2;
    
else  % progressive
    
    for j = 1:nframes
        yout(:,:,j) = yin(:,:,results(j));
    end
    
end

return

