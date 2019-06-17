function rr_original (file_list, video_standard, model, results_log, flag)
% wrapper for RRNR-TV test: interface required for compilation
% call with no arguments for help information.

if nargin == 0,
    fprintf('RR_ORIGINAL Version 1.2\n');
    fprintf('  Take original test sequences in uncompressed big-YUV file format.\n');
    fprintf('  Calculate intermediate data file for NTIA low-bandwidth model \n');
    fprintf('  (or NTIA fast low-bandwidth model) and calibration.\n');
    fprintf('SYNTAX\n');
    fprintf('  rr_original file_list video_standard model results_log\n');
    fprintf('  rr_original file_list video_standard model results_log ''fr''\n');
    fprintf('DESCRIPTION\n');
    fprintf('  ''file_list'' is a text file containing original and processed file\n');
    fprintf('                names in pairs, one pair on each line.  Paths are okay.\n');
    fprintf('                Each file must be in big-YUV format.  After the second file\n');
    fprintf('                name, optional (i.e., manual) calibration values may be listed.\n');
    fprintf('                RR_ORIGINAL does not use these values.  See "EXAMPLE LIST".\n');
    fprintf('                Optional calibration values are listed in the following order:\n');
    fprintf('         luma_gain luma_offset horiz_shift vert_shift delay\n'); 
    fprintf('                 luma_gain is luminance gain, double precision\n');
    fprintf('                 luma_offset is luminance offset, double precision\n');
    fprintf('                 horiz_shift is horizontal shift, integer; positive means\n');
    fprintf('                      processed has been moved right with respect to original\n');
    fprintf('                 vert_shift is vertical shift, integer; positive means\n');
    fprintf('                      processed has been moved down with respect to original\n');
    fprintf('                      Odd values mean that the processed video has been\n');
    fprintf('                      reframed (i.e., 1st field in time of original, aligns\n');
    fprintf('                      with 2nd field in time of processed -- i.e. +0.5 frame\n');
    fprintf('                      delay)\n');
    fprintf('                 delay is time delay in frames, integer; this value adjusts\n');
    fprintf('                       the start frame used by RR_PROCESSED -- it adds "delay"\n');
    fprintf('                       to the 0.8sec starting frame for the processed segment\n');
    fprintf('                       used.\n');
    fprintf('  ''video_standard'' indicates the frame rate and video size:\n');
    fprintf('    ''525''          525-line, 30fps video (720 pixels by 486 rows), "NTSC"\n');
    fprintf('                   Interlaced fields, lower field presented earlier in time\n');
    fprintf('    ''625''          625-line, 25fps video (720 pixels by 576 rows), "PAL"\n');
    fprintf('                   Interlaced fields, upper field presented earlier in time\n');
    fprintf('  ''model''          The name of the video quality model desired.  Must\n');
    fprintf('                   be one of the following:\n');
    fprintf('    ''lowbw''        Low Bandwidth Model\n');
    fprintf('    ''fastlowbw''    Fast Low Bandwidth Model, ITU-T Recommendation J.244\n');
    fprintf('    ''general''      General Model (FR-TV Phase II), ITU-T Recommendation J.144\n');
    fprintf('    ''developers''   Developers Model (approximates the FR-TV Phase II model)\n\n');
    fprintf('  ''results_log'' is the prefix (with path) for text files, where results will be\n');
    fprintf('                  written.  If results_log is ''c:\\temp\\525log'', then \n');
    fprintf('                  errors will be appended to ''c:\\temp\\525log_error.txt''\n');
    fprintf('  ''fr''         Optional flag, indicating that "full reference bandwidth" should\n');
    fprintf('               be used for calibration features.  Intended for validation.\n');
    fprintf('\n');
    fprintf('  Compressed reduced reference calibration and model features will be written to\n');
    fprintf('  files named after the original video sequence.  The calibration features will\n');
    fprintf('  have "_calibration.mat" appended, in the directory that contains that original\n');
    fprintf('  video sequence.  The model features will have "_features.dat" appended, in the\n');
    fprintf('  directory that contains the original video sequence.\n');
    fprintf('EXAMPLE CALL:\n');
    fprintf('  rr_original ''list_525.txt'' ''525'' ''lowbw'' ''log525''\n');
    fprintf('  rr_original ''list_625.txt'' ''625'' ''fastlowbw'' ''log625''\n');
    fprintf('EXAMPLE LIST:\n');
    fprintf('  c:\\v525\\SRC_08__525.yuv  c:\\v525\\SRC_08_MPEG2_m2@1000_525.yuv\n');
    fprintf('  c:\\v525\\SRC_15__525.yuv  c:\\v525\\SRC_15_v1@3000_525.yuv\n');
    fprintf('  c:\\SRC_16__525.yuv  c:\\SRC_16_H264_h4@6000_525.yuv 1.01 -3.2 +1 -1 2\n');
    fprintf('RESTRICTIONS:\n');
    fprintf('  All video sequences must be exactly 8-seconds in duration.\n');
    fprintf('\n');
    fprintf('  Test plan and model constraints taken together produces a maximum temporal\n');
    fprintf('  segment of 7-seconds for VQM and calibration.  The first 0.8 sec\n');
    fprintf('  and last 0.2s of the original video sequence will be ignored.\n');
    fprintf('  This software assumes valid video for the following region:\n');
    fprintf('     525-line/NTSC: top=21, left=31, bottom=466, right=690\n');
    fprintf('     625-line/PAL: top=21, left=31, bottom=556, right=690\n');
    fprintf('  RRNR-TV test plan constraints demand that the random algorithms use\n');
    fprintf('  a pseudo-random sequence (e.g., output the same VQM score when run twice).\n');
    fprintf('  This impacts the FastLowbw model and spatial shift registration.\n');
    fprintf('NOTES:\n');
    fprintf('  The models and calibration used by RR_ORIGINAL and RR_PROCESSED are unchanged\n');
    fprintf('  from those previously released to the public (i.e., CVQM, BVQM).\n');
    fprintf('  Source code and binary executables are available for download at\n');
    fprintf('  www.its.bldrdoc.gov  These algorithms can be freely used for commercial\n');
    fprintf('  and non-commercial applications.\n');
    fprintf('\n');
    return;
end


% strip off the extra single quotes '' for Windows compile, comment these
% eval lines out for Linux compile
file_list = eval(file_list); 
video_standard = eval(video_standard);
model = eval(model);
results_log = eval(results_log);

error_log = [results_log '_error.txt'];

% read list of original video files to be processed.
fid=fopen(file_list,'r');
warning off;
control_list = textscan(fid, '%s %s %f %f %d %d %d');
warning on;
fclose(fid);

if strcmp(video_standard, '525'),
    video_standard = 'interlace_lower_field_first';
    rows = 486;
    cols = 720;
    fps = 30;
    is_ntsc = 1;
elseif strcmp(video_standard, '625'),
    video_standard = 'interlace_upper_field_first';
    rows = 576;
    cols = 720;
    fps = 25;
    is_ntsc = 0;
else
    error('rr_original called with an invalid string for ''video_standard''.  Aborting');
end

% note whether to quantize 
is_quantize = 1;
if exist('flag','var'),
    flag = eval(flag);
    if strcmp(flag,'fr'),
        is_quantize = 0;
    else
        cvqm_error(error_log, 1, 'Flag on command line not recognized.  Ignored');
    end
end

% loop through all unique original video sequences.
original_list = unique(control_list{:,1});
for loop=1:length(original_list),
    original_file = original_list{loop};
    
    % note lack of quantizer!
    if ~is_quantize,
        fid = fopen(error_log,'a');
        fprintf(fid, '\r\nRR_ORIGINAL: %s %s\r\n', original_file);
        fprintf(fid, '  Warning: uncompressed original calibration features used\r\n');
        pause(0.2);
        fclose(fid);
    end

    % create names for output files & temporary file
    data_file = sprintf('%s_calibration.mat', original_file);
    model_file = sprintf('%s_features.dat', original_file);
    model_mat_file = sprintf('%s_features.mat', original_file);

    % write name of these files to the error log.
    fid = fopen(error_log,'a');
    fprintf(fid, 'RR_ORIGINAL: %s\r\n', original_file);
    error_ftell = ftell(fid);
    pause(0.2);
    fclose(fid);

    % check model requested is available
    if strcmp(model,'lowbw')
        % fine!
        is_model = 'l';
    elseif strcmp(model,'fastlowbw'),
        % fine!
        is_model = 'f';
    elseif strcmp(model,'general'),
        % fine!
        is_model = 'g';
    elseif strcmp(model,'developers'),
        % fine!
        is_model = 'd';
    else
        cvqm_error(error_log, 1, 'Model request string not recognized.');
        return;
    end
    

    try
        if ~exist(original_file,'file'),
            cvqm_error(error_log, 3, sprintf('Fatal Error. Original video file "%s" does not exist.  Clip skipped.\r\n', original_file));
            continue;
        end
%     fprintf('COMMENT IN TRY!\r\n\n');

        % initialize file read
        dll_video('initialize', 1, original_file, 'uyvy', video_standard, rows, cols, fps);
        [code8s] = dll_video('delay_8s', 1, 0);
        dll_calib_video('initialize', 1);
        
        if code8s == 1,
            % print warning.
            cvqm_error(error_log, 2, sprintf('Original video file "%s" is longer than 8sec. First 8sec used, only.\r\n', original_file));
        elseif code8s == 2,
            % abort!  cannot use this clip.
            cvqm_error(error_log, 3, sprintf('Fatal Error. Original video file "%s" is shorter than 8sec. Clip skipped.\r\n', original_file));
            continue;
        end

        % mandate seconds of video to be used.
        [total_sec] = 7;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calibration
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % default low bandwidth valid region
        [vr] = dll_default_vr(1);

        % perform low bandwidth temporal registration
        [ti2_orig_A, ti10_orig_A, ymean_orig_A, orig_is_white_clip, orig_is_black_clip] = ...
            dll_lowbw_temporal_features(1, total_sec, vr);

        % quantize
        if is_quantize,
            [ti2_index_A, ti10_index_A, y_index_A] = ...
              dll_lowbw_temporal_quant(1, ti2_orig_A, ti10_orig_A, ymean_orig_A);
        end

        % note white & black level clipping, if any
        if orig_is_white_clip,
            cvqm_error(error_log, 2, 'White level clipping detected on original video sequence may cause VQM errors.');
        end
        if orig_is_black_clip,
            cvqm_error(error_log, 2, 'Black level clipping detected on original video sequence may cause VQM errors.');
        end

        % calculate random seed, to initialize spatial shift & scaling
        % preferably this should be randomly seeded using the current time,
        % as per the following function:
        % % [seed_state] = dll_lowbw_calib_initialize;
        % however, this program needs to be predictable (same results each
        % time) so instead, the following altorithm will be used:  use
        % pixel (200,200) from the luminance plane of the first original image.
        dll_video('set_rewind', 1);
        y = dll_video('sec', 1, 0, 1.0/fps);
        dll_video('rewind', 1);
        seed_state = y(200,200);
% fprintf('Seed state = %d\n', seed_state);
        clear y;
        
        % calculate shift, & scaling 
        [orig_pixels, orig_horiz_profile, orig_vert_profile] = dll_lowbw_calib_original(1, seed_state, total_sec);

        % to quantize:
        if is_quantize,
            [index_horiz_profile, index_vert_profile] = ...
              dll_lowbw_calib_quant(1, orig_horiz_profile, orig_vert_profile);
        end


        % perform low bandwidth valid region
        [ovr] = dll_orig_valid_region;

        % estimate gain & offset
        [orig_y, orig_cb, orig_cr, yesno] = dll_lowbw_gain_v2_original(1, total_sec);

        % quantize gain & offset features
        if is_quantize,
            [index_y, index_cb, index_cr] = dll_lowbw_gain_v2_quant(1, orig_y, orig_cb, orig_cr); % quantize
        end

        % use constant PVR
        cvr.top = 21;
        cvr.left = 31;
        cvr.bottom = rows - 20;
        cvr.right = cols - 30;

        % input improved PVR estimate 
        dll_calib_video('calibration', 0, 0, cvr, ...
            1, 0, 1000, 1000);

        % perform low bandwidth temporal registration a second time!
        [ti2_orig_B, ti10_orig_B, ymean_orig_B, orig_is_white_clip, orig_is_black_clip] = ...
            dll_lowbw_temporal_features(1, total_sec, cvr);

        % quantize
        if is_quantize,
            [ti2_index_B, ti10_index_B, y_index_B] = ...
              dll_lowbw_temporal_quant(1, ti2_orig_B, ti10_orig_B, ymean_orig_B);
        end


        if is_quantize,
            eval(sprintf('save %s ti2_index_A ti10_index_A y_index_A orig_pixels index_horiz_profile index_vert_profile ovr cvr ti2_index_B ti10_index_B y_index_B seed_state total_sec cvr index_y index_cb index_cr yesno is_model is_quantize is_ntsc', data_file) );
        else
            eval(sprintf('save %s ti2_orig_A ti10_orig_A ymean_orig_A orig_pixels orig_horiz_profile orig_vert_profile ovr cvr ti2_orig_B ti10_orig_B ymean_orig_B seed_state total_sec cvr orig_y orig_cb orig_cr yesno is_model is_quantize is_ntsc', data_file) );
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % done calculating original video calibration statistics
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        if strcmp(model,'lowbw'),

            dll_features('Low', 1, total_sec, model_file);

        elseif strcmp(model,'fastlowbw'),

            dll_features('Fast', 1, total_sec, model_file);

        elseif strcmp(model,'general'),

            orig_features = dll_features('General', 1, total_sec);
            eval( sprintf('save %s orig_features', model_mat_file));

        elseif strcmp(model,'developers'),

            orig_features = dll_features('Developers', 1, total_sec);
            eval( sprintf('save %s orig_features', model_mat_file));

        end

        % note if no errors happened
        fid = fopen(error_log,'a');
        if error_ftell == ftell(fid),
            fprintf(fid, '  No errors encountered.\r\n');
        end
        pause(0.2);
        fclose(fid);

    catch
        cvqm_error(error_log, 1, lasterr);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cvqm_error(file_name, code, message);
% cvqm_error(file_name, code, message); 
%   % Append an error message to file 'file_name'.
%   % each line should start with a number (in "code") and then contain
%   % message describing the problem or issue. 
%   % if code == 0, delete the previous file (if any) and ignore message.
%   % code=1 	Data Input/Output invalid.  Operation impossible.
%   % code=2    Calibration values should be examined; a problem may exist.
%   % code=3    Fatal error


% open file & remember pointer
if code == 0,
    try
        warning off;
        delete(file_name);
        warning on;
    catch
    end
    return;
end

fid = fopen(file_name,'a');
fprintf(fid, '%d %s\r\n', code, message);
fclose(fid);








