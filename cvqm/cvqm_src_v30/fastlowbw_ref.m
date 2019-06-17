function fastlowbw_ref (original_file, processed_file, file_type, calibration);
% wrapper for fastlowbw_ref.  This is not the typical MATLAB interface, but rather
% an alternate interface required for compilation.
% 1) call with no arguments for help information.
% 2) if running inside MATLAB, add three quotes (''') both before and after
% every argument. 

if nargin == 0,
    fprintf('FASTLOWBW_REF -- Version 1.1\n');
    fprintf('  Take original and processed test sequences in raw BIG-YUV \n');
    fprintf('  file format.  Calculate NTIA Fast Low Bandwidth model and calibration.\n');
    fprintf('  This is reference software.\n');
    fprintf('SYNTAX\n');
    fprintf('  fastlowbw_ref original_file processed_file video_standard calibration model\n');
    fprintf('DESCRIPTION\n');
    fprintf('  Takes the name of a UYVY formatted file (''original_file'') containing an\n');
    fprintf('  original video sequence and UYVY formatted file (''processed_file'') containing\n');
    fprintf('  a processed video sequence.  Both files must contain uncompressed video in\n');
    fprintf('  either the UYVY or RGB color space.\n');
    fprintf('  ''video_standard'' indicates the frame rate and video size:\n');
    fprintf('    ''525''          525-line, 30fps video (720 pixels by 486 rows), "NTSC"\n');
    fprintf('                   Interlaced fields, lower field presented earlier in time\n');
    fprintf('    ''625''          625-line, 25fps video (720 pixels by 576 rows), "PAL"\n');
    fprintf('                   Interlaced fields, upper field presented earlier in time\n');
    fprintf('    ''VGA30''        VGA (480 lines x 640 pixels), 30fps video, progressive.\n');
    fprintf('    ''VGA25''        VGA (480 lines x 640 pixels), 25fps video, progressive.\n');
    fprintf('    ''CIF30''        CIF (288 lines x 352 pixels), 30fps video, progressive.\n');
    fprintf('    ''CIF25''        CIF (288 lines x 352 pixels), 25fps video, progressive.\n');
    fprintf('    ''QCIF30''        CIF (144 lines x 176 pixels), 30fps video, progressive.\n');
    fprintf('    ''QCIF25''        CIF (144 lines x 176 pixels), 25fps video, progressive.\n');
    fprintf('  ''calibration'' indicates the calibration options desired, and must\n');
    fprintf('  be one of the following:\n');
    fprintf('    ''none''       No calibration will be performed.\n');
    fprintf('    ''rrcal2''      Reduced Reference Bandwidth Calibration Version 2.0 (J.244)\n');
    fprintf('    ''rrcal2scale'' Reduced Reference Bandwidth Calibration Version 2.0 (J.244),\n');
    fprintf('                  including estimation of spatial scaling\n');
    fprintf('EXAMPLE CALL:\n');
    fprintf('  fastlowbw_ref ''original.yuv'' ''processed.yuv'' ''525'' ''rrcal2''\n');
    fprintf('RESTRICTIONS:\n');
    fprintf('  If the video sequences (after calibration) are longer than 15 seconds, then\n');
    fprintf('  only the first 15 seconds will be used for model calculation\n');
    fprintf('  Temporal registration uncertainty will +/- 1 sec.\n');
    fprintf('NOTES:\n');
    fprintf('  These algorithms are unchanged from those previously released to the public\n');
    fprintf('  (i.e., CVQM, BVQM).  Source code and binary executables are available for\n');
    fprintf('  download at www.its.bldrdoc.gov  These algorithms can be freely used for\n');
    fprintf('  commercial and non-commercial applications.\n');
    
    return;
end

% NOTES FOR PROGRAMMERS:
% 
% 1. To modify the temporal delay search range, change variable "uncert"
% passed into function dll_lowbw_temporal.m
% 
% 2. To modify spatial shift and scaling search limits, change variables
% "max_shift_horiz", "max_shift_vert", "max_scale_horiz" and
% "max_scale_vert" at the beginning of function dll_lowbw_calib_processed.m
% and dll_lowbw_calib_original.m.  These values must be identical for both
% functions, and are currently set automatically based on image size.
% Choosing larger values may make the algorithm unstable. 
%
% 3. The original file is referred to by an ID of "1" throughout, and the 
% processed file is referred to by an ID of "2".  These are referred to
% within functions by the variable name "fn" (i.e., "F" for "file" and "N" for
% "number"). 
%
% 4. To separate into upstream and downstream parts, simply call each
% function requiring the original file ("1") into one piece, and the
% functions requiring the processed file ("2") into another piece.  This
% code presumes downstream monitoring or bi-directional communication. 
% Some of the functions require both processed video file and the original
% features and results from the previous steps. 
% 
% 5. If implementing a system that is entirely down-stream (i.e., no
% communication from processed to original), then the original will need to
% assume a valid region.  Discarding the overscan is recommended (see also
% dll_default_vr.m). 
%
% 6. The reduced reference calibration routines and model presented herein
% use some randomized processies.  As a result, the results of this program
% will differ slightly from one run to the next.  This variability can be
% prevented by initializeing the random number generator with a constant
% value.  This will aid debugging. 
%
% 7. When comparing results from this implementation to results from 
% BVQM, slight differences will be obverved.  These stem from the random
% number generator, yet also because BVQM performs filtering over the
% results from multiple video sequences.  The places where filtering should
% be performed is mentioned in comments, below. 

  
        
% strip off the extra single quotes '' for Windows compile, comment these
% eval lines out for Linux compile
original_file = eval(original_file); 
processed_file = eval(processed_file); 
file_type = eval(file_type);
calibration = eval(calibration); 

% create names for output files & temporary file
temporary_file = sprintf('%s_temp.txt', processed_file);
model_file = sprintf('%s_model.txt', processed_file);
calibration_file = sprintf('%s_calibration.txt', processed_file);
error_file = sprintf('%s_errors.txt', processed_file);
cvqm_error(error_file, 0);

if strcmpi(file_type, '525'),
    video_standard = 'interlace_lower_field_first';
    rows = 486;
    cols = 720;
    fps = 30;
elseif strcmpi(file_type, '625'),
    video_standard = 'interlace_upper_field_first';
    rows = 576;
    cols = 720;
    fps = 25;
elseif strcmpi(file_type, 'VGA30'),
    video_standard = 'progressive';
    rows = 480;
    cols = 640;
    fps = 30;
elseif strcmpi(file_type, 'VGA25'),
    video_standard = 'progressive';
    rows = 480;
    cols = 640;
    fps = 25;
elseif strcmpi(file_type, 'CIF30'),
    video_standard = 'progressive';
    rows = 288;
    cols = 352;
    fps = 30;
elseif strcmpi(file_type, 'CIF25'),
    video_standard = 'progressive';
    rows = 288;
    cols = 352;
    fps = 25;
elseif strcmpi(file_type, 'QCIF30'),
    video_standard = 'progressive';
    rows = 144;
    cols = 176;
    fps = 30;
elseif strcmpi(file_type, 'QCIF25'),
    video_standard = 'progressive';
    rows = 144;
    cols = 176;
    fps = 25;
else
    error('fastlowbw_ref called with an invalid string for ''video_standard''.  Aborting');
end

if strcmpi(calibration,'rrcal2') || strcmpi(calibration,'rrcal2scale') || strcmpi(calibration,'none'),
    % fine!
else
    cvqm_error(error_file, 1, 'Calibration request string not recognized.');
    return;
end

model = 'fastlowbw';

try
    % for debugging, you'll want to comment out this try/catch loop.
    % the following print is to remind you to re-insert the try/catch back.
    % fprintf('COMMENT IN TRY!\n\n');

    % initialize file read
    dll_video('initialize', 1, original_file, 'uyvy', video_standard, rows, cols, fps);
    dll_video('initialize', 2, processed_file, 'uyvy', video_standard, rows, cols, fps);
    dll_calib_video('initialize', 1);

    % calculate seconds of video to be used.
    [total_sec] = min( dll_video('total_sec',1), dll_video('total_sec',2) );
    if total_sec > 15,
        total_sec = 15;
        cvqm_error(error_file, 3, 'File is longer than 15 seconds.');
        cvqm_error(error_file, 3, 'Results will be calculated using first 15 seconds only.');
    end
    if total_sec < 4,
        cvqm_error(error_file, 3, 'Video files are too brief.  Aborting.');
        return;
    end

    warning off;
    delete(calibration_file);
    warning on;

    if strcmpi(calibration, 'none'),
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use "no calibration" values
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        shift.horizontal = 0;
        shift.vertical = 0;
        pvr = dll_default_vr(1);
        y_gain = 1.0;
        y_offset = 0.0;
        scale.horizontal = 1000;
        scale.vertical = 1000;
        delay = 0;

        % We don't need to do anything with the calibration values (above), 
        % because these are the defaults set by dll_calib_video.m.  

        % write results.
        [sucess] = cvqm_save_calibration(calibration_file, calibration, ...
            shift, pvr, y_gain, y_offset, scale, delay); 
        if sucess == 0,
            cvqm_error(error_file, 1, 'Cannot open file to write calibration results.');
            return;
        end
        
    else
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Start Calibration
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Determine whether or not to estimate scaling.
        if strcmpi(calibration,'rrcal2'),
            no_scaling = 1;
        elseif strcmpi(calibration,'rrcal2scale'),
            no_scaling = 0;
        end

        total_sec = floor(total_sec);

        % Assume the default low bandwidth valid region as a starting point
        [vr] = dll_default_vr(1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Initial temporal registration
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % calcualte low bandwidth temporal registration features on the
        % original video sequence.
        [ti2_orig, ti10_orig, ymean_orig, orig_is_white_clip, orig_is_black_clip] = ...
            dll_lowbw_temporal_features(1, total_sec, vr);

        % Quantize the original temporal registration features and encode.
        % Return variables are suitable for low bandwidth transmission 
        % (e.g., when saved to a *.mat file).  This save/load step is not 
        % demonstrated.
        [ti2_index, ti10_index, y_index] = ...
            dll_lowbw_temporal_quant(1, ti2_orig, ti10_orig, ymean_orig);

        % Reverse: go from encoded variables back into quantized original features.
        [ti2_orig, ti10_orig, ymean_orig] = ...
            dll_lowbw_temporal_quant(0, ti2_index, ti10_index, y_index);

        % calculate low bandwidth temporal registration features on the
        % processed video sequence.
        [ti2_proc, ti10_proc, ymean_proc, proc_is_white_clip, proc_is_black_clip] = ...
            dll_lowbw_temporal_features(2, total_sec, vr);

        % Print errors and warnings.  Note white & black level clipping, if any
        if orig_is_white_clip,
            cvqm_error(error_file, 2, ...
                'White level clipping detected on original video sequence may cause VQM errors.');
        end
        if orig_is_black_clip,
            cvqm_error(error_file, 2,...
                'Black level clipping detected on original video sequence may cause VQM errors.');
        end
        if proc_is_white_clip,
            cvqm_error(error_file, 2, ...
                'White level clipping detected on processed video sequence may cause VQM errors.');
        end
        if proc_is_black_clip,
            cvqm_error(error_file, 2, ...
                'Black level clipping detected on processed video sequence may cause VQM errors.');
        end


        % Compute temporal registration, using original and processed features.
        % variable "uncert" is the uncertainty in seconds that should be
        % searched by this temporal registration algorithmm.
        uncert = 1;
        [delay, sucess, is_still] = ...
            dll_lowbw_temporal (1, ti2_orig, ti2_proc, ti10_orig, ti10_proc, ymean_orig, ...
            ymean_proc, uncert, 'field');

        if sucess == 0 && is_still,
            cvqm_error(error_file, 2, ...
                'Still or nearly still sequence; initial temporal registration cannot be computed.');
        elseif sucess == 0,
            cvqm_error(error_file, 2, 'Initial temporal registration algorithm failed.');
        end

        % note if re-framing is indicated
        if mod(delay,1) == 0.5,
            dll_calib_video('set_reframe',1);
            delay = delay - 0.5;
        end

        % re-calculate seconds of video to be used (may have changed)
        [total_sec] = min( dll_calib_video('total_sec',1), dll_calib_video('total_sec',2) );
        if total_sec > 15,
            total_sec = 15;
        end

        % apply the temporal registration to all future video read calls.
        dll_lowbw_temporal_original(1, delay);
        dll_lowbw_temporal_processed(2, delay);

        [total_sec] = min( dll_calib_video('total_sec',1), dll_calib_video('total_sec',2) );
        if total_sec > 15,
            total_sec = 15;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Spatial Shift and [optional] spatial scaling
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [seed_state] = dll_lowbw_calib_initialize;

        % compute original features.
        [orig_pixels, orig_horiz_profile, orig_vert_profile] = ...
            dll_lowbw_calib_original(1, seed_state, total_sec);

        % Quantize the original spatial shift/scaling features and encode.
        % Return variables are suitable for low bandwidth transmission 
        % (e.g., when saved to a *.mat file).  This save/load step is not 
        % demonstrated.
        [index_horiz_profile, index_vert_profile] = ...
            dll_lowbw_calib_quant(1, orig_horiz_profile, orig_vert_profile);

        % Reverse: go from encoded variables back into quantized original features.
        [orig_horiz_profile, orig_vert_profile] = ...
            dll_lowbw_calib_quant(0, index_horiz_profile, index_vert_profile);

        % Compute processed features   
        [shift, scale,status] = ...
            dll_lowbw_calib_processed(2, seed_state, total_sec, orig_pixels, ...
            orig_horiz_profile, orig_vert_profile, no_scaling);

        % print errors and warnings.
        if status.scale && ~no_scaling,
            cvqm_error(error_file, 2, 'Actual spatial scaling may be beyond search limits.');
        end
        if status.shift,
            cvqm_error(error_file, 2, 'Actual spatial shift may be beyond search limits.');
        end
        if (strcmp(video_standard,'interlace_lower_field_first') || ...
                strcmp(video_standard,'interlace_upper_field_first')) && mod(shift.vertical,2),
            cvqm_error(error_file, 2, 'Processed video will be reframed.');
            cvqm_error(error_file, 2, 'Actual temporal delay is 0.5 frames greater than reported value.');
        end
        if abs(shift.horizontal) > 8 || abs(shift.vertical) > 5,
            cvqm_error(error_file, 2, 'Extreme spatial shift detected.');
        end
        if scale.vertical ~= 1000 || scale.horizontal ~= 1000,
            cvqm_error(error_file, 2, 'Video scaling detected; please examine other scenes.');
        end


        % Apply above estimates, to all future video read calls. 
        %
        % If 2+ video sequences are available for the same system, ideally the
        % spatial shift and spatial scaling values would be filtered across the 
        % results from all video sequences (e.g., sort horizontal shifts and use
        % the median horizontal shift for all video sequences).  This increases 
        % the shift and scaling estimation accuracy, and is particularly
        % important for scaling factors.
        %
        % Since this reference code uses one original/processed video 
        % sequence pair, this cannot be demonstrated.
        y_gain = 1.0;
        y_offset = 0.0;
        dll_calib_video('calibration', shift.horizontal, shift.vertical, vr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);

        % re-calculate seconds of video to be used (may have changed)
        [total_sec] = min( dll_calib_video('total_sec',1), dll_calib_video('total_sec',2) );
        if total_sec > 15,
            total_sec = 15;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Valid Video Estimation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % calculate original video valid region.  This returns 4 integers, so
        % quantization is unnecessary.  
        [ovr] = dll_orig_valid_region;

        % calculate processed video valid region.
        [pvr] = dll_proc_valid_region(ovr);

        % print errors and warnings.
        if (pvr.bottom - pvr.top + 1) / rows < 0.55 || ...
                (pvr.right - pvr.left + 1) / cols < 0.80,
            cvqm_error(error_file, 2, 'Greatly reduced valid region detected.');
        end

        % Apply processed valid region estimate to all future video read calls.
        %
        % If 2+ systems will be compared, ideally the same valid region should
        % be used for that scene for all systems (PVSs).  Choose the smallest
        % valid region, and apply to all PVSs.  
        %
        % Since this reference code uses one original/processed video 
        % sequence pair, this cannot be demonstrated.

        dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Luminance, Cb and Cr gain and level offset estimation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Calculate low bandwidth gain/offset features on the original video.
        [orig_y, orig_cb, orig_cr, yesno] = dll_lowbw_gain_v2_original(1, total_sec);

        % Quantize the original gain/offset features and encode.
        % Return variables are suitable for low bandwidth transmission 
        % (e.g., when saved to a *.mat file).  This save/load step is not 
        % demonstrated.
        [index_y, index_cb, index_cr] = dll_lowbw_gain_v2_quant(1, orig_y, orig_cb, orig_cr); % quantize

        % Reverse: go from encoded variables back into quantized original features.
        [orig_y, orig_cb, orig_cr] = dll_lowbw_gain_v2_quant(0, index_y, index_cb, index_cr); % reconstruct

        % calculate Y, Cb, and Cr gain and level offset from the original video
        % features and the processed video sequence file.
        [y_gain, y_offset, cb_gain, cb_offset, cr_gain, cr_offset, sucess] = ...
            dll_lowbw_gain_v2_processed(2, total_sec, orig_y, orig_cb, orig_cr, yesno);
        cvqm_error(error_file, 2, 'Color Gain & Offset estimated but not removed.');

        % print errors and warnings. 
        if y_gain < 0.9 || y_gain > 1.1,
            cvqm_error(error_file, 2, 'Extreme Luminance Gain detected.');
        end
        if y_offset < -20 || y_offset > 20,
            cvqm_error(error_file, 2, 'Extreme Luminance Offset detected.');
        end
        if sucess == 0,
            cvqm_error(error_file, 2, ...
                'Warning: algorithm used to estimate Luminance Gain & Offset may have failed.');
        end
        if sucess == -1,
            cvqm_error(error_file, 2, ...
                'Luminance gain & offset algorithm detected extreme values or failed; results discarded.');
        end
        if isnan(cb_gain),
            cvqm_error(error_file, 2, 'Warning: algorithm used to estimate Cb Gain & Offset failed.');
        end
        if isnan(cr_gain),
            cvqm_error(error_file, 2, 'Warning: algorithm used to estimate Cr Gain & Offset failed.');
        end

        % Apply the luminance gain & offset estimates to all future video read
        % calls.  Do NOT apply the Cb and Cr estimates -- these are considered
        % errors that the viewer may object to.  
        %
        % If 2+ video sequences are available for the same system, ideally the
        % luma gain and offset values would be filtered across the results from
        % all video sequences (e.g., sort values and use the median gain for
        % all video sequences).  This increases the luma gain/offset accuracy.
        %
        % Since this reference code uses one original/processed video 
        % sequence pair, this cannot be demonstrated.
        dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Second temporal registration -- improved estimate
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % calcualte low bandwidth temporal registration features on the
        % original video sequence.
        [ti2_orig, ti10_orig, ymean_orig, orig_is_white_clip, orig_is_black_clip] = ...
            dll_lowbw_temporal_features(1, total_sec, pvr);

        % Quantize the original temporal registration features and encode.
        % Return variables are suitable for low bandwidth transmission 
        % (e.g., when saved to a *.mat file).  This save/load step is not 
        % demonstrated.
        [ti2_index, ti10_index, y_index] = ...
            dll_lowbw_temporal_quant(1, ti2_orig, ti10_orig, ymean_orig);

        % Reverse: go from encoded variables back into quantized original features.
        [ti2_orig, ti10_orig, ymean_orig] = ...
            dll_lowbw_temporal_quant(0, ti2_index, ti10_index, y_index);

        % calcualte low bandwidth temporal registration features on the
        % processed video sequence.
        [ti2_proc, ti10_proc, ymean_proc, proc_is_white_clip, orig_is_black_clip] = ...
            dll_lowbw_temporal_features(2, total_sec, pvr);

        % calculate temporal registration from the original and processed
        % features.
        uncert = 1;
        [delay2, sucess, is_still] = ...
            dll_lowbw_temporal (1, ti2_orig, ti2_proc, ti10_orig, ti10_proc, ymean_orig, ymean_proc, uncert);

        % print errors and warnings.
        if sucess == 0 && is_still,
            cvqm_error(error_file, 2, ...
                'Still or nearly still sequence; temporal registration cannot be computed.');
        elseif sucess == 0,
            cvqm_error(error_file, 2, 'Final temporal registration algorithm failed.');
        end

        % apply the improved delay estimate to all future video read calls.
        dll_lowbw_temporal_original(1, delay2);
        dll_lowbw_temporal_processed(2, delay2);

        % print errors and warnings. 
        if abs(delay+delay2) >= round(fps),
            cvqm_error(error_file, 2, 'Temporal mis-registration exceeds 1 second uncertainty limit.');
        end

        % write results.
        [sucess] = cvqm_save_calibration(calibration_file, calibration, shift, pvr, y_gain, ...
            y_offset, scale, delay+delay2, cb_gain, cb_offset, cr_gain, cr_offset); 
        if sucess == 0,
            cvqm_error(error_file, 1, 'Cannot open file to write lowbw calibration results.');
            return;
        end

    end


    % re-calculate seconds of video to be used (may have changed)
    [total_sec] = min( dll_calib_video('total_sec',1), dll_calib_video('total_sec',2) );
    if total_sec > 15,
        total_sec = 15;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate FastLowBW model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    warning off;
    delete(model_file);
    warning on;

    % calculate features on the original video sequences, and save to file
    % "temporary_file"
    dll_features('Fast', 1, total_sec, temporary_file);
    
    % calculate features on the processed video sequence, and hold in
    % variable "proc_features"
    proc_features = dll_features('Fast', 2, total_sec);
    
    % Compute video quality metric, and return in variable "VQM"
    [vqm, pars, par_names] = dll_model('vqm',temporary_file, proc_features);
    
    % delete the temporary file with original features.
    delete(temporary_file);
    
    % Save model results to file.
    [sucess] = cvqm_save_model(model_file, model, vqm, pars, par_names);

    
catch
    cvqm_error(error_file, 1, lasterr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sucess] = cvqm_save_calibration(file_name, calibration,...
    shift, pvr, y_gain, y_offset,scale, delay, cb_gain, cb_offset, cr_gain, cr_offset);
% [sucess] = cvqm_save_calibration(file_name, calibration, shift, pvr, y_gain, y_offset,scale, delay); 
%   % Open for writing 'file_name' & record type of calibration requested.
%   % 'sucess' is 1 if this function can write; 0 if fail to write.
%   % Record to file values produced by calibration 


% open file & remember pointer
fid = fopen(file_name,'w');
if fid <= 0,
    sucess = 0;
    return;
else
    sucess = 1;
end

% write out type of calibration performed.
fprintf(fid, '%s\r\n', calibration);

% write out calibration values
fprintf(fid,'%5d  Horizontal Shift\r\n', shift.horizontal);
fprintf(fid,'%5d  Vertical Shift\r\n', shift.vertical);
fprintf(fid,'%5d  Valid Region Top\r\n', pvr.top);
fprintf(fid,'%5d  Valid Region Left\r\n', pvr.left);
fprintf(fid,'%5d  Valid Region Bottom\r\n', pvr.bottom);
fprintf(fid,'%5d  Valid Region Right\r\n', pvr.right);
fprintf(fid,'%5.3f Luminance Gain\r\n', y_gain);
fprintf(fid,'%5.3f Luminance Offset\r\n', y_offset);
fprintf(fid,'%5d  Horizontal Scale\r\n', scale.horizontal);
fprintf(fid,'%5d  Vertical Scale\r\n', scale.vertical);
fprintf(fid,'%5d  Temporal Delay\r\n', delay);

if exist('cr_offset'),
    fprintf(fid,'\r\n%5.3f Cb Gain\r\n', cb_gain);
    fprintf(fid,'%5.3f Cb Offset\r\n', cb_offset);
    fprintf(fid,'%5.3f Cr Gain\r\n', cr_gain);
    fprintf(fid,'%5.3f Cr Offset\r\n', cr_offset);
end

fclose(fid);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sucess, shift, pvr, y_gain, y_offset,scale, delay] = cvqm_load_calibration(file_name);
% [sucess, shift, pvr, y_gain, y_offset,scale, delay] = cvqm_save_calibration(file_name); 
%   % Open 'file_name' & read calibration values.
%   % 'sucess' is 2 if this function can read sucessfully; 0 if failed.
%   % 'sucess' is 1 if values look to be unreasonably extreme.


% open file & remember pointer
fid = fopen(file_name,'r');
if fid <= 0,
    sucess = 0;
    return;
end
sucess = 2;

% write out type of calibration performed.
fgets(fid);

% write out calibration values
shift.horizontal = fscanf(fid,'%d');
fgets(fid);
shift.vertical = fscanf(fid,'%d');
fgets(fid);
pvr.top = fscanf(fid,'%d');
fgets(fid);
pvr.left = fscanf(fid,'%d');
fgets(fid);
pvr.bottom = fscanf(fid,'%d');
fgets(fid);
pvr.right = fscanf(fid,'%d');
fgets(fid);
y_gain = fscanf(fid,'%f');
fgets(fid);
y_offset = fscanf(fid,'%f');
fgets(fid);
scale.horizontal = fscanf(fid,'%d');
fgets(fid);
scale.vertical = fscanf(fid,'%d');
fgets(fid);
delay = fscanf(fid,'%d');
fgets(fid);

fclose(fid);

% error check
if shift.horizontal < -50 || shift.horizontal > 50 || shift.vertical < -50 || shift.vertical > 50,
    % shift really unreasonable.
    status = 1;
end
[rows,cols,fps,durration] = dll_video('size', 1);
if pvr.top < 1 || pvr.left < 1 || pvr.bottom > rows || pvr.right > cols,
    % pvr really unreasonable
    status = 1;
end
if y_gain < 0.2 || y_gain > 3.0 || y_offset < -100 || y_offset > 100,
    % gain and/or offset really unreasonable
    status = 1;
end
if scale.horizontal < 500 || scale.horizontal > 2000 || scale.vertical < 500 || scale.vertical > 2000,
    % scale really unreasonable.
    status = 1;
end
delay_sec = delay / fps;
if abs(delay_sec) > durration/2,
    % delay is longer than half the file durration!  really unreasonable.
    status = 1;
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sucess] = cvqm_save_model(file_name, model,vqm, pars, par_names);
% [sucess] = cvqm_save_model(file_name, model,vqm, pars, par_names); 
%   % Open for writing 'file_name' & record type of model requested.
%   % 'sucess' is 1 if this function can write; 0 if fail to write.
%   % Record to file values produced by model 


% open file & remember pointer
fid = fopen(file_name,'w');
if fid <= 0,
    sucess = 0;
    return;
else
    sucess = 1;
end

% write out type of model performed.
fprintf(fid, '%f %s\r\n', vqm, model);

% write out calibration values
for cnt=1:length(pars)
    fprintf(fid, '%f %s\r\n', pars(cnt), par_names{cnt});
end

fclose(fid);





