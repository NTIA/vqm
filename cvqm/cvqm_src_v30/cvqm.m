function cvqm (original_file, processed_file, video_standard, calibration, model, varargin)
% Call with no arguments for help information.

if nargin == 0,
    fprintf('CVQM -- Version 3.0\n');
    fprintf('  Given an original and processed video sequence in uncompressed AVI file format, perform\n');
    fprintf('  an NTIA video calibration option and then calculate an NTIA video quality model (VQM).\n');
    fprintf('\n');
    fprintf('SYNTAX\n');
    fprintf('  cvqm ''original_file'' ''processed_file'' ''video_standard'' ''calibration'' ''model''\n');
    fprintf('\n');
    fprintf('DESCRIPTION\n');
    fprintf('  Processes one original AVI formatted file (''original_file'') and a corresponding processed\n');
    fprintf('  AVI formatted file (''processed_file'') of the specified ''video_standard'', and performs\n');
    fprintf('  the requested video quality ''calibration'' and ''model''. The AVI formatted files must be\n');
    fprintf('  uncompressed YV12, YUY2, UYVY, V210, or RGB color space.\n');
    fprintf('\n');
    fprintf('  ''video_standard'' indicates field or frame ordering, and must be one of the following:\n');
    fprintf('    ''progressive''                  Progressive frames\n');
    fprintf('    ''interlace_lower_field_first''  Interlaced fields, lower field presented earlier in time\n');
    fprintf('    ''interlace_upper_field_first''  Interlaced fields, upper field presented earlier in time\n');
    fprintf('\n');
    fprintf('  ''calibration'' indicates the calibration option desired, and must be one of the following:\n');
    fprintf('    ''none''       No calibration will be performed.\n');
    fprintf('    ''manual''     Load calibration results from previous CVQM output.\n');
    fprintf('    ''rrcal''      Reduced Reference Bandwidth Calibration\n');
    fprintf('    ''rrcalscale'' Reduced Reference Bandwidth Calibration, including \n');
    fprintf('                 estimation of spatial scaling\n');
    fprintf('    ''rrcal2''      Reduced Reference Bandwidth Calibration Version 2.0 (ITU-T J.244)\n');
    fprintf('                  [Color Gain & Offset estimates; slightly improved Y gain & offset.]\n');
    fprintf('    ''rrcal2scale'' Reduced Reference Bandwidth Calibration Version 2.0 (ITU-T J.244),\n');
    fprintf('                 including estimation of spatial scaling\n');
    fprintf('    ''frcal''      Full Reference Bandwidth Calibration,\n');
    fprintf('                 The General Model''s ANSI & ITU standard calibration (ANSI T1.801.03, ITU-T J.144)\n');
    fprintf('    ''frtime''     Full Reference temporal registraion and valid region, only; no other\n');
    fprintf('                 calibration performed\n');
    fprintf('    ''rrtime''     Reduced Reference temporal registraion and valid region, only; no other\n');
    fprintf('                 calibration performed\n');
    fprintf('    ''frtimemanual''     Load calibration results from previous CVQM output.  Delay ignored;\n');
    fprintf('                       Full Reference temporal registration performed.\n');
    fprintf('    ''rrtimemanual''     Load calibration results from previous CVQM output.  Delay ignored;\n');
    fprintf('                       Reduced Reference temporal registration performed.\n');
    fprintf('    ''psnr_search''        Peak Signal to Noise Ratio (PSNR) exhaustive search based on ITU-T J.340,\n');
    fprintf('                         spatial search is +/- 3 pixels, temporal search is +/- 1 second.\n');
    fprintf('    ''fast_psnr_search1''   Reduced Reference Bandwidth Calibration Version 2.0 (ITU-T J.244),\n');
    fprintf('                         then PSNR exhaustive search based on ITU-T J.340 (spatial search\n');
    fprintf('                         is +/- 1 pixel, temporal search is +/- 0.5 seconds).\n');
    fprintf('    ''fast_psnr_search2''   Full Reference Bandwidth Calibration, (ANSI T1.801.03, ITU-T J.144),\n');
    fprintf('                          then PSNR exhaustive search based on ITU-T J.340 (spatial search\n');
    fprintf('                          is +/- 1 pixel, temporal search is +/- 0.5 seconds).\n');
    fprintf('\n');
    fprintf('  ''model'' is the video quality model (VQM) desired, and must be one of the following:\n');
    fprintf('    ''none''       No model will be calculated,\n');
    fprintf('    ''general''    NTIA General Model, ANSI T1.801.03 & ITU-T J.144 standard\n');
    fprintf('    ''developers'' Developers Model\n');
    fprintf('    ''lowbw''      Low Bandwidth Model\n');
    fprintf('    ''fastlowbw''  Fast Low Bandwidth Model, ITU-T J.249\n');
    fprintf('    ''psnr''       Peak Signal to Noise Ratio (PSNR) \n');
    fprintf('    ''psnr_vfd''   Peak Signal to Noise Ratio (PSNR) with variable frame delay (VFD) processing.\n');
    fprintf('                 This model option will recalculate gain and level offset after the VFD correction\n');
    fprintf('                 and the calibration file will be updated accordingly.\n');
    fprintf('    ''vqm_vfd''    Neural Network (NN) based VQM with variable frame delay (VFD) processing\n');
    fprintf('\n');
    fprintf('  Any or all of the following optional properties may be requested:\n');
    fprintf('  ''viewing_distance'' vd   This sets the viewing distance (vd) for the model vqm_vfd. The units\n');
    fprintf('                          for viewing distance are picture heights (H). A default viewing distance\n');
    fprintf('                          is set based on the number of image rows if option is not specified.\n');
    fprintf('                          Default viewing distances are:\n');
    fprintf('                          vd = 8: image_rows > 72 and image_rows < 216\n');
    fprintf('                          vd = 7: image_rows >= 216 && image_rows < 384\n');
    fprintf('                          vd = 5: image_rows >= 384 && image_rows < 648\n');
    fprintf('                          vd = 4: image_rows >= 648 && image_rows < 900\n');
    fprintf('                          vd = 3: image_rows >= 900 && image_rows < 1260\n');
    fprintf('\n');
    fprintf('EXAMPLE CALL:\n');
    fprintf('  cvqm ''original.avi'' ''processed.avi'' ''progressive'' ''none'' ''general'' \n');
    fprintf('\n');
    fprintf('RESTRICTIONS:\n');
    fprintf('  If the video sequences (after calibration) are longer than 15 seconds, then only the first 15\n');
    fprintf('  seconds will be used. Temporal registration uncertainty is set to +/- 1 sec. The number of\n');
    fprintf('  vertical pixels in the image must be between 72 and 1260.\n');
    return;
end

% strip off the extra single quotes ''
original_file = eval(original_file); 
processed_file = eval(processed_file); 
video_standard = eval(video_standard);
calibration = eval(calibration); 
model = eval(model);

% create names for output files & temporary file
temporary_file = sprintf('%s_temp.txt', processed_file);
model_file = sprintf('%s_model.txt', processed_file);
calibration_file = sprintf('%s_calibration.txt', processed_file);
error_file = sprintf('%s_errors.txt', processed_file);
cvqm_error(error_file, 0);

% Evaluate extra input arguments
cnt = 1;
while cnt <= length(varargin)
    if strcmpi(eval(char(varargin(cnt))),'viewing_distance')
        viewing_distance = str2double(varargin{cnt+1});
        if viewing_distance <= 2
            cvqm_error(error_file, 1, 'Viewing distance too small, must be greater than 2');
            return;
        elseif viewing_distance >= 12
            cvqm_error(error_file, 1, 'Viewing distance too large, must be less than 12');
            return;
        end
        cnt = cnt + 2;
    else
        cvqm_error(error_file, 1, 'Property value passed into cvqm not recognized');
        return;
    end
end

% look for command line errors.
if strcmp(video_standard,'progressive')|| strcmp(video_standard,'interlace_lower_field_first') || ...
        strcmp(video_standard,'interlace_upper_field_first'),
    % fine!
else
    cvqm_error(error_file, 1, 'Video Standard request string not recognized.');
    return;
end

if strcmp(calibration,'none')|| strcmp(calibration,'rrcal') || strcmp(calibration,'rrtime') || ...
        strcmp(calibration,'none')|| strcmp(calibration,'rrcal2') || ...
        strcmp(calibration,'rrcalscale') || strcmp(calibration,'rrcal2scale') || ...
        strcmp(calibration,'manual') || ...
        strcmp(calibration,'rrtimemanual') || strcmp(calibration,'frtimemanual') || ...
        strcmp(calibration,'frtime') || strcmp(calibration,'frcal') || ...
        strcmpi(calibration, 'psnr_search') || strcmpi(calibration, 'fast_psnr_search1') || strcmpi(calibration, 'fast_psnr_search2')
    % fine!
else
    cvqm_error(error_file, 1, 'Calibration request string not recognized.');
    return;
end

if strcmp(model,'general') || strcmp(model,'developers') || ...
	strcmp(model,'lowbw') || strcmp(model,'fastlowbw') || ...
    strcmp(model,'none') || strcmpi(model,'psnr') || strcmp(model, 'vqm_vfd') || strcmpi(model, 'psnr_vfd')
    % fine!
else
    cvqm_error(error_file, 1, 'Model request string not recognized.');
    return;
end

try
% fprintf('COMMENT IN TRY!\n\n');

    % initialize file read
    dll_video('initialize', 1, original_file, 'avi', video_standard);
    dll_video('initialize', 2, processed_file, 'avi', video_standard);
    dll_calib_video('initialize', 1);

    % error checks on above files -- does this make sense?
    [rows1,cols1,fps1] = dll_video('size', 1);
    [rows2,cols2,fps2] = dll_video('size', 2);
    if rows1 ~= rows2 || cols1 ~= cols2,
        cvqm_error(error_file, 1, 'Original and Processed video files contain different image sizes.');
        return;
    elseif fps1 ~= fps2,
        cvqm_error(error_file, 1, 'Original and Processed video files contain different frame rates.');
        return;
    end
    
    % Set viewing_distance (if not already set by the user)
    if ~(exist('viewing_distance', 'var') == 1)
        image_rows = rows1;
        if (image_rows > 72 && image_rows < 216)  % QCIF
            viewing_distance = 8;
        elseif (image_rows >= 216 && image_rows < 384)  % CIF
            viewing_distance = 7;
        elseif (image_rows >= 384 && image_rows < 648)  % VGA and SD
            viewing_distance = 5;
        elseif (image_rows >= 648 && image_rows < 900)  % HD 720
            viewing_distance = 4;
        elseif (image_rows >= 900 && image_rows < 1260)  % HD 1080
            viewing_distance = 3;
        else  % Unsupported resolution
            cvqm_error(error_file, 1, 'Unsupported image resolution (image_rows<=72 or image_rows>=1260).');
            return;
        end
    end

    % calculate seconds of video to be used.
    [total_sec] = min( dll_video('total_sec',1), dll_video('total_sec',2) );
    if total_sec > 15,
        total_sec = 15;
        cvqm_error(error_file, 3, 'File is longer than 15 seconds.  Results will be calculated using first 15 seconds only.');
    end
    if total_sec < 4,
        cvqm_error(error_file, 3, 'Video files is too brief.  CVQM requires at least 4-seconds of video to run.');
        return;
    end
    
    % If a PSNR calibration option is selected, set the fraction_sampled
    % value (based on the size of the video). The psnr_vfd model must be
    % included here also as the gain and offset must be recalculated for
    % this model (this calculation uses the fraction_sampled value).
    if strcmpi(calibration, 'psnr_search') || strcmpi(calibration, 'fast_psnr_search1') || strcmpi(calibration, 'fast_psnr_search2') || strcmpi(model, 'psnr_vfd')
        if cols1 <= 400
            fs = 1;
        elseif cols1 > 400 && cols1 <= 800
            fs = .3;
        else
            fs = .1;
        end
    end

    if strcmp(calibration,'none'),
        warning off;
        delete(calibration_file);
        warning on;

        % Use default "no calibration" values.
        shift.horizontal = 0;
        shift.vertical = 0;
        pvr = dll_default_vr(1);
        y_gain = 1.0;
        y_offset = 0.0;
        scale.horizontal = 1000;
        scale.vertical = 1000;
        delay = 0;

        % write results.
        [sucess] = cvqm_save_calibration(calibration_file, calibration, ...
            shift, pvr, y_gain, y_offset, scale, delay); 
        if sucess == 0,
            cvqm_error(error_file, 1, 'Cannot open file to write calibration results.');
            return;
        end

    elseif strcmp(calibration,'manual'),
        % load previous results
        [sucess, shift, pvr, y_gain, y_offset,scale, delay] = cvqm_load_calibration(calibration_file);
        if sucess == 0,
            cvqm_error(error_file, 1, 'Cannot open file to read manual calibration.');
            return;
        end

        % apply temporal registration
        dll_lowbw_temporal_original(1, delay);
        dll_lowbw_temporal_processed(2, delay);

        % apply other calibration settings
        dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);

    elseif strcmp(calibration,'frtimemanual'),
        % load previous results
        [sucess, shift, pvr, y_gain, y_offset,scale, delay] = cvqm_load_calibration(calibration_file);
        if sucess == 0,
            cvqm_error(error_file, 1, 'Cannot open file to read manual calibration.');
            return;
        end

        % apply other calibration settings
        dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);
        
        % perform full bandwidth temporal registration
        [delay, sucess, is_still, is_ambiguous, is_uncert] = dll_itu_temporal_registration;
        
        if sucess == 0 && is_still,
            cvqm_error(error_file, 2, 'Still or nearly still sequence; temporal registration cannot be computed.');
        end
        if sucess == 0 && is_ambiguous,
            cvqm_error(error_file, 2, 'Temporal registration results ambiguous.');
        end
        if sucess == 0 && is_uncert,
            cvqm_error(error_file, 2, 'Temporal mis-registration exceeds 1 second uncertainty limit.');
        end

        % write results.
        [sucess] = cvqm_save_calibration(calibration_file, calibration, shift, pvr, y_gain, y_offset, scale, delay); 
        if sucess == 0,
            cvqm_error(error_file, 1, 'Cannot open file to write full reference calibration results.');
            return;
        end

    elseif strcmp(calibration,'rrtimemanual'),
        % load previous results
        [sucess, shift, pvr, y_gain, y_offset,scale, delay] = cvqm_load_calibration(calibration_file);
        if sucess == 0,
            cvqm_error(error_file, 1, 'Cannot open file to read manual calibration.');
            return;
        end

        % apply other calibration settings
        dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);
        
        % perform RR bandwidth temporal registration
        % perform low bandwidth temporal registration
        [ti2_orig, ti10_orig, ymean_orig, orig_is_white_clip, orig_is_black_clip] = ...
            dll_lowbw_temporal_features(1, total_sec, pvr);

        [ti2_proc, ti10_proc, ymean_proc, proc_is_white_clip, proc_is_black_clip] = ...
            dll_lowbw_temporal_features(2, total_sec, pvr);

        if orig_is_white_clip, 
            cvqm_error(error_file, 2, 'White level clipping detected on original video sequence may cause VQM errors.');
        end
        if orig_is_black_clip,
            cvqm_error(error_file, 2, 'Black level clipping detected on original video sequence may cause VQM errors.');
        end
        if proc_is_white_clip,
            cvqm_error(error_file, 2, 'White level clipping detected on processed video sequence may cause VQM errors.');
        end
        if proc_is_black_clip,
            cvqm_error(error_file, 2, 'Black level clipping detected on processed video sequence may cause VQM errors.');
        end

        uncert = 1;
        [delay, sucess, is_still] = ...
            dll_lowbw_temporal (1, ti2_orig, ti2_proc, ti10_orig, ti10_proc, ymean_orig, ymean_proc, uncert);
        if sucess == 0 && is_still,
            cvqm_error(error_file, 2, 'Still or nearly still sequence; temporal registration cannot be computed.');
        elseif sucess == 0,
            cvqm_error(error_file, 2, 'Lowbw temporal registration algorithm failed.');
        end

        dll_lowbw_temporal_original(1, delay);
        dll_lowbw_temporal_processed(2, delay);

        if abs(delay) >= round(fps1),
            cvqm_error(error_file, 2, 'Temporal mis-registration exceeds 1 second uncertainty limit.');
        end

        % write results.
        [sucess] = cvqm_save_calibration(calibration_file, calibration, shift, pvr, y_gain, y_offset, scale, delay); 
        if sucess == 0,
            cvqm_error(error_file, 1, 'Cannot open file to write full reference calibration results.');
            return;
        end

    elseif strcmp(calibration,'rrtime'),

        warning off;
        delete(calibration_file);
        warning on;

        % Use default "no calibration" values for everything except delay
        shift.horizontal = 0;
        shift.vertical = 0;
        pvr = dll_default_vr(1);
        y_gain = 1.0;
        y_offset = 0.0;
        scale.horizontal = 1000;
        scale.vertical = 1000;

        % perform low bandwidth temporal registration
        [ti2_orig, ti10_orig, ymean_orig, orig_is_white_clip, orig_is_black_clip] = ...
            dll_lowbw_temporal_features(1, total_sec, pvr);

        [ti2_proc, ti10_proc, ymean_proc, proc_is_white_clip, proc_is_black_clip] = ...
            dll_lowbw_temporal_features(2, total_sec, pvr);

        if orig_is_white_clip, 
            cvqm_error(error_file, 2, 'White level clipping detected on original video sequence may cause VQM errors.');
        end
        if orig_is_black_clip,
            cvqm_error(error_file, 2, 'Black level clipping detected on original video sequence may cause VQM errors.');
        end
        if proc_is_white_clip,
            cvqm_error(error_file, 2, 'White level clipping detected on processed video sequence may cause VQM errors.');
        end
        if proc_is_black_clip,
            cvqm_error(error_file, 2, 'Black level clipping detected on processed video sequence may cause VQM errors.');
        end

        uncert = 1;
        [delay, sucess, is_still] = ...
            dll_lowbw_temporal (1, ti2_orig, ti2_proc, ti10_orig, ti10_proc, ymean_orig, ymean_proc, uncert);
        if sucess == 0 && is_still,
            cvqm_error(error_file, 2, 'Still or nearly still sequence; temporal registration cannot be computed.');
        elseif sucess == 0,
            cvqm_error(error_file, 2, 'Lowbw temporal registration algorithm failed.');
        end

        dll_lowbw_temporal_original(1, delay);
        dll_lowbw_temporal_processed(2, delay);

        if abs(delay) >= round(fps1),
            cvqm_error(error_file, 2, 'Temporal mis-registration exceeds 1 second uncertainty limit.');
        end

        % perform low bandwidth valid region
        [ovr] = dll_orig_valid_region;
        [pvr] = dll_proc_valid_region(ovr);
        
        if (pvr.bottom - pvr.top + 1) / rows1 < 0.55 || ...
                (pvr.right - pvr.left + 1) / cols1 < 0.80,
            cvqm_error(error_file, 2, 'Greatly reduced valid region detected.');
        end

        % write results.
        [sucess] = cvqm_save_calibration(calibration_file, calibration, shift, pvr, y_gain, y_offset, scale, delay); 
        if sucess == 0,
            cvqm_error(error_file, 1, 'Cannot open file to write lowbw calibration results.');
            return;
        end


    elseif strcmp(calibration,'frtime'),

        warning off;
        delete(calibration_file);
        warning on;

        % Use default "no calibration" values for everything except delay
        shift.horizontal = 0;
        shift.vertical = 0;
        pvr = dll_default_vr(1);
        y_gain = 1.0;
        y_offset = 0.0;
        scale.horizontal = 1000;
        scale.vertical = 1000;

        % perform full bandwidth temporal registration
        [delay, sucess, is_still, is_ambiguous, is_uncert] = dll_itu_temporal_registration;
        
        if sucess == 0 && is_still,
            cvqm_error(error_file, 2, 'Still or nearly still sequence; temporal registration cannot be computed.');
        end
        if sucess == 0 && is_ambiguous,
            cvqm_error(error_file, 2, 'Temporal registration results ambiguous.');
        end
        if sucess == 0 && is_uncert,
            cvqm_error(error_file, 2, 'Temporal mis-registration exceeds 1 second uncertainty limit.');
        end

        % perform full bandwidth valid region
        [ovr] = dll_itu_orig_valid_region;
        [pvr] = dll_itu_proc_valid_region(ovr);
        
        if (pvr.bottom - pvr.top + 1) / rows1 < 0.55 || ...
                (pvr.right - pvr.left + 1) / cols1 < 0.80,
            cvqm_error(error_file, 2, 'Greatly reduced valid region detected.');
        end

        % input improved PVR estimate 
        dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);

        % write results.
        [sucess] = cvqm_save_calibration(calibration_file, calibration, shift, pvr, y_gain, y_offset, scale, delay); 
        if sucess == 0,
            cvqm_error(error_file, 1, 'Cannot open file to write full reference calibration results.');
            return;
        end

    elseif strcmp(calibration,'rrcal') || strcmp(calibration,'rrcalscale') || strcmp(calibration,'rrcal2') || strcmp(calibration,'rrcal2scale'),
        warning off;
        delete(calibration_file);
        warning on;

        % note whether or not to estimate scaling.
        if strcmp(calibration,'rrcal2'),
            no_scaling = 1;
            is_rrcal2 = 1;
        elseif strcmp(calibration,'rrcal'),
            no_scaling = 1;
            is_rrcal2 = 0;
        elseif strcmp(calibration,'rrcal2scale'),
            no_scaling = 0;
            is_rrcal2 = 1;
        elseif strcmp(calibration,'rrcalscale'),
            no_scaling = 0;
            is_rrcal2 = 0;
        end

        total_sec = floor(total_sec);

        % default low bandwidth valid region
        [vr] = dll_default_vr(1);

        % perform low bandwidth temporal registration
        [ti2_orig, ti10_orig, ymean_orig, orig_is_white_clip, orig_is_black_clip] = ...
            dll_lowbw_temporal_features(1, total_sec, vr);

        [ti2_proc, ti10_proc, ymean_proc, proc_is_white_clip, proc_is_black_clip] = ...
            dll_lowbw_temporal_features(2, total_sec, vr);

        % note white & black level clipping, if any
        if orig_is_white_clip,
            cvqm_error(error_file, 2, 'White level clipping detected on original video sequence may cause VQM errors.');
        end
        if orig_is_black_clip,
            cvqm_error(error_file, 2, 'Black level clipping detected on original video sequence may cause VQM errors.');
        end
        if proc_is_white_clip,
            cvqm_error(error_file, 2, 'White level clipping detected on processed video sequence may cause VQM errors.');
        end
        if proc_is_black_clip,
            cvqm_error(error_file, 2, 'Black level clipping detected on processed video sequence may cause VQM errors.');
        end
        
        % example code to use quantizer. 
%         [ti2_index, ti10_index, y_index] =
%           dll_lowbw_temporal_quant(1, ti2_orig, ti10_orig, ymean_orig);
%         [ti2_orig, ti10_orig, ymean_orig] =
%           dll_lowbw_temporal_quant(0, ti2_index, ti10_index, y_index);

        % complete temporal registration.
        uncert = 1;
        [delay, sucess, is_still] = ...
            dll_lowbw_temporal (1, ti2_orig, ti2_proc, ti10_orig, ti10_proc, ymean_orig, ...
            ymean_proc, uncert, 'field');
        if sucess == 0 && is_still,
            cvqm_error(error_file, 2, 'Still or nearly still sequence; initial temporal registration cannot be computed.');
        elseif sucess == 0,
            cvqm_error(error_file, 2, 'Initial temporal registration algorithm failed.');
        end
        
        if mod(delay,1) == 0.5,
            dll_calib_video('set_reframe',1);
            delay = delay - 0.5;
        end

        % re-calculate seconds of video to be used (may have changed)
        [total_sec] = min( dll_calib_video('total_sec',1), dll_calib_video('total_sec',2) );
        if total_sec > 15,
            total_sec = 15;
        end
        
        dll_lowbw_temporal_original(1, delay);
        dll_lowbw_temporal_processed(2, delay);
        
        [total_sec] = min( dll_calib_video('total_sec',1), dll_calib_video('total_sec',2) );
        if total_sec > 15,
            total_sec = 15;
        end

        % calculate shift, & scaling 
        [seed_state] = dll_lowbw_calib_initialize;
        [orig_pixels, orig_horiz_profile, orig_vert_profile] = dll_lowbw_calib_original(1, seed_state, total_sec);

        % to quantize:
        %   [index_horiz_profile, index_vert_profile] =
        %       dll_lowbw_calib_quant(1, orig_horiz_profile, orig_vert_profile);
        %   [orig_horiz_profile, orig_vert_profile] =
        %       dll_lowbw_calib_quant(0, index_horiz_profile, index_vert_profile);

        [shift, scale,status] = ...
            dll_lowbw_calib_processed(2, seed_state, total_sec, orig_pixels, ...
            orig_horiz_profile, orig_vert_profile, no_scaling);
        if status.scale && ~no_scaling,
            cvqm_error(error_file, 2, 'Actual spatial scaling may be beyond search limits.');
        end
        if status.shift,
            cvqm_error(error_file, 2, 'Actual spatial shift may be beyond search limits.');
        end
        if (strcmp(video_standard,'interlace_lower_field_first') || ...
                strcmp(video_standard,'interlace_upper_field_first')) && mod(shift.vertical,2),
            cvqm_error(error_file, 2, 'Processed video will be reframed. Actual temporal delay is 0.5 frames greater than reported value.');
        end
        if abs(shift.horizontal) > 8 || abs(shift.vertical) > 5,
            cvqm_error(error_file, 2, 'Extreme spatial shift detected.');
        end
        if scale.vertical ~= 1000 || scale.horizontal ~= 1000,
            cvqm_error(error_file, 2, 'Video scaling detected; please examine other scenes.');
        end


        % input above estimates 
        y_gain = 1.0;
        y_offset = 0.0;
        dll_calib_video('calibration', shift.horizontal, shift.vertical, vr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);

        % re-calculate seconds of video to be used (may have changed)
        [total_sec] = min( dll_calib_video('total_sec',1), dll_calib_video('total_sec',2) );
        if total_sec > 15,
            total_sec = 15;
        end
        
        % perform low bandwidth valid region
        [ovr] = dll_orig_valid_region;
        [pvr] = dll_proc_valid_region(ovr);
        
        if (pvr.bottom - pvr.top + 1) / rows1 < 0.55 || ...
                (pvr.right - pvr.left + 1) / cols1 < 0.80,
            cvqm_error(error_file, 2, 'Greatly reduced valid region detected.');
        end

        % input improved PVR estimate 
        dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);
        
        % estimate luminance gain & offset
        if is_rrcal2,
            [orig_y, orig_cb, orig_cr, yesno] = dll_lowbw_gain_v2_original(1, total_sec);

            % How to quantize
            % [index_y, index_cb, index_cr] = dll_lowbw_gain_v2_quant(1, orig_y, orig_cb, orig_cr); % quantize
            % [orig_y, orig_cb, orig_cr] = dll_lowbw_gain_v2_quant(0, index_y, index_cb, index_cr); % reconstruct
            
            [y_gain, y_offset, cb_gain, cb_offset, cr_gain, cr_offset, sucess] = ...
                dll_lowbw_gain_v2_processed(2, total_sec, orig_y, orig_cb, orig_cr, yesno);
            cvqm_error(error_file, 2, 'Color Gain & Offset estimated but not removed.');
        else
            [orig_blocks] = dll_lowbw_gain_original(1, total_sec);
            [y_gain, y_offset, sucess] = dll_lowbw_gain_processed(2, total_sec, orig_blocks);
        end
        
        if y_gain < 0.9 || y_gain > 1.1,
            cvqm_error(error_file, 2, 'Extreme Luminance Gain detected.');
        end
        if y_offset < -20 || y_offset > 20,
            cvqm_error(error_file, 2, 'Extreme Luminance Offset detected.');
        end
        if sucess == 0,
            cvqm_error(error_file, 2, 'Warning: algorithm used to estimate Luminance Gain & Offset may have failed.');
        end
        if sucess == -1,
            cvqm_error(error_file, 2, 'Luminance gain & offset algorithm detected extreme values or failed; results discarded.');
        end
        if is_rrcal2 && isnan(cb_gain),
            cvqm_error(error_file, 2, 'Warning: algorithm used to estimate Cb Gain & Offset failed.');
        end
        if is_rrcal2 && isnan(cr_gain),
            cvqm_error(error_file, 2, 'Warning: algorithm used to estimate Cr Gain & Offset failed.');
        end

        % input luminance gain & offset estimates
        dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);

        % perform low bandwidth temporal registration a second time!
        [ti2_orig, ti10_orig, ymean_orig, orig_is_white_clip, orig_is_black_clip] = ...
            dll_lowbw_temporal_features(1, total_sec, pvr);

        [ti2_proc, ti10_proc, ymean_proc, proc_is_white_clip, orig_is_black_clip] = ...
            dll_lowbw_temporal_features(2, total_sec, pvr);

        uncert = 1;
        [delay2, sucess, is_still] = ...
            dll_lowbw_temporal (1, ti2_orig, ti2_proc, ti10_orig, ti10_proc, ymean_orig, ymean_proc, uncert);

        if sucess == 0 && is_still,
            cvqm_error(error_file, 2, 'Still or nearly still sequence; temporal registration cannot be computed.');
        elseif sucess == 0,
            cvqm_error(error_file, 2, 'Final temporal registration algorithm failed.');
        end

        dll_lowbw_temporal_original(1, delay2);
        dll_lowbw_temporal_processed(2, delay2);

        if abs(delay+delay2) >= round(fps1),
            cvqm_error(error_file, 2, 'Temporal mis-registration exceeds 1 second uncertainty limit.');
        end

        % write results.
        if is_rrcal2,
            [sucess] = cvqm_save_calibration(calibration_file, calibration, shift, pvr, y_gain, y_offset, scale, delay+delay2, cb_gain, cb_offset, cr_gain, cr_offset); 
        else
            [sucess] = cvqm_save_calibration(calibration_file, calibration, shift, pvr, y_gain, y_offset, scale, delay+delay2); 
        end
        if sucess == 0,
            cvqm_error(error_file, 1, 'Cannot open file to write lowbw calibration results.');
            return;
        end

    elseif strcmp(calibration,'frcal'),
        warning off;
        delete(calibration_file);
        warning on;

        % Use default "no calibration" values for scaling
        shift.horizontal = 0;
        shift.vertical = 0;
        pvr = dll_default_vr(1);
        y_gain = 1.0;
        y_offset = 0.0;
        
        scale.horizontal = 1000;
        scale.vertical = 1000;
        
        % perform spatial registration
        [shift.horizontal, shift.vertical, status] = dll_itu_spatial; 
        if status == 0,
            cvqm_error(error_file, 2, 'Spatial Registration algorithm failed.  Assume no shift.');
            shift.horizontal = 0;
            shift.vertical = 0;
        end
        if (strcmp(video_standard,'interlace_lower_field_first') || ...
                strcmp(video_standard,'interlace_upper_field_first')) && mod(shift.vertical,2),
            cvqm_error(error_file, 2, 'Processed video will be reframed.  Actual temporal delay is 0.5 frames greater than reported value.');
        end
        if abs(shift.horizontal) > 8 || abs(shift.vertical) > 5,
            cvqm_error(error_file, 2, 'Extreme spatial shift detected.');
        end

        % Apply spatial registration settings
        dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);

        % perform full bandwidth valid region
        [ovr] = dll_itu_orig_valid_region;
        [pvr] = dll_itu_proc_valid_region(ovr);
        
        if (pvr.bottom - pvr.top + 1) / rows1 < 0.55 || ...
                (pvr.right - pvr.left + 1) / cols1 < 0.80,
            cvqm_error(error_file, 2, 'Greatly reduced valid region detected.');
        end
        
        % Apply Valid Region settings
        dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);

        % gain/offset
        [y_gain, y_offset, status] = dll_itu_gain_offset;
        if y_gain < 0.9 || y_gain > 1.1,
            cvqm_error(error_file, 2, 'Extreme luminance gain detected.');
        end
        if y_offset < -20 || y_offset > 20,
            cvqm_error(error_file, 2, 'Extreme luminance offset detected.');
        end
        if status == 0,
            cvqm_error(error_file, 2, 'Luminance gain & offset algorithm detected extreme values or failed; results discarded.');
        end

        % Apply Luminance Gain & Offset settings
        dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);

        % perform full bandwidth temporal registration
        [delay, sucess, is_still, is_ambiguous, is_uncert] = dll_itu_temporal_registration;
        
        if sucess == 0 && is_still,
            cvqm_error(error_file, 2, 'Still or nearly still sequence; temporal registration cannot be computed.');
        end
        if sucess == 0 && is_ambiguous,
            cvqm_error(error_file, 2, 'Temporal registration results ambiguous.');
        end
        if sucess == 0 && is_uncert,
            cvqm_error(error_file, 2, 'Temporal mis-registration exceeds 1 second uncertainty limit.');
        end

        % write results.
        [sucess] = cvqm_save_calibration(calibration_file, calibration, shift, pvr, y_gain, y_offset, scale, delay); 
        if sucess == 0,
            cvqm_error(error_file, 1, 'Cannot open file to write full reference calibration results.');
            return;
        end
        
    elseif strcmpi(calibration, 'psnr_search')
        warning off;
        delete(calibration_file);
        warning on;

        % set video standard
        if strcmp(dll_video('get_video_standard', 1), 'interlace_lower_field_first') || ...
                strcmp(dll_video('get_video_standard', 1), 'interlace_upper_field_first');
            dll_calib_video('image_mode', 'field');
        elseif strcmp(dll_video('get_video_standard', 1), 'progressive')
            dll_calib_video('image_mode', 'frame');
        else
            error('Incorrect video standard type.');
        end
        
        % Use default "no calibration" values for scaling
        shift.horizontal = 0;
        shift.vertical = 0;
        y_gain = 1.0;
        y_offset = 0.0;
        scale.horizontal = 1000;
        scale.vertical = 1000;
        
        [ovr] = dll_orig_valid_region;
        [pvr] = dll_proc_valid_region(ovr);  % calculate pvr
        
        dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);
        
        % set 'spatial_uncertainty' to 3
        su = 3;
        
        % set 'temporal_uncertainty' to number of frames in 1 sec, round up
        tu = ceil(dll_video('fps', 1));
        
        % call function dll_psnr_search.m
        results = dll_psnr_search({'spatial_uncertainty', su, su}, {'temporal_uncertainty', tu}, {'fraction_sampled', fs});
        
        % Set calibration values in dll_calib_video.m with results
        pvr = dll_calib_video('pvr');
        shift.horizontal = -results.xshift;
        shift.vertical = -results.yshift;
        y_gain = 1/results.gain;
        y_offset = -results.offset*y_gain;
        delay = -results.tshift;
        dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);
        
        % crop video based on time delay
        if mod(delay, 1) == 0.5
            delay = delay - 0.5;
        end
        if delay > 0,
            dll_video('discard', 2, delay / dll_video('fps'));
        elseif delay < 0,
            dll_video('discard', 1, -delay / dll_video('fps'));
        end

        % Calculate new pvr and include it in the calibration variable
        [pvr] = dll_proc_valid_region(ovr);
        dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);
        
        % save returned PSNR value in case model selected is also
        % exhaustive search psnr
        savedPSNR = results.psnr;
        
        % write results.
        [sucess] = cvqm_save_calibration(calibration_file, calibration, shift, pvr, y_gain, y_offset, scale, delay); 
        if sucess == 0,
            cvqm_error(error_file, 1, 'Cannot open file to write full reference calibration results.');
            return;
        end
    elseif strcmpi(calibration, 'fast_psnr_search1')
        warning off;
        delete(calibration_file);
        warning on;

        no_scaling = 1;
        
        total_sec = floor(total_sec);

        % default low bandwidth valid region
        [vr] = dll_default_vr(1);

        % perform low bandwidth temporal registration
        [ti2_orig, ti10_orig, ymean_orig, orig_is_white_clip, orig_is_black_clip] = ...
            dll_lowbw_temporal_features(1, total_sec, vr);

        [ti2_proc, ti10_proc, ymean_proc, proc_is_white_clip, proc_is_black_clip] = ...
            dll_lowbw_temporal_features(2, total_sec, vr);

        % note white & black level clipping, if any
        if orig_is_white_clip,
            cvqm_error(error_file, 2, 'White level clipping detected on original video sequence may cause VQM errors.');
        end
        if orig_is_black_clip,
            cvqm_error(error_file, 2, 'Black level clipping detected on original video sequence may cause VQM errors.');
        end
        if proc_is_white_clip,
            cvqm_error(error_file, 2, 'White level clipping detected on processed video sequence may cause VQM errors.');
        end
        if proc_is_black_clip,
            cvqm_error(error_file, 2, 'Black level clipping detected on processed video sequence may cause VQM errors.');
        end

        % complete temporal registration.
        uncert = 1;
        [delay, sucess, is_still] = ...
            dll_lowbw_temporal (1, ti2_orig, ti2_proc, ti10_orig, ti10_proc, ymean_orig, ...
            ymean_proc, uncert, 'field');
        if sucess == 0 && is_still,
            cvqm_error(error_file, 2, 'Still or nearly still sequence; initial temporal registration cannot be computed.');
        elseif sucess == 0,
            cvqm_error(error_file, 2, 'Initial temporal registration algorithm failed.');
        end
        
        if mod(delay,1) == 0.5,
            dll_calib_video('set_reframe',1);
            delay = delay - 0.5;
        end

        % re-calculate seconds of video to be used (may have changed)
        [total_sec] = min( dll_calib_video('total_sec',1), dll_calib_video('total_sec',2) );
        if total_sec > 15,
            total_sec = 15;
        end
        
        dll_lowbw_temporal_original(1, delay);
        dll_lowbw_temporal_processed(2, delay);
        
        [total_sec] = min( dll_calib_video('total_sec',1), dll_calib_video('total_sec',2) );
        if total_sec > 15,
            total_sec = 15;
        end

        % calculate shift, & scaling 
        [seed_state] = dll_lowbw_calib_initialize;
        [orig_pixels, orig_horiz_profile, orig_vert_profile] = dll_lowbw_calib_original(1, seed_state, total_sec);


        [shift, scale,status] = ...
            dll_lowbw_calib_processed(2, seed_state, total_sec, orig_pixels, ...
            orig_horiz_profile, orig_vert_profile, no_scaling);
        if status.scale && ~no_scaling,
            cvqm_error(error_file, 2, 'Actual spatial scaling may be beyond search limits.');
        end
        if status.shift,
            cvqm_error(error_file, 2, 'Actual spatial shift may be beyond search limits.');
        end
        if (strcmp(video_standard,'interlace_lower_field_first') || ...
                strcmp(video_standard,'interlace_upper_field_first')) && mod(shift.vertical,2),
            cvqm_error(error_file, 2, 'Processed video will be reframed. Actual temporal delay is 0.5 frames greater than reported value.');
        end
        if abs(shift.horizontal) > 8 || abs(shift.vertical) > 5,
            cvqm_error(error_file, 2, 'Extreme spatial shift detected.');
        end
        if scale.vertical ~= 1000 || scale.horizontal ~= 1000,
            cvqm_error(error_file, 2, 'Video scaling detected; please examine other scenes.');
        end

        prevshift = shift;
        
        % input above estimates 
        y_gain = 1.0;
        y_offset = 0.0;
        dll_calib_video('calibration', shift.horizontal, shift.vertical, vr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);

        % re-calculate seconds of video to be used (may have changed)
        [total_sec] = min( dll_calib_video('total_sec',1), dll_calib_video('total_sec',2) );
        if total_sec > 15,
            total_sec = 15;
        end
        
        % perform low bandwidth valid region
        [ovr] = dll_orig_valid_region;
        [pvr] = dll_proc_valid_region(ovr);
        
        if (pvr.bottom - pvr.top + 1) / rows1 < 0.55 || ...
                (pvr.right - pvr.left + 1) / cols1 < 0.80,
            cvqm_error(error_file, 2, 'Greatly reduced valid region detected.');
        end

        % input improved PVR estimate 
        dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);
        

        % perform low bandwidth temporal registration a second time!
        [ti2_orig, ti10_orig, ymean_orig, orig_is_white_clip, orig_is_black_clip] = ...
            dll_lowbw_temporal_features(1, total_sec, pvr);

        [ti2_proc, ti10_proc, ymean_proc, proc_is_white_clip, orig_is_black_clip] = ...
            dll_lowbw_temporal_features(2, total_sec, pvr);

        uncert = 1;
        [delay2, sucess, is_still] = ...
            dll_lowbw_temporal (1, ti2_orig, ti2_proc, ti10_orig, ti10_proc, ymean_orig, ymean_proc, uncert);

        if sucess == 0 && is_still,
            cvqm_error(error_file, 2, 'Still or nearly still sequence; temporal registration cannot be computed.');
        elseif sucess == 0,
            cvqm_error(error_file, 2, 'Final temporal registration algorithm failed.');
        end

        dll_lowbw_temporal_original(1, delay2);
        dll_lowbw_temporal_processed(2, delay2);

        if abs(delay+delay2) >= round(fps1),
            cvqm_error(error_file, 2, 'Temporal mis-registration exceeds 1 second uncertainty limit.');
        end

        prevdelay = delay + delay2;
        
        % set video standard
        if strcmp(dll_video('get_video_standard', 1), 'interlace_lower_field_first') || ...
                strcmp(dll_video('get_video_standard', 1), 'interlace_upper_field_first');
            dll_calib_video('image_mode', 'field');
        elseif strcmp(dll_video('get_video_standard', 1), 'progressive')
            dll_calib_video('image_mode', 'frame');
        else
            error('Incorrect video standard type.');
        end

        
        % Set 'spatial_uncertainty' to 1
        su = 1;
        
        % Set 'temporal_uncertainty'
        tu = ceil(dll_video('fps', 1)*.5);
        
        % call function dll_psnr_search.m
        results = dll_psnr_search({'spatial_uncertainty', su, su}, {'temporal_uncertainty', tu}, {'fraction_sampled', fs});
        
        % Set calibration values in dll_calib_video.m with results
        pvr = dll_calib_video('pvr');
        shift.horizontal = prevshift.horizontal -results.xshift;
        shift.vertical = prevshift.vertical - results.yshift;
        y_gain = 1/results.gain;
        y_offset = -results.offset*y_gain;
        delay_temp = -results.tshift;
        dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);
        
        % crop video based on the time delay
        if mod(delay_temp, 1) == 0.5
            delay_temp = delay_temp - 0.5;
        end
        if delay > 0,
            dll_video('discard', 2, delay_temp / dll_video('fps'));
        elseif delay < 0,
            dll_video('discard', 1, -delay_temp / dll_video('fps'));
        end
        delay = delay_temp + prevdelay;

        % Calculate new pvr and include it in the calibration variable
        [pvr] = dll_proc_valid_region(ovr);
        dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);
        
        % save returned PSNR value in case model selected is also
        % exhaustive search psnr
        savedPSNR = results.psnr;
        
        % write results.
        [sucess] = cvqm_save_calibration(calibration_file, calibration, shift, pvr, y_gain, y_offset, scale, delay); 
        if sucess == 0,
            cvqm_error(error_file, 1, 'Cannot open file to write full reference calibration results.');
            return;
        end

        
    elseif strcmpi(calibration, 'fast_psnr_search2')
                warning off;
        delete(calibration_file);
        warning on;

        % Use default "no calibration" values for scaling
        shift.horizontal = 0;
        shift.vertical = 0;
        pvr = dll_default_vr(1);
        y_gain = 1.0;
        y_offset = 0.0;
        
        scale.horizontal = 1000;
        scale.vertical = 1000;
        
        % perform spatial registration
        [shift.horizontal, shift.vertical, status] = dll_itu_spatial; 
        if status == 0,
            cvqm_error(error_file, 2, 'Spatial Registration algorithm failed.  Assume no shift.');
            shift.horizontal = 0;
            shift.vertical = 0;
        end
        if (strcmp(video_standard,'interlace_lower_field_first') || ...
                strcmp(video_standard,'interlace_upper_field_first')) && mod(shift.vertical,2),
            cvqm_error(error_file, 2, 'Processed video will be reframed.  Actual temporal delay is 0.5 frames greater than reported value.');
        end
        if abs(shift.horizontal) > 8 || abs(shift.vertical) > 5,
            cvqm_error(error_file, 2, 'Extreme spatial shift detected.');
        end

        prevshift = shift;
        
        % Apply spatial registration settings
        dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);

        % perform full bandwidth valid region
        [ovr] = dll_itu_orig_valid_region;
        [pvr] = dll_itu_proc_valid_region(ovr);
        
        if (pvr.bottom - pvr.top + 1) / rows1 < 0.55 || ...
                (pvr.right - pvr.left + 1) / cols1 < 0.80,
            cvqm_error(error_file, 2, 'Greatly reduced valid region detected.');
        end
        
        % Apply Valid Region settings
        dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);

        % perform full bandwidth temporal registration
        [delay, sucess, is_still, is_ambiguous, is_uncert] = dll_itu_temporal_registration;
        
        if sucess == 0 && is_still,
            cvqm_error(error_file, 2, 'Still or nearly still sequence; temporal registration cannot be computed.');
        end
        if sucess == 0 && is_ambiguous,
            cvqm_error(error_file, 2, 'Temporal registration results ambiguous.');
        end
        if sucess == 0 && is_uncert,
            cvqm_error(error_file, 2, 'Temporal mis-registration exceeds 1 second uncertainty limit.');
        end
        
        prevdelay = delay;
        
        % set video standard
        if strcmp(dll_video('get_video_standard', 1), 'interlace_lower_field_first') || ...
                strcmp(dll_video('get_video_standard', 1), 'interlace_upper_field_first');
            dll_calib_video('image_mode', 'field');
        elseif strcmp(dll_video('get_video_standard', 1), 'progressive')
            dll_calib_video('image_mode', 'frame');
        else
            error('Incorrect video standard type.');
        end

        
        % Set 'spatial_uncertainty' to 1
        su = 1;
        
        % Set 'temporal_uncertainty'
        tu = ceil(dll_video('fps', 1)*.5);
        
        % call function dll_psnr_search.m
        results = dll_psnr_search({'spatial_uncertainty', su, su}, {'temporal_uncertainty', tu}, {'fraction_sampled', fs});
        
        % Set calibration values in dll_calib_video.m with results
        pvr = dll_calib_video('pvr');
        shift.horizontal = prevshift.horizontal -results.xshift;
        shift.vertical = prevshift.vertical - results.yshift;
        y_gain = 1/results.gain;
        y_offset = -results.offset*y_gain;
        delay_temp = -results.tshift;
        dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);
        
        % crop video based on the time delay
        if mod(delay_temp, 1) == 0.5
            delay_temp = delay_temp - 0.5;
        end
        if delay > 0,
            dll_video('discard', 2, delay_temp / dll_video('fps'));
        elseif delay < 0,
            dll_video('discard', 1, -delay_temp / dll_video('fps'));
        end
        delay = delay_temp + prevdelay;

        % Calculate new pvr and include it in the calibration variable
        [pvr] = dll_proc_valid_region(ovr);
        dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
            y_gain, y_offset, scale.horizontal, scale.vertical);
        
        % save returned PSNR value in case model selected is also
        % exhaustive search psnr
        savedPSNR = results.psnr;
        
        % write results.
        [sucess] = cvqm_save_calibration(calibration_file, calibration, shift, pvr, y_gain, y_offset, scale, delay); 
        if sucess == 0,
            cvqm_error(error_file, 1, 'Cannot open file to write full reference calibration results.');
            return;
        end

        
        
    end
    
    % performs variable frame delay calibration for any model option that
    % includes vfd calibration
    if strcmp(model, 'psnr_vfd') || strcmp(model, 'vqm_vfd')
        if strcmp(calibration, 'none')
            dll_calib_video('max_roi');
        end
        tu = ceil(dll_video('fps', 1));
        
        % VFD calibration
        if strcmp(dll_video('get_video_standard', 1), 'progressive')
            [results results_fuzzy] = dll_vfd(delay, 't_uncert', tu, 'causal'); % 8/5/11 added 'delay' as an input
        else
            [results results_fuzzy] = dll_vfd(delay, 't_uncert', tu, 'causal', 'reframe'); % 8/5/11 added 'delay' as an input
        end
       
        dll_calib_video('set_vfd', results, results_fuzzy);
        
        if ~exist('delay', 'var')
            delay = 0;
        end
        
        % Write out the results of the vfd calibration to a .txt file
        vfd_file_name = sprintf('%s_vfd.csv', processed_file);
        dll_vfd_print(results, vfd_file_name, delay);
        
        % If the user selected the psnr_vfd model, the new luminance
        % gain/offset calculated can greatly affect the PSNR value.  This
        % code is to update those values and update the calibration file.
        if strcmp(model, 'psnr_vfd')
            dll_video('set_rewind', 1);
            dll_video('set_rewind', 2);

            tlength_original = dll_calib_video('total_sec', 1);
            tlength_processed = dll_calib_video('total_sec', 2);
            if delay < 0   % 8/5/11 if delay is -, take frames off the end of the processed
                tlength_processed = tlength_processed + delay/dll_video('fps', 2);
            end
            dll_video('set_tslice', 1, tlength_original);
            dll_video('set_tslice', 2, tlength_processed);
            y_orig = dll_calib_video('tslice', 1);
            y_proc = dll_calib_video('tslice', 2);
            
            [nrows, ncols, nsamps] = size(y_proc);
            y_proc = reshape(y_proc,nrows*ncols*nsamps,1);
            y_orig = reshape(y_orig,nrows*ncols*nsamps,1);
            rand_nums = round(randperm((nrows*ncols*nsamps))); %Randomizes numbers from 1 to nrows*ncols*nframes
            this_fit = polyfit(y_proc(rand_nums(1:round(nrows*ncols*nsamps*fs))),...
                                    y_orig(rand_nums(1:round(nrows*ncols*nsamps*fs))),1);
            peak = 255.0;
            % Calculate the new PSNR value
            savedPSNR = 10*(log10(peak*peak)-log10(sum(sum(sum(((this_fit(1)*y_proc+this_fit(2))-y_orig).^2)))/(nrows*ncols*nsamps)));
            clear y_orig y_proc;
            y_gain = y_gain/this_fit(1);
            y_offset = y_offset - this_fit(2);
            
            dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
                y_gain, y_offset, scale.horizontal, scale.vertical);
            
            % Now that we have the new y_gain and y_offset, write the
            % calibration file to reflect the updated values.
            warning off;
            delete(calibration_file);
            warning on;
            [sucess] = cvqm_save_calibration(calibration_file, calibration, shift, pvr, y_gain, y_offset, scale, delay); 
            if sucess == 0,
                cvqm_error(error_file, 1, 'Cannot open file to write full reference calibration results.');
                return;
            end

            dll_video('rewind', 1);
            dll_video('rewind', 2);
        end
    end
    
    % re-caluclate seconds of video to be used (may have changed)
    [total_sec] = min( dll_calib_video('total_sec',1), dll_calib_video('total_sec',2) );
    if total_sec > 15,
        total_sec = 15;
    end

    % calculate model
    if strcmp(model,'none'),
        warning off;
        delete(model_file);
        warning on;

        % do nothing
        sucess = 1;

    elseif strcmp(model,'general'),
        warning off;
        delete(model_file);
        warning on;

        orig_features = dll_features('General', 1, total_sec);
        proc_features = dll_features('General', 2, total_sec);
        [vqm, pars, par_names] = dll_model('vqm',orig_features, proc_features);
        [sucess] = cvqm_save_model(model_file, model, vqm, pars, par_names);

    elseif strcmp(model,'developers'),
        warning off;
        delete(model_file);
        warning on;

        orig_features = dll_features('Developers', 1, total_sec);
        proc_features = dll_features('Developers', 2, total_sec);
        [vqm, pars, par_names] = dll_model('vqm',orig_features, proc_features);
        [sucess] = cvqm_save_model(model_file, model, vqm, pars, par_names);

    elseif strcmp(model,'lowbw'),
        warning off;
        delete(model_file);
        warning on;

        dll_features('Low', 1, total_sec, temporary_file);
        proc_features = dll_features('Low', 2, total_sec);
        [vqm, pars, par_names] = dll_model('vqm',temporary_file, proc_features);
        delete(temporary_file);
        [sucess] = cvqm_save_model(model_file, model, vqm, pars, par_names);

    elseif strcmp(model,'fastlowbw'),
        warning off;
        delete(model_file);
        warning on;

        dll_features('Fast', 1, total_sec, temporary_file);
        proc_features = dll_features('Fast', 2, total_sec);
        [vqm, pars, par_names] = dll_model('vqm',temporary_file, proc_features);
        delete(temporary_file);
        [sucess] = cvqm_save_model(model_file, model, vqm, pars, par_names);
        
    elseif strcmpi(model,'psnr')
        warning off;
        delete(model_file);
        warning on;
        
        % If the user did not select anything that already calculated the
        % PSNR value for the clip pair...
        if ~(strcmpi(calibration, 'psnr_search') || strcmpi(calibration, 'fast_psnr_search1') || strcmpi(calibration, 'fast_psnr_search2'));
            if strcmp(calibration, 'none')
                pvr_temp = dll_calib_video('pvr');
                dll_calib_video('sroi', pvr_temp, 0);
            end
            tslice_orig = dll_calib_video('total_sec', 1);
            tslice_proc = dll_calib_video('total_sec', 2);
            if tslice_orig > tslice_proc
                tslice_orig = tslice_proc;
            elseif tslice_proc > tslice_orig
                tslice_proc = tslice_orig;
            end
            dll_video('set_tslice', 1, tslice_orig);
            dll_video('set_tslice', 2, tslice_proc);
            y_orig = dll_calib_video('tslice', 1);
            y_proc = dll_calib_video('tslice', 2);
            [nrows, ncols, nsamps] = size(y_proc);
            peak = 255.0;
            savedPSNR = 10*(log10(peak*peak)-log10(sum(sum(sum((y_proc-y_orig).^2)))/(nrows*ncols*nsamps))); 
            clear y_orig y_proc;
        end
        
        [sucess] = cvqm_save_model(model_file, model, savedPSNR, {}, {});
        
    elseif strcmp(model,'psnr_vfd')
        warning off;
        delete(model_file);
        warning on;
        
        % Use previously calculated PSNR value
        [success] = cvqm_save_model(model_file, model, savedPSNR, {}, {});
        
    elseif strcmp(model, 'vqm_vfd')
        warning off;
        delete(model_file);
        warning on;
        
        % Collect the calibration data to pass into dll_vqm_vfd
        luminance.gain = y_gain;
        luminance.offset = y_offset;
        nn_model = dll_vqm_vfd(shift, luminance, scale, viewing_distance, delay); % 8/5/11 added delay as an input
        
        [sucess] = cvqm_save_model(model_file, model, nn_model, {}, {});
        
    end
        

    if sucess == 0,
        cvqm_error(error_file, 1, 'Cannot open file to write model results.');
        return;
    end
    
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
    fprintf(fid, '%.15f %s\r\n', pars(cnt), par_names{cnt});
end

fclose(fid);




