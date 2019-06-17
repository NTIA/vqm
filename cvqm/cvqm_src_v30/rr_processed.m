function rr_processed (file_list, video_standard, model, results_log, flag)
% wrapper for RRNR-TV test: interface required for compilation
% call with no arguments for help information.

if nargin == 0,
    fprintf('RR_PROCESSED Version 1.2\n');
    fprintf('  Take processed test sequences in uncompressed big-YUV file format.\n');
    fprintf('  Calculate NTIA low-bandwidth model (or NTIA fast low-bandwidth model) \n');
    fprintf('  and calibration.\n');
    fprintf('SYNTAX\n');
    fprintf('  rr_processed file_list video_standard model results_log\n');
    fprintf('DESCRIPTION\n');
    fprintf('  ''file_list'' is a text file containing original and processed file names in\n');
    fprintf('              pairs, one pair on each line.  Paths are okay.  Each file must be\n');
    fprintf('              in big-YUV format.  After the second file name, optional (i.e.,\n');
    fprintf('              manual) calibration values may be listed. See "EXAMPLE LIST".\n');
    fprintf('         -->  When an original-processed file pair list includes calibration\n');
    fprintf('              values, then automated calibration will be skipped, and manual\n');
    fprintf('              values used instead.\n');
    fprintf('              Optional calibration values are listed in the following order:\n');
    fprintf('       luma_gain luma_offset horiz_shift vert_shift delay\n'); 
    fprintf('               luma_gain is luminance gain, double precision\n');
    fprintf('               luma_offset is luminance offset, double precision\n');
    fprintf('               horiz_shift is horizontal shift, integer; positive means\n');
    fprintf('                        processed has been moved right with respect to original\n');
    fprintf('               vert_shift is vertical shift, integer; positive means\n');
    fprintf('                        processed has been moved down with respect to original\n');
    fprintf('                        Odd values mean that the processed video has been\n');
    fprintf('                        reframed (i.e., 1st field in time of original, aligns\n');
    fprintf('                        with 2nd field in time of processed -- i.e., +0.5 frame\n');
    fprintf('                        delay)\n');
    fprintf('               delay is time delay in frames, integer; this value adjusts the\n');
    fprintf('                        start frame used by RR_PROCESSED -- it adds "delay" to\n');
    fprintf('                        the 0.8sec starting frame for the processed segment\n');
    fprintf('                        used.\n');
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
    fprintf('  ''results_log'' is the prefix (with path) for text files, where results will\n');
    fprintf('              be written.  If results_log is ''c:\\temp\\525log'', then \n');
    fprintf('              VQM results will be appended to ''c:\\temp\\525log_vqm.txt''\n');
    fprintf('              Errors will append to ''c:\\temp\\525log_error.txt''\n');
    fprintf('              Calibration values append to''c:\\temp\\525log_calibration.txt''\n');
    fprintf('              (Note: Cb and Cr gain & offset are calculated but not removed.)\n');
    fprintf('              Calibration limit exceeded warnings append to\n');
    fprintf('              ''c:\\temp\\525log_limitwarnings.txt''\n');
    fprintf('              Model Parameters append to ''c:\\temp\\525log_parameters.txt''\n');
    fprintf('\n');
    fprintf('  Compressed reduced reference calibration and model features will be read\n');
    fprintf('  from files named after the original video sequence.  The calibration\n');
    fprintf('  features must have "_calibration.mat" appended, in the directory that\n');
    fprintf('  contains that original video sequence.  The model features must have\n');
    fprintf('  "_features.dat" appended, in the directory that contains the original\n');
    fprintf('  video sequence.\n');
    fprintf('EXAMPLE CALL:\n');
    fprintf('  rr_processed ''list_525.txt'' ''525'' ''lowbw'' ''log525''\n');
    fprintf('  rr_processed ''list_625.txt'' ''625'' ''fastlowbw'' ''log625''\n');
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
    fprintf('  The processed segment used is the segment that best aligns with the original,\n');
    fprintf('  within the constraints of the RRNR-TV test plan (or optionally input).\n');
    fprintf('  This software assumes valid video for the following region:\n');
    fprintf('     525-line/NTSC: top=21, left=31, bottom=466, right=690\n');
    fprintf('     625-line/PAL: top=21, left=31, bottom=556, right=690\n');
    fprintf('NOTES:\n');
    fprintf('  The models and calibration used by RR_ORIGINAL and RR_PROCESSED are unchanged\n');
    fprintf('  from those previously released to the public (i.e., CVQM, BVQM).\n');
    fprintf('  Source code and binary executables are available for download\n');
    fprintf('  at www.its.bldrdoc.gov  These algorithms can be freely used for commercial\n');
    fprintf('  and non-commercial applications.\n');
    
    return;
end

% strip off the extra single quotes '' for Windows compile, comment these
% eval lines out for Linux compile
file_list = eval(file_list); 
video_standard = eval(video_standard);
model = eval(model);
results_log = eval(results_log);

calib_log = [results_log '_calibration.txt'];
limit_warning_log = [results_log '_limitwarning.txt'];
error_log = [results_log '_error.txt'];
parameter_log = [results_log '_parameters.txt'];
results_log = [results_log '_vqm.txt'];

fid=fopen(file_list,'r');
fseek(fid,-1,'eof');
tmp = fscanf(fid,'%c');
if tmp ~= sprintf('\n'),
    fclose(fid);
    fid=fopen(file_list,'a');
    fprintf(fid,'\n');
    fclose(fid);
    fid=fopen(file_list,'r');
else
    frewind(fid);
end
warning off;
control_list = textscan(fid, '%s %s %f %f %d %d %d');
warning on;
fclose(fid);

if strcmp(video_standard, '525'),
    video_standard = 'interlace_lower_field_first';
    rows = 486;
    cols = 720;
    fps = 30;
elseif strcmp(video_standard, '625'),
    video_standard = 'interlace_upper_field_first';
    rows = 576;
    cols = 720;
    fps = 25;
else
    error('rr_original called with an invalid string for ''video_standard''.  Aborting');
end

original_list = control_list{:,1};
processed_list = control_list{:,2};
[row,col] = size(control_list);
calib_list = zeros(row, 5);
calib_list(:,:) = nan;
for j=1:5,
    curr_col = control_list{:,j+2};
    for i=1:length(curr_col),
        calib_list(i,j) = curr_col(i);
    end
end
% set delay to 0 if not specified.
for i=1:row,
    if isnan(calib_list(i,5)),
        calib_list(i,5) = 0;
    end
end

for loop=1:length(original_list),
    original_file = original_list{loop};
    processed_file = processed_list{loop};
    calib_values = calib_list(loop,:);

    data_file = sprintf('%s_calibration.mat', original_file);
    model_file = sprintf('%s_features.dat', original_file);
    model_mat_file = sprintf('%s_features.mat', original_file);

    % write name of these files to the error log.
    fid = fopen(error_log,'a');
    fprintf(fid, 'RR_PROCESSED:  %s %s\r\n', original_file, processed_file);
    error_ftell = ftell(fid);
    pause(0.2);
    fclose(fid);

    if strcmp(model,'lowbw'),
        % fine
    elseif strcmp(model,'fastlowbw'),
        % fine!
    elseif strcmp(model,'general'),
        % fine!
    elseif strcmp(model,'developers'),
        % fine!
    else
        cvqm_error(error_log, 1, 'Model request string not recognized.');
        return;
    end


    try
% fprintf('COMMENT IN TRY!\n\n');

        if ~exist(processed_file,'file'),
            cvqm_error(error_log, 3, sprintf('Fatal Error. Processed video file "%s" does not exist.  Clip skipped.\r\n', processed_file));
            continue;
        end
        if ~exist(data_file,'file'),
            cvqm_error(error_log, 3, ...
                sprintf('Fatal Error. Original calibration file "%s" does not exist.  Clip skipped.\r\n', ...
                data_file));
            continue;
        end
        if strcmp(model,'lowbw') || strcmp(model,'fastlowbw'),
            if ~exist(model_file,'file'),
                cvqm_error(error_log, 3, ...
                    sprintf('Fatal Error. Original feature file "%s" does not exist.  Clip skipped.\r\n', ...
                    model_file));
                continue;
            end        
        end
        if strcmp(model,'general') || strcmp(model,'developers'),
            if ~exist(model_mat_file,'file'),
                cvqm_error(error_log, 3, ...
                    sprintf('Fatal Error. Original feature file "%s" does not exist.  Clip skipped.\r\n', ...
                    model_mat_file));
                continue;
            end        
        end

        % load data from original file.
        eval(sprintf('load %s', data_file));
        if ~is_quantize,
            fid = fopen(limit_warning_log,'a');
            fprintf(fid, '\r\nRR_PROCESSED: %s %s\r\n', original_file, processed_file);
            fprintf(fid, '  Warning: uncompressed original calibration features used\r\n');
            pause(0.2);
            fclose(fid);

            fid = fopen(error_log,'a');
            fprintf(fid, '\r\nRR_PROCESSED: %s %s\r\n', original_file, processed_file);
            fprintf(fid, '  Warning: uncompressed original calibration features used\r\n');
            pause(0.2);
            fclose(fid);
        end
        
        if strcmp(model,'lowbw'),
            if is_model ~= 'l',
                cvqm_error(error_log, 3, 'rr_original ran a different model.  Aborting.');
                return;
            end
        elseif strcmp(model,'fastlowbw'),
            if is_model ~= 'f',
                cvqm_error(error_log, 3, 'rr_original ran a different model.  Aborting.');
                return;
            end
        elseif strcmp(model,'general'),
            if is_model ~= 'g',
                cvqm_error(error_log, 3, 'rr_original ran a different model.  Aborting.');
                return;
            end
        elseif strcmp(model,'developers'),
            if is_model ~= 'd',
                cvqm_error(error_log, 3, 'rr_original ran a different model.  Aborting.');
                return;
            end
        end
        if (is_ntsc && rows ~= 486)|| (~is_ntsc && rows ~= 576),
            cvqm_error(error_log, 3, 'rr_original ran a different video standard.  Aborting.');
            return;
        end

        % initialize file read

        dll_video('initialize', 2, processed_file, 'uyvy', video_standard, rows, cols, fps);
        code8s = dll_video('delay_8s', 2, 0);
        dll_calib_video('initialize', 2);
        
        if code8s == 1,
            % print warning.
            cvqm_error(error_log, 2, sprintf('Processed video file "%s" is longer than 8sec. First 8sec used, only.\r\n', processed_file));
        elseif code8s == 2,
            % abort!  cannot use this clip.
            cvqm_error(error_log, 3, sprintf('Fatal Error. Processed video file "%s" is shorter than 8sec. Clip skipped.\r\n', processed_file));
            continue;
        end

        if ~isnan(calib_values(1)),

            % input above estimates 
            y_gain = calib_values(1);
            y_offset = calib_values(2);
            shift.horizontal = calib_values(3);
            shift.vertical = calib_values(4); 
            delay2 = calib_values(5);
            scale.horizontal = 1000;
            scale.vertical = 1000;
            
            dll_calib_video('calibration', shift.horizontal, shift.vertical, cvr, ...
                y_gain, y_offset, scale.horizontal, scale.vertical);

            dll_video('delay_8s', 2, delay2);

%             % write results.
%             [sucess] = cvqm_save_calibration(calibration_file, 'manual', ...
%                 shift, cvr, y_gain, y_offset, scale, delay2, nan, nan, nan, nan); 
%             if sucess == 0,
%                 cvqm_error(error_log, 1, 'Cannot open file to write manual calibration results.');
%                 return;
%             end
            delay1_changed = 0;
            delay2_changed = 0;
            shift_changed = 0;
            y_gain_changed = 0;
            y_offset_changed = 0;
            cb_gain = nan;
            cb_offset = nan; 
            cr_gain = nan; 
            cr_offset = nan;

        else
            % default low bandwidth valid region
            [vr] = dll_default_vr(2);

            % perform low bandwidth temporal registration
            [ti2_proc, ti10_proc, ymean_proc, proc_is_white_clip, proc_is_black_clip] = ...
                dll_lowbw_temporal_features(2, total_sec, vr);

            % note white & black level clipping, if any
            if proc_is_white_clip,
                cvqm_error(error_log, 2, 'White level clipping detected on processed video sequence may cause VQM errors.');
            end
            if proc_is_black_clip,
                cvqm_error(error_log, 2, 'Black level clipping detected on processed video sequence may cause VQM errors.');
            end

            if is_quantize,
                % reconstruct temporal features from quantizer
                [ti2_orig_A, ti10_orig_A, ymean_orig_A] = ...
                  dll_lowbw_temporal_quant(0, ti2_index_A, ti10_index_A, y_index_A);
            end
            
            % complete temporal registration.
            uncert = 0.2; % seconds
            [delay, sucess, is_still] = ...
                dll_lowbw_temporal (2, ti2_orig_A, ti2_proc, ti10_orig_A, ti10_proc, ymean_orig_A, ...
                ymean_proc, uncert, 'field');
            if sucess == 0 && is_still,
                cvqm_error(error_log, 2, 'Still or nearly still sequence; initial temporal registration cannot be computed.');
            elseif sucess == 0,
                cvqm_error(error_log, 2, 'Initial temporal registration algorithm failed.');
            end
            
            % move uncertainty within +/- 3F
            delay1_changed = 0;
            delay1_was = sprintf('Delay before clipping was %d Frames',  delay);
            while delay < -2,
                delay = delay + 1;
                delay1_changed = 1;
            end
            while delay > 2,
                delay = delay - 1;
                delay1_changed = 1;
            end

            % is reframing needed?
            if mod(delay,1) == 0.5,
                dll_calib_video('set_reframe',1);
                delay = delay - 0.5;
            end

            % set delay on processed file.
            dll_video('delay_8s', 2, delay);

            % calculate shift, & scaling 
            if is_quantize,
                [orig_horiz_profile, orig_vert_profile] = ...
                    dll_lowbw_calib_quant(0, index_horiz_profile, index_vert_profile);
            end

            no_scaling = 1;
            [shift, scale,status] = ...
                dll_lowbw_calib_processed(2, seed_state, total_sec, orig_pixels, ...
                orig_horiz_profile, orig_vert_profile, no_scaling);
            if status.scale && ~no_scaling,
                cvqm_error(error_log, 2, 'Actual spatial scaling may be beyond search limits.');
            end
            if status.shift,
                cvqm_error(error_log, 2, 'Actual spatial shift may be beyond search limits.');
            end
            if (strcmp(video_standard,'interlace_lower_field_first') || ...
                    strcmp(video_standard,'interlace_upper_field_first')) && mod(shift.vertical,2),
                cvqm_error(error_log, 2, 'Processed video will be reframed. Actual temporal delay is 0.5 frames greater than reported value.');
            end
            if abs(shift.horizontal) > 8 || abs(shift.vertical) > 5,
                cvqm_error(error_log, 2, 'Extreme spatial shift detected.');
            end
            if scale.vertical ~= 1000 || scale.horizontal ~= 1000,
                cvqm_error(error_log, 2, 'Video scaling detected; please examine other scenes.');
            end
            
            % check limits imposed by RRNR-TV test plan
            shift_changed = 0;
            shift_was = sprintf('Shift before clipping was Horizontal = %d, Vertical = %d\n', ...
                shift.horizontal, shift.vertical);
            if shift.horizontal < -1,
                shift.horizontal = -1;
                shift_changed = 1;
            elseif shift.horizontal > 1,
                shift.horizontal = 1;
                shift_changed = 1;
            end
            if shift.vertical < -1,
                shift.vertical = -1;
                shift_changed = 1;
            elseif shift.vertical > 1,
                shift.vertical = 1;
                shift_changed = 1;
            end

            % input above estimates 
            y_gain = 1.0;
            y_offset = 0.0;
            dll_calib_video('calibration', shift.horizontal, shift.vertical, vr, ...
                y_gain, y_offset, scale.horizontal, scale.vertical);

            % perform low bandwidth valid region
            [pvr] = dll_proc_valid_region(ovr);

            if (pvr.bottom - pvr.top + 1) / rows < 0.55 || ...
                    (pvr.right - pvr.left + 1) / cols < 0.80,
                cvqm_error(error_log, 2, 'Greatly reduced valid region detected.');
            end
            % look for anywhere that PVR is outside of region decided to
            % use (constant)
            if pvr.bottom < cvr.bottom || pvr.right < cvr.right || ...
                    pvr.top > cvr.top || pvr.left > cvr.left,
                cvqm_error(error_log, 2, 'Automatically computed valid region is smaller than region used by model.');
                cvqm_error(error_log, 2, sprintf('Computed valid region: (%d,%d) (%d,%d)', ...
                    pvr.top, pvr.left, pvr.bottom, pvr.right));
                cvqm_error(error_log, 2, sprintf('Constant valid region: (%d,%d) (%d,%d)', ...
                    cvr.top, cvr.left, cvr.bottom, cvr.right));
            end

            %%%%%%!!!!!!!%%%%%%%%%%
            % discard PVR and use CVR instead
            %%%%!!!!!!!%%%%%%%%%%%%%%
            pvr = cvr;

            % input improved PVR estimate 
            dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
                y_gain, y_offset, scale.horizontal, scale.vertical);

            % estimate luminance gain & offset
            if is_quantize,
                [orig_y, orig_cb, orig_cr] = dll_lowbw_gain_v2_quant(0, index_y, index_cb, index_cr); % reconstruct
            end
            
            [y_gain, y_offset, cb_gain, cb_offset, cr_gain, cr_offset, sucess] = ...
                dll_lowbw_gain_v2_processed(2, total_sec, orig_y, orig_cb, orig_cr, yesno);
%             cvqm_error(error_log, 2, 'Color Gain & Offset estimated but not removed.');

            if y_gain < 0.9 || y_gain > 1.1,
                cvqm_error(error_log, 2, 'Extreme Luminance Gain detected.');
            end
            if y_offset < -20 || y_offset > 20,
                cvqm_error(error_log, 2, 'Extreme Luminance Offset detected.');
            end
            if sucess == 0,
                cvqm_error(error_log, 2, 'Warning: algorithm used to estimate Luminance Gain & Offset may have failed.');
            end
            if sucess == -1,
                cvqm_error(error_log, 2, 'Luminance gain & offset algorithm detected extreme values or failed; results discarded.');
            end
            if isnan(cb_gain),
                cvqm_error(error_log, 2, 'Warning: algorithm used to estimate Cb Gain & Offset failed.');
            end
            if isnan(cr_gain),
                cvqm_error(error_log, 2, 'Warning: algorithm used to estimate Cr Gain & Offset failed.');
            end
            
            % check RRNR-TV test plan limits.
            y_gain_was = sprintf('Y Gain before clipping was %f\n', y_gain);
            y_gain_changed = 0;
            if y_gain < 0.97,
                y_gain_changed = 1;
                y_gain = 0.97;
            elseif y_gain > 1.03,
                y_gain_changed = 1;
                y_gain = 1.03;
            end
            
            y_offset_changed = 0;
            y_offset_was = sprintf('Y Offset before clipping was %f\n', y_offset);
            if y_offset < -10,
                y_offset_changed = 1;
                y_offset = -10;
            elseif y_offset > 10,
                y_offset_changed = 1;
                y_offset = 10;
            end

            % input luminance gain & offset estimates
            dll_calib_video('calibration', shift.horizontal, shift.vertical, pvr, ...
                y_gain, y_offset, scale.horizontal, scale.vertical);

            if is_quantize,
                % reconstruct temporal features from quantizer
                [ti2_orig_B, ti10_orig_B, ymean_orig_B] = ...
                  dll_lowbw_temporal_quant(0, ti2_index_B, ti10_index_B, y_index_B);
            end

            % re-set delay on processed file to zero
            dll_video('delay_8s', 2, 0);

            % perform low bandwidth temporal registration a second time!
            [ti2_proc, ti10_proc, ymean_proc, proc_is_white_clip, orig_is_black_clip] = ...
                dll_lowbw_temporal_features(2, total_sec, pvr);

            uncert = 1;
            [delay2, sucess, is_still] = ...
                dll_lowbw_temporal (2, ti2_orig_B, ti2_proc, ti10_orig_B, ti10_proc, ymean_orig_B, ymean_proc, uncert);

            if sucess == 0 && is_still,
                cvqm_error(error_log, 2, 'Still or nearly still sequence; temporal registration cannot be computed.');
            elseif sucess == 0,
                cvqm_error(error_log, 2, 'Final temporal registration algorithm failed.');
            end

            % move uncertainty within +/- 3F
            delay2_changed = 0;
            delay2_was = sprintf('Delay before clipping was %d Frames',  delay2);
            while delay2 < -2,
                delay2 = delay2 + 1;
                delay2_changed = 1;
            end
            while delay2 > 2,
                delay2 = delay2 - 1;
                delay2_changed = 1;
            end
            
            % set delay on processed file.
            dll_video('delay_8s', 2, delay2);

%             % write results.
%             [sucess] = cvqm_save_calibration(calibration_file, 'rrcal2', shift, actual_pvr, y_gain, y_offset, scale, delay2, cb_gain, cb_offset, cr_gain, cr_offset); 
%             if sucess == 0,
%                 cvqm_error(error_log, 1, 'Cannot open file to write lowbw calibration results.');
%                 return;
%             end
        end


        if strcmp(model,'lowbw'),

            proc_features = dll_features('Low', 2, total_sec);
            [vqm, pars, par_names] = dll_model('vqm',model_file, proc_features);
%             [sucess] = cvqm_save_model(sprintf('%s_model.txt', processed_file), model, vqm, pars, par_names);

        elseif strcmp(model,'fastlowbw'),
            proc_features = dll_features('Fast', 2, total_sec);
            [vqm, pars, par_names] = dll_model('vqm',model_file, proc_features);
%             [sucess] = cvqm_save_model(sprintf('%s_model.txt', processed_file), model, vqm, pars, par_names);

        elseif strcmp(model,'general'),
            eval( sprintf('load %s', model_mat_file));
            proc_features = dll_features('General', 2, total_sec);
            [vqm, pars, par_names] = dll_model('vqm',orig_features, proc_features);
%             [sucess] = cvqm_save_model(sprintf('%s_model.txt', processed_file), model, vqm, pars, par_names);

        elseif strcmp(model,'developers'),
            eval( sprintf('load %s', model_mat_file));
            proc_features = dll_features('Developers', 2, total_sec);
            [vqm, pars, par_names] = dll_model('vqm',orig_features, proc_features);
%             [sucess] = cvqm_save_model(sprintf('%s_model.txt', processed_file), model, vqm, pars, par_names);

        end

        % write features to feature log
        fid = fopen(parameter_log,'a');
        if loop == 1,
            fprintf(fid, 'Processed_file VQM-%s ', model);
            for i=1:length(par_names),
                fprintf(fid, '%s ', par_names{i});
            end
            fprintf(fid, '\r\n');
        end
        fprintf(fid, '%s %f ', processed_file, vqm);
        for i=1:length(pars),
            fprintf(fid, '%f ', pars(i));
        end
        fprintf(fid, '\r\n');
        pause(0.2);
        fclose(fid);
            
        
% fprintf('remove next lines after debugged\r\n');
% [orig_features.si_orig, orig_features.part_si_min, orig_features.part_si_max, ...
%      orig_features.hv_feat_orig, orig_features.part_hv_min, orig_features.part_hv_max, orig_features.y_orig, ...
%      orig_features.cb_orig, orig_features.cr_orig, orig_features.part_c_min, orig_features.part_c_max, orig_features.part_c, ...
%      orig_features.ati_orig, orig_features.part_ati_min, orig_features.part_ati_max, orig_features.part_ati, orig_features.code_ati ] ...
%      = model_lowbw_compression('uncompress', model_file);
% save 'features.mat' proc_features orig_features;

        fid = fopen(results_log,'a');
        if loop == 1,
            fprintf(fid, 'Processed_file VQM-%s ', model);
            fprintf(fid, '\r\n');
        end

        fprintf(fid, '%s %f\r\n', processed_file, vqm);
        pause(0.2);
        fclose(fid);

        fid = fopen(calib_log,'a');
        if loop == 1,
            fprintf(fid, '\r\nProcessed_file VQM-%s Y_gain Y_offset Shift-Horizontal Shift-Vertical Delay-frames  Cb-Gain Cb-Offset Cr-Gain Cr-Offset\r\n', model);
        end
        fprintf(fid, '%s %f  %f %f %d %d %d   %f %f %f %f\r\n', processed_file, vqm, ...
            y_gain, y_offset, shift.horizontal, shift.vertical, delay2, ...
            cb_gain, cb_offset, cr_gain, cr_offset);
        pause(0.2);
        fclose(fid);

        fid = fopen(error_log,'a');
        if error_ftell == ftell(fid),
            fprintf(fid, '  No errors encountered.\r\n');
        end
        pause(0.2);
        fclose(fid);

        
        if delay1_changed || delay2_changed || shift_changed || y_gain_changed || y_offset_changed,
            fid = fopen(limit_warning_log,'a');
            fprintf(fid, '\r\nRR_PROCESSED: %s\r\n', processed_file);
            if delay1_changed,
                fprintf(fid, '  1st pass calculation of delay was beyond 2F; delay clipped\r\n');
                fprintf(fid, '  %s\r\n', delay1_was);
            end
            if delay2_changed,
                fprintf(fid, '  1st pass calculation of delay was beyond 2F; delay clipped\r\n');
                fprintf(fid, '  %s\r\n', delay2_was);
            end
            if shift_changed,
                fprintf(fid, '  Spatial shift was beyond +/- 1; shift clipped\r\n');
                fprintf(fid, '  %s\r\n', shift_was);
            end
            if y_gain_changed,
                fprintf(fid, '  Y gain beyond [0.97,1.03]; gain clipped.\r\n');
                fprintf(fid, '  %s\r\n', y_gain_was);
            end
            if y_offset_changed,
                fprintf(fid, '  Y offset beyond [-10,+10]; offset clipped.\r\n');
                fprintf(fid, '  %s\r\n', y_offset_was);
            end
            pause(0.2);
            fclose(fid);
        end

    catch
        cvqm_error(error_log, 1, lasterr);
    end
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





