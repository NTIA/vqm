function [new_clip_structs, status] = temporal_registration_lowbw(test_structs, clip_structs,feature_base_dir, varargin)
% TEMPORAL_REGISTRATION
%   Automatic computation of temporal registration, 2005 low-bandwidth
%   algorithm.
% SYNTAX
%  [new_clip_structs, status] = temporal_registration_lowbw(test_structs, clip_structs, feature_base_dir);
%  [...] = temporal_registration(...,'PropertyName',...)
% DESCRIPTION
%  Perform temporal registration using a low bandwidth.
%
%  'feature_base_dir' is the directory containing relevant field-based features.
%  (Warning:  if features pre-computed, then temporal registration & other
%  calibration values in clip_structs will be ignored!)
%
%  Temporal registration (etc.) within clip_structs will be used as the initial
%  temporal registration for the algorithm.
%
%  Optional properties:
%   'frame_select'           Default.  Align to frame/reframe indicated by
%                            spatial registration. 
%   'field_select'           Align to field accuracy, ignoring the 
%                            frame/reframe spatial registration.  No impact
%                            on progressive sequences.
%   'Uncertainty',value,     Assume temporal registration uncertainty of
%                            'value' seconds.  By default, 1/2 second.
%   'verbose',               Print status information to screen
%   'quiet',                 Print no information to screen.
%   'log',file1,             Create log files.  Echo to 'file1' the
%                            default screen summary.  Append '.mat' and
%                            write out matlab variable holding alignment results.
%   'aligned'                Use current temporal registration as the
%                            starting point.
%   'unaligned'              Default.  No temporal registration available.
%                            Use 'first frames align' as an initial guess. 
%   'Quantize'               Quantize original video features, as if
%                            operating in-service.
%
%  Return value is new_clip_struct, reflects clip_struct but with modified
%  temporal registration.  Also return status,
%  which has the following fields defined:
%     'status.error'              0 if okay; 1 if a fatal error was encountered.
%                                 2 if original is missing for a scene.  
%                                 3 if multiple original sequences are
%                                 defined for one scene.  If status.error is
%                                 not 0, then the other structure fields will 
%                                 be unavailable.  See also function "lasterr"
%     'status.uncertainty'        Number of video clips that should be re-run
%                                 with a larger uncertainty.
%     'status.still'              Number of video clips that came from still scenes

try

    % handle optional properties.
    uncert_sec = 0.5;
    is_verbose = 0;
    aligned_flag = 0;
    log = 0;
    frame_select = 1;
    is_quantize = 0;

    cnt = 1;
    while cnt <= nargin-3,
        if strcmpi(varargin{cnt},'uncertainty'),
            uncert_sec = varargin{cnt+1};
            cnt = cnt + 2;
        elseif strcmpi(varargin{cnt},'quantize'),
            is_quantize = 1;
            cnt = cnt + 1;
        elseif strcmpi(varargin{cnt},'verbose'),
            is_verbose = 1;
            cnt = cnt + 1;
        elseif strcmpi(varargin{cnt},'quiet'),
            is_verbose = 0;
            cnt = cnt + 1;
        elseif strcmpi(varargin{cnt},'log'),
            if nargin-1 < cnt+1,
                error('Specify log file name');
            end
            log = 1;
            logfile1 = fopen(varargin{cnt+1},'w');
            logfile2_name = [varargin{cnt+1} '.mat'];
            if logfile1 == -1,
                error('Could not open log file "%s"', varargin{cnt+1});
            end
            cnt = cnt + 2;
        elseif strcmpi(varargin{cnt},'frame_select'),
            frame_select = 1;
            cnt = cnt + 1;
        elseif strcmpi(varargin{cnt},'field_select'),
            frame_select = 0;
            cnt = cnt + 1;
        elseif strcmpi(varargin{cnt},'aligned'),
            aligned_flag = 1;
            cnt = cnt + 1;
        elseif strcmpi(varargin{cnt},'unaligned'),
            aligned_flag = 0;
            cnt = cnt + 1;
        else
            error('Property Name not recognized.  Aborting.');
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Define the quantizers for the three temporal registration features.
    %  These designs are for 10 bit linear quantizers.
    
    %  cont feature
    start = 0.0;  % first code
    last = 255.0;  % 252 is the maximum observed in the training data
    high_codes = 4096;  % number of codes for 12-bit quantizer
    code_cont = start:(last-start)/(high_codes-1):last;
    %  Generate the partitions, halfway between codes
    code_cont_size = size(code_cont,2);
    part_cont = (code_cont(2:code_cont_size)+code_cont(1:code_cont_size-1))/2;
    
    %  ti2 feature
    start = 0.0;  % first code
    last = 210.0;  % 209 is the maximum observed in the training data
    high_codes = 4096;  % number of codes for 12-bit quantizer
    code_ti2 = start:(last-start)/(high_codes-1):last;
    %  Generate the partitions, halfway between codes
    code_ti2_size = size(code_ti2,2);
    part_ti2 = (code_ti2(2:code_ti2_size)+code_ti2(1:code_ti2_size-1))/2;
    
    %  ti10 feature
    start = 0.0;  % first code
    last = 210.0;  % 209 is the maximum observed in the training data
    high_codes = 4096;  % number of codes for 12-bit quantizer
    code_ti10 = start:(last-start)/(high_codes-1):last;
    %  Generate the partitions, halfway between codes
    code_ti10_size = size(code_ti10,2);
    part_ti10 = (code_ti10(2:code_ti10_size)+code_ti10(1:code_ti10_size-1))/2;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute features.  This may be fast, if the features are already there.
    if is_verbose,
        verbose_string = 'verbose';
    else
        verbose_string = 'quiet';
    end

    % check incoming structure.
    [bool] = check_clips(clip_structs, test_structs, verbose_string);
    if bool,
        error('clip structure invalid, aborting');
    end
    
    if ~aligned_flag,
        tmp=sort_clips_by('scene',clip_structs, test_structs);
        for cnt=1:length(tmp),
            hold = tmp{cnt};
            min_len = inf;
            for loop = hold,
                min_len = min(min_len, clip_structs(loop).loc_stop-clip_structs(loop).loc_start);
            end
            for loop = hold,
                clip_structs(loop).align_start = clip_structs(loop).loc_start;
                clip_structs(loop).align_stop = clip_structs(loop).loc_start + min_len;
            end
        end
    end

    feature_temporal_registration_lowbw(test_structs, clip_structs, feature_base_dir, verbose_string);

    % sort clips by scene.
    offsets = sort_clips_by('scene', clip_structs, test_structs);

%     % clear any luminance gain/offset settings.
%     for i=1:length(clip_structs),
%         clip_structs(i).luminance_gain = 1;
%         clip_structs(i).luminance_offset = 0;
%     end

    % replicate clip structure.
    new_clip_structs = clip_structs;

    % keep track of statistics
    uncert_too_small = 0;
    align_ambiguous = 0;
    still_scene = 0;
    invalid_scene = 0;
    reframe_issue = 0;

    % Loop through each scene.
    clip_number = 1;
    for cnt = 1:length(offsets),
        curr_offsets = offsets{cnt};
        clip = curr_offsets(1);

        % Check if this is an interlace or progressive sequence
        if strcmp(clip_structs(clip).video_standard,'progressive'),
            is_progressive = 1;
        elseif strcmp(clip_structs(clip).video_standard,'interlace_lower_field_first') || strcmp(clip_structs(clip).video_standard,'interlace_upper_field_first'),
            is_progressive = 0;
        else
            error('video standard not recognized');
        end

        % read in original features.
        [orig_ti2] = read_clip_feature(clip_structs(curr_offsets(1)), feature_base_dir, 'feature_Y_ti2_field_rms');
        [orig_ti10] = read_clip_feature(clip_structs(curr_offsets(1)), feature_base_dir, 'feature_Y_ti10_field_rms');
        [orig_y] = read_clip_feature(clip_structs(curr_offsets(1)), feature_base_dir, 'feature_Y_cont_field_mean');

        uncert = round(uncert_sec * clip_structs(curr_offsets(1)).fps);
        if ~is_progressive,
            uncert = uncert * 2;
        end

        if is_quantize,
            %  Quantize original features
            [index_ti2] = quantiz_fast(squeeze(orig_ti2)',part_ti2);
            orig_ti2 = code_ti2(1+index_ti2);
            orig_ti2 = reshape(orig_ti2,1,1,length(orig_ti2));
            
            [index_ti10] = quantiz_fast(squeeze(orig_ti10)',part_ti10);
            orig_ti10 = code_ti10(1+index_ti10);
            orig_ti10 = reshape(orig_ti10,1,1,length(orig_ti10));

            [index_cont] = quantiz_fast(squeeze(orig_y)',part_cont);
            orig_y = code_cont(1+index_cont);
            orig_y = reshape(orig_y,1,1,length(orig_y));
        end
        
        % loop through each processed version of this scene.
        for loop = 2:length(curr_offsets),
            clip = curr_offsets(loop);

            % read in features.
            [proc_ti2] = read_clip_feature(clip_structs(clip), feature_base_dir, 'feature_Y_ti2_field_rms');
            [proc_ti10] = read_clip_feature(clip_structs(clip), feature_base_dir, 'feature_Y_ti10_field_rms');
            [proc_y] = read_clip_feature(clip_structs(clip), feature_base_dir, 'feature_Y_cont_field_mean');

            [valid, still, delay] = align_with_three_features(...
                orig_ti2, proc_ti2, orig_y, proc_y, orig_ti10, proc_ti10, ...
                uncert, is_progressive, frame_select, ...
                clip_structs(clip).test{1}, clip_structs(clip).scene{1}, clip_structs(clip).hrc{1}, ...
                clip_number);

            if still,
                delay = 0;
                still_scene = still_scene + 1;
            elseif ~valid,
                delay = 0;
                invalid_scene = invalid_scene + 1;
            end
            if abs(delay) == uncert || abs(delay) == uncert-1,
                uncert_too_small = uncert_too_small + 1;
            end
            if mod(delay,1) > 0,
                reframe_issue = reframe_issue + 1;
            end

            % echo results to log
            if log,         
               fprintf(logfile1,'%s:%s(%s) delay=%d valid=%d still=%d\n', clip_structs(clip).test{1}, ...
                   clip_structs(clip).scene{1}, clip_structs(clip).hrc{1}, delay, valid, still );
               results(clip_number).test = clip_structs(clip).test;
               results(clip_number).scene = clip_structs(clip).scene;
               results(clip_number).hrc = clip_structs(clip).hrc;
               results(clip_number).delay = delay;
               results(clip_number).valid = valid;
               results(clip_number).still = still;
            end
            clip_number = clip_number + 1;

            if is_verbose,
                fprintf('\t%s:%s(%s)  delay %f',  ...
                    clip_structs(clip).test{1}, clip_structs(clip).scene{1}, clip_structs(clip).hrc{1}, delay );
                if ~valid,
                    if still,
                        fprintf('   Still or near still');
                    else
                        fprintf('   Cannot compute delay');
                    end
                end
                fprintf('\n');
            end

            clip_delay(loop) = delay;
        end

        %
        clip_delay(1) = 0;
        clip_delay = -clip_delay;
        orig_extra_start = max(clip_delay(2:length(curr_offsets)));
        new_clip_structs(curr_offsets(1)).align_start = new_clip_structs(curr_offsets(1)).align_start + orig_extra_start;
        clip_delay = orig_extra_start - clip_delay;
        
        if mod(clip_delay(1),1),
            clip_delay = clip_delay + 0.5;
        end

        % record overlapping range for all clips of this scene.
        % Cut off the beginning of the clip, to cause clips to align
        for loop = 1:length(curr_offsets),
            clip = curr_offsets(loop);
            new_clip_structs(clip).align_start = clip_structs(clip).align_start + clip_delay(loop);
            new_clip_structs(clip).align_stop = clip_structs(clip).align_stop + clip_delay(loop);
            is_length(loop) = new_clip_structs(clip).loc_stop - new_clip_structs(clip).align_start;
        end

        % shorten clip as needed by file length
        use_length = floor(min(is_length));
        for loop = 1:length(curr_offsets),
            clip = curr_offsets(loop);
            new_clip_structs(clip).align_stop = new_clip_structs(clip).align_start + use_length;
            if mod(new_clip_structs(clip).align_start,1) ~= 0,
                new_clip_structs(clip).align_stop = new_clip_structs(clip).align_stop + 0.5;
            end
        end
        
        % if original's align_start is after loc_start, move all clips
        % inward.
        clip_start = new_clip_structs(curr_offsets(1)).loc_start - new_clip_structs(curr_offsets(1)).align_start;
        clip_start = max(clip_start, 0);
        for loop = curr_offsets,
            new_clip_structs(loop).align_start = new_clip_structs(loop).align_start + clip_start;
        end


    end
    
    if log,
        save(logfile2_name,'results');
    end

    status.error = 0;
    status.uncertainty = uncert_too_small;
    status.still = still_scene;
    status.reframe = reframe_issue;

    if is_verbose,
        fprintf('\n\n%d clips should be re-run with a larger uncertainty\n', uncert_too_small);
        fprintf('%d clips came from still or nearly still scenes\n', still_scene);
        fprintf('%d clips contained motion but could not be aligned\n', invalid_scene);
        fprintf('Warning: %d clips indicated a reframing contradictory to spatial shift.  Alignment invalid\n', reframe_issue);
    end
    if log,
        fprintf(logfile1, '\n\n%d clips should be re-run with a larger uncertainty\n', uncert_too_small);
        fprintf(logfile1,'%d clips came from still or nearly still scenes\n', still_scene);
        fprintf(logfile1,'%d clips contained motion but could not be aligned\n', invalid_scene);
        fprintf(logfile1, 'Warning: %d clips indicated a reframing contradictory to spatial shift.  Alignment invalid\n', reframe_issue);
    end

    if is_verbose,
        fprintf('\n\nWarning:  segments selected will imperfectly match those originally specified.\n');
        fprintf('A manual check is recommended, if the intent was to select an exact portion\n');
        fprintf('of the available video.\n');
    end
    if log,
        fprintf(logfile1, '\n\nWarning:  segments selected will imperfectly match those originally specified.\n');
        fprintf(logfile1, 'A manual check is recommended, if the intent was to select an exact portion\n');
        fprintf(logfile1, 'of the available video.\n');
    end

    if log,
        fclose(logfile1);
    end

    if is_verbose,
        fprintf('\n\nWARNING:  Run check_clips on the returned clip structure, to ensure validity.\n');
        fprintf('This temoral registration function does not check or validate reframing.\n');
    end
    
    % Make sure that the end point is correct for both frame & reframe
    % alignments.
    if frame_select,
        new_clip_structs = fix_temporal(test_structs, new_clip_structs, 'Spatial');
    end
catch
    status.error = 1;
    if is_verbose,
        fprintf('%s\n', lasterr);
    end
    new_clip_structs = [];
end

%********************************************************************************
function [diff, is_valid] = correlate_one_feature(feature_orig, feature_proc, threshold, uncert, is_progressive, frame_select, hold_name, file_name);

    diff = NaN;

    is_valid = 1;
    
    hold_length = length(feature_proc)-2*uncert;
    proc = squeeze(feature_proc(uncert+1:uncert+hold_length));
    hold_std = std(proc);
    
    if hold_std < threshold,
        is_valid = 0;
        hold_std = 1;
    end
    proc = proc ./ hold_std;
    
    if is_progressive,
        for delay = 1:uncert*2+1,
            src = squeeze(feature_orig(delay:delay+hold_length-1));
            hold_std = std(src);
            if hold_std < threshold,
                is_valid = 0;
                return;
            end
            diff(delay) = std(src ./ hold_std - proc);
        end
    else
        if frame_select,
            for delay = 1:2:uncert*2+1,
                src = squeeze(feature_orig(delay:delay+hold_length-1));
                hold_std = std(src);
                if hold_std < threshold,
                    is_valid = 0;
                    hold_std = 1;;
                end
                diff( (delay+1)/2 ) = std(src ./ hold_std - proc);
                
            end
        else
            % align to nearest field -- may indicate different spatial
            % shift.
            for delay = 1:uncert*2+1,
                src = squeeze(feature_orig(delay:delay+hold_length-1));
                hold_std = std(src);
                if hold_std < threshold,
                    is_valid = 0;
                    hold_std = 1;
                end
                diff( delay ) = std(src ./ hold_std - proc);
            end
        end
    end
    
    

%********************************************************************************
function [is_valid, is_still, delay] = align_with_three_features(ti2_orig, ti2_proc, ...
    y_orig, y_proc, ti10_orig, ti10_proc, uncert,...
    is_progressive, frame_select, test, scene, hrc, clip_number);


 	STILL_TI = 0.15;
 	STILL_Y = 0.25;

    [diff_ti2,valid_ti2] = correlate_one_feature(ti2_orig, ti2_proc, STILL_TI, uncert, is_progressive, frame_select);
    [diff_y,valid_y] = correlate_one_feature(y_orig, y_proc, STILL_Y, uncert, is_progressive, frame_select);
    [diff_ti10,valid_ti10] = correlate_one_feature(ti10_orig, ti10_proc, STILL_TI, uncert, is_progressive, frame_select);

    [is_valid, is_still, delay] = align_combine(diff_ti2, diff_y, diff_ti10, uncert, ...
        valid_ti2, valid_y, valid_ti10, ...
        is_progressive, frame_select, test, scene, hrc, clip_number);

%********************************************************************************
function [is_valid, is_still, delay] = align_combine(diff_ti2, diff_y, diff_ti10, uncert, ...
    passvalid_ti2, passvalid_y, passvalid_ti10,  ...
    is_progressive, frame_select, test, scene, hrc, clip_number);

    
    CONST_VALID = 0.25;
    CONST_INVALID = 1.4;
    CONST_DELTA = 0.04;
    CONST_TI_RANGE = 3;
    CONST_Y_RANGE = 4;
    
    % judge whether each of the three features is valid
    if passvalid_ti2 & min(diff_ti2) < CONST_INVALID,
        range_high = find(diff_ti2 <= min(diff_ti2) + CONST_DELTA);
        range_high = range_high(length(range_high)) - range_high(1) + 1;
        if min(diff_ti2) < CONST_VALID | range_high <= CONST_TI_RANGE,
            valid_ti2 = 1;
        else
            valid_ti2 = 0;
        end
    else
        valid_ti2 = 0;
    end

    if passvalid_y & min(diff_y) < CONST_INVALID,
        range_high = find(diff_y <= min(diff_y) + CONST_DELTA);
        range_high = range_high(length(range_high)) - range_high(1) + 1;
        if min(diff_y) < CONST_VALID | range_high <= CONST_Y_RANGE,
            valid_y = 1;
        else
            valid_y = 0;
        end
    else valid_y = 0;
    end

    if passvalid_ti10 & min(diff_ti10) < CONST_INVALID,
        range_high = find(diff_ti10 <= min(diff_ti10) + CONST_DELTA);
        range_high = range_high(length(range_high)) - range_high(1) + 1;
        if min(diff_ti10) < CONST_VALID | range_high <= CONST_TI_RANGE,
            valid_ti10 = 1;
        else
            valid_ti10 = 0;
        end
    else
        valid_ti10 = 0;
    end


    % make composite plot of valid differences/correlations
    is_valid = 1;
    if valid_ti2 & valid_y & valid_ti10,
        diff = (diff_ti2+diff_y+diff_ti10)/3;
    elseif valid_ti2 & valid_y,
        diff = (diff_ti2+diff_y)/2;
    elseif valid_y & valid_ti10,
        diff = (diff_y+diff_ti10)/2;
    elseif valid_ti2 & valid_ti10,
        diff = (diff_ti2+diff_ti10)/2;
    elseif valid_ti2,
        diff = diff_ti2;
    elseif valid_y,
        diff = diff_y;
    elseif valid_ti10,
        diff = diff_ti10;
    else
        is_valid = 0;
        delay = NaN;
    end


    if is_valid,
        % find delay that minimizes correlation
        [mindiff,mindelay] = min(diff);
        if is_progressive,
            delay = - (mindelay - 1 - uncert);
        else
            if frame_select,
                delay = - (mindelay - 1 - uncert/2);
            else
                mindelay = (mindelay+1)/2;
                delay = - (mindelay - 1 - uncert/2);
            end
        end
    end    
    %
    

%     % plot correlations
%     plot(diff_ti2,'r'); 
%     hold on; 
%     plot(diff_y,'g'); 
%     plot(diff_ti10,'m');  
%     if is_valid, plot(diff,'--k'); end
%     hold off;
%     if is_valid,
%         title(sprintf('Delay = %d',delay)), 
%     else
%         title('Delay Ambiguous -- Invalid');
%     end
%     xlabel(sprintf('red=ti2(%d), green=ti10(%d), magenta=Y(%d), black=valid\n', valid_ti2,valid_y,valid_ti10));
%     0; 

    if ~passvalid_ti2 & ~passvalid_y & ~passvalid_ti10,
        is_still = 1;
    else
        is_still = 0;
    end


