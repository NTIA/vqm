function [delay, sucess, is_still] = dll_lowbw_temporal (fn, ti2_orig, ti2_proc, ...
    ti10_orig, ti10_proc, ymean_orig, ymean_proc, uncert, varargin);
% DLL_LOWBW_TEMPORAL
%   Step 2: Compute delay
% SYNTAX
%   [delay, sucess, is_still] = dll_lowbw_temporal (ti2_orig, ti2_proc,
%       ti10_orig, ti10_proc, ymean_orig, ymean_proc, progressive, uncert);
%   [...] = dll_lowbw_temporal(...,'Flag',...);
% DESCRIPTION
%   Input arguments are 'fn' from dll_video, file ID for either original or
%   processed -- it doesn't matter which.  The next six input arguments are
%   the three features computed on the original video by function
%   dll_lowbw_temporal_features (named ti2_orig, ti10_orig & ymean_orig); 
%   and the features computed on the processed video by function
%   dll_lowbw_temporal_features (named ti2_proc, ti10_proc & ymean_proc).  
%   And finally, the temporal registration uncertainty, 'uncert', in seconds.
%
%   Optional Flags:
%       'field'     For interlaced systems, align to field accuracy.
%       'frame'     For interlaced systems, align to frame accuracy.  Default.
%
%   Return 'delay', the temporal registration delay in frames; 'sucess',
%   which is 1 if the algorithm succeeded & 0 if the algorithm failed; and
%   'is_still' which contains 1 if the video sequence appears to be still
%   or nearly still (thus temporal registration will always fail), and 0
%   otherwise. 

frame_select = 1;

cnt = 1;
while cnt <= length(varargin),
    if strcmpi(varargin{cnt},'field'),
        frame_select = 0;
        cnt = cnt + 1;
    elseif strcmpi(varargin{cnt},'frame'),
        frame_select = 1;
        cnt = cnt + 1;
    else
        error('optional flag not recognized');
    end
end

[video_standard] = dll_video('get_video_standard', fn);
if strcmp(video_standard,'interlace_upper_field_first'), % e.g., if fps == 25,
    fld_num(1) = 1;
    fld_num(2) = 2;
    progressive = 0;
elseif strcmp(video_standard,'interlace_lower_field_first'),
    fld_num(1) = 2;
    fld_num(2) = 1;
    progressive = 0;
elseif strcmp(video_standard,'progressive'),
    fld_num(1) = 1;
    progressive = 1;
else
    warning('video standard not recognized');
    fld_num(1) = 1;
    progressive = 1;
end

uncert = ceil(uncert * dll_video('fps', fn));
if ~progressive,
    uncert = uncert * 2;
end

[sucess, is_still, delay] = trc_align_with_three_features(ti2_orig, ti2_proc, ...
    ymean_orig, ymean_proc, ti10_orig, ti10_proc, uncert, progressive, frame_select);


%********************************************************************************
%********************************************************************************
%********************************************************************************
function [diff, is_valid] = trc_correlate_one_feature(feature_orig, feature_proc, ...
    threshold, uncert, is_progressive, frame_select, hold_name, file_name);

    diff = NaN;

    is_valid = 1;
    
    hold_length = min( length(feature_proc), length(feature_proc)) -2*uncert;
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
function [is_valid, is_still, delay] = trc_align_with_three_features(ti2_orig, ti2_proc, ...
    y_orig, y_proc, ti10_orig, ti10_proc, uncert, is_progressive, frame_select);

 	STILL_TI = 0.15;
 	STILL_Y = 0.25;

    [diff_ti2,valid_ti2] = trc_correlate_one_feature(ti2_orig, ti2_proc, STILL_TI, ...
        uncert, is_progressive, frame_select);
    [diff_y,valid_y] = trc_correlate_one_feature(y_orig, y_proc, STILL_Y, uncert, ...
        is_progressive, frame_select);
    [diff_ti10,valid_ti10] = trc_correlate_one_feature(ti10_orig, ti10_proc, STILL_TI, ...
        uncert, is_progressive, frame_select);

    [is_valid, is_still, delay] = trc_align_combine(diff_ti2, diff_y, diff_ti10, uncert, ...
        valid_ti2, valid_y, valid_ti10, ...
        is_progressive, frame_select);

%********************************************************************************
function [is_valid, is_still, delay] = trc_align_combine(diff_ti2, diff_y, diff_ti10, uncert, ...
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
        delay = 0;
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
    
    if ~passvalid_ti2 & ~passvalid_y & ~passvalid_ti10,
        is_still = 1;
    else
        is_still = 0;
    end


