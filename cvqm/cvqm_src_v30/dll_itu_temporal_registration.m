function [delay, sucess, is_still, is_ambiguous, is_uncert] = dll_itu_temporal_registration(uncert_sec);
% DLL_ITU_TEMPORAL_REGISTRATION 
%   Compute frame-based temporal registration.
% SYNTAX
%  [delay, sucess, is_still, is_ambiguous, is_uncert] = dll_itu_temporal_registration;
%  [...] = dll_itu_temporal_registration(uncert);
% DESCRIPTION
%  Perform temporal registration.
%  Use standard temporal registration algorithm from 2003 NTIA report 
%  & ANSI Rec. 801.03, frame based algorithm.
%
%  Precondition:  dll_video must have been initialized with fn=1 for the
%  original video file, and fn=2 for the processed video file.
%  Additionally, dll_calib_video must have been initialized and set to
%  calibration values for fn=2, particularly processed valid region (PVR),
%  luminance gain & offset, and spatial registration.
%
%  Optional input argument 'uncert' contains the alignment uncertanty, in
%  seconds.  By default 'uncert' = 1.0 (one second).
%
%  Return values are defined as follows:
%   'delay'  Delay between original and processed video, in frames.  This
%   return value is for informational purposes only; dll_video has already
%   been adjusted appropriately.
%   'sucess' contains 1 if the algorithm succeeded & 0 if the algorithm failed.
%   'is_still' which contains 1 if the video sequence appears to be still
%   or nearly still (thus temporal registration will always fail), and 0
%   otherwise.
%   'is_ambiguous' contains 1 if temporal registration is ambiguous, 0
%   otherwise.
%   'is_uncert' contains 1 if temporal registration should be re-run with a
%   larger temporal uncertainty.
%
%  Temporal delay calculated but not removed from sequences.

warning('off','MATLAB:max:mixedSingleDoubleInputs');

sucess = 0;
% try
    % handle optional properties.
    if ~exist('uncert_sec','var'),
        uncert_sec = 1.0;
    end

    % set subsampling size & other constants.
    hsize = 16;
    vsize = 16;
    STILL_THRESHOLD = 0.002;
    HFW = 3;
    BELOW_WARN = 0.9;
    DELTA = 4;

% 'unnecesary, but a difference'
% %     % clear any luminance gain/offset settings.
% %     for i=1:length(clip_structs),
% %         clip_structs(i).luminance_gain = 1;
% %         clip_structs(i).luminance_offset = 0;
% %     end

    fn_orig = 1;
    fn_proc = 2;

    % keep track of statistics
    is_uncert = 0;
    is_ambiguous = 0;
    is_still = 0;

    % Compute the default adjusted SROI and number of blocks available
    [hold.image_size.rows,hold.image_size.cols, fps]  = dll_video('size', fn_orig); 
    hold.cvr = dll_calib_video('pvr');
    [sroi,vert,horiz] = adjust_requested_sroi (hold, 'vsize',vsize,'hsize',hsize, 'evenodd');
    dll_calib_video('sroi', sroi, 0); % 0 extra pixels needed
    tslice_length_sec = 1.0 / fps;
    number_tslices = min( dll_video('total_frames',fn_orig), dll_video('total_frames',fn_proc) );

    % Check if this is an interlace or progressive sequence
    % Allocate memory to hold all features for this clip.
    if strcmp(dll_video('get_video_standard', 1),'progressive'),
        is_progressive = 1;
        src = zeros(vert*horiz, number_tslices);
    elseif strcmp(dll_video('get_video_standard', 1),'interlace_lower_field_first') | strcmp(dll_video('get_video_standard', 1),'interlace_upper_field_first'),
        is_progressive = 0;
        src = zeros(vert*horiz, number_tslices*2);
    else
        error('video standard not recognized');
    end

    % read in each original field (or frame), and sub-sample it.
    dll_video('set_rewind', fn_orig);
    dll_video('set_tslice', fn_orig, tslice_length_sec);
    num = 0;
    for frame=1:number_tslices,
        y = dll_calib_video('tslice', fn_orig);
        if is_progressive,
            src_one(:,frame) = reshape( block_statistic(y,hsize,vsize,'mean'), vert*horiz, 1);
        else
            [one,two] = split_into_fields(y);

            src_one(:,frame) = reshape( block_statistic(one,hsize/2,vsize,'mean'), vert*horiz, 1);
            src_two(:,frame) = reshape( block_statistic(two,hsize/2,vsize,'mean'), vert*horiz, 1);
        end
    end
    dll_video('rewind', fn_orig);

    % normalize original images.  Don't amplify noise.
    for frame=1:number_tslices,
        src_one(:,frame) = src_one(:,frame) ./ max( single(1.0), std(src_one(:,frame)));
        if ~is_progressive,
            src_two(:,frame) = src_two(:,frame) ./ max( single(1.0), std(src_two(:,frame)));
        end
    end

    % read in processed fields (or frames), at the requested interval, and sub-sample them.
    all_align = [];

    uncert = round(fps * uncert_sec);
    all_std_of_diff = zeros(1,uncert*2+1);

    dll_video('set_rewind', fn_proc);
    dll_video('discard', fn_proc, tslice_length_sec * uncert);  
    dll_video('set_tslice', fn_proc, tslice_length_sec);
    for frame=1+uncert:number_tslices-uncert,
        % stop if run out of source images.
        if frame + uncert > number_tslices,
            break;
        end

        % read in next frame
        y = dll_calib_video('tslice', fn_proc);

        % split into fields if needed.  Record number /offset relative
        % to source array of images.
        if is_progressive,
            one = reshape( block_statistic(y,hsize,vsize,'mean'), vert*horiz, 1);
        else
            [one,two] = split_into_fields(y);
            one = reshape( block_statistic(one,hsize/2,vsize,'mean'), vert*horiz, 1);
            two = reshape( block_statistic(two,hsize/2,vsize,'mean'), vert*horiz, 1);
        end


        % normalize & find temporal registration for this image.
        one =one ./ max(single(1.0), std(one));
        [still, src_num, std_of_diff] = temporal_registration_search(src_one, one, frame, uncert, STILL_THRESHOLD);
        if ~still,
            all_align = [ all_align (frame-src_num) ];
            all_std_of_diff = all_std_of_diff + std_of_diff;
        end

        if ~is_progressive,
            % normalize & find temporal registration for this image.
            two = two ./ max(single(1.0), std(two));
            [still, src_num, std_of_diff] = temporal_registration_search(src_two, two, frame, uncert, STILL_THRESHOLD);
            if ~still,
                all_align = [ all_align (frame-src_num) ];
                all_std_of_diff = all_std_of_diff + std_of_diff;
            end
            % fprintf('\tframe %d field1 offset %d; field2 offset %d\n', frame, all_align(length(all_align)-1), frame-src_num);
        end
    end
    dll_video('rewind', fn_proc);
    
    % still check.  If still scene, record default alignment & go to
    % next processed clip.  Take mean of all_std_of_diff AFTER the
    % min/max, so that you only need to divide once (not 2*uncert+1 times).
    if length(all_align) == 0 | (max(all_std_of_diff) - min(all_std_of_diff)) / length(all_align) < STILL_THRESHOLD,
        delay = 0;
        is_still = is_still + 1;
        return;
    end

    % histogram, 
    hist_unsmooth = histc(all_align,-uncert:uncert);

    % smooth
    k=0:(2*HFW);
    fk = 0.5 + 0.5 * cos(pi * (k - HFW) / (1+HFW));
    fk = fk ./ sum(fk);
    hist_smooth = conv(fk,hist_unsmooth);

    % take off extra length, so left only with smoothed portion,
    % length is now (uncert-HFW)*2+1
    hist_smooth = hist_smooth((HFW*2+1):(length(hist_smooth) - HFW*2));

    %
    max_H_value = max(hist_unsmooth);
    [max_SH_value,max_SH_offset] = max(hist_smooth);
    delay = (max_SH_offset + HFW) - uncert - 1;

    % check if uncert was large enough.
    bins_missing = [1:HFW (length(hist_unsmooth) - HFW):length(hist_unsmooth)];
    if max(hist_unsmooth(bins_missing)) > max_H_value * BELOW_WARN,
        is_uncert = is_uncert + 1;
    end
    
    bins_away = [1:(max_SH_offset-DELTA-1) (max_SH_offset + DELTA+1):length(hist_smooth)];
    if hist_smooth(bins_away) > max_SH_value * BELOW_WARN,
        is_ambiguous = is_ambiguous + 1;
    end

    % Adjust read point (current frame) in video files.
    if delay > 0,
        dll_video('discard', fn_proc, delay * tslice_length_sec);
    elseif delay < 0,
        dll_video('discard', fn_orig, -delay * tslice_length_sec);
    end

    % clear variables
    clear src_one src_two extra_start extra_stop;
    clear one two clip curr_offsets loop too_few_start too_few_stop;
    clear clip_delay number_tslices number_tslices std_of_diff all_std_of_diff;


    sucess = 1;

% catch
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [still, src_num, std_of_diff] = temporal_registration_search(src, deg, num, uncert, still_threshold);
%

% compare processed & original images.
std_of_diff = std(src(:,(num-uncert):(num+uncert)) - repmat(deg,1,2*uncert+1));

% find minimum stdev of difference.
[a,where] = min(std_of_diff);

% detect whether still
if max(std_of_diff) - min(std_of_diff) < still_threshold,
    still = 1;
else
    still = 0;
end

src_num = (where - uncert - 1) + num;